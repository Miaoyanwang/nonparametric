library(rTensor)
library(Matrix)
library(quadprog)
min.thresh=10^(-3)
max.thresh=10^10

binaryloss=function(Ybar,W,Yfit){
  return(mean(W*abs(Ybar-sign(Yfit)),na.rm=TRUE))
}

loss=function(y,type=c("logistic","hinge")){
  if(type=="hinge") return(ifelse(1-y>0,1-y,0))
  if(type=="logistic") return(log(1+exp(-y)))
}



likelihood = function(data,theta){
  index=which(is.na(data)==F & is.na(theta)==F)
  return(sqrt(sum((data[index]-theta[index])^2)))
}

gradientm=function(A1,A2,mode,Ybar,W,type=c("logistic","hinge")){
  margin = Ybar*(A1%*%t(A2))
  d = dim(Ybar)
  if(type=="logistic"){
    tem=-W*Ybar*exp(-margin)/(1+exp(-margin))
  }else if(type=="hinge"){
    tem=-W*Ybar*(margin<1)
  }
  scale = length(tem)-sum(is.na(tem))
  tem[is.na(tem)]=0
  
  if(mode==1){
    Grad = tem%*%A2/scale
  }else if(mode==2){
    Grad = t(tem)%*%A1/scale
  }
  return(Grad)
}       

costm = function(A1,A2,Ybar,W,type = c("logistic","hinge")){
  return(mean(W*loss((A1%*%t(A2))*Ybar,type),na.rm=TRUE))
}


SignM=function(Y,truer,H=5,Lmin,Lmax,option=1){
  B_fitted=result=list()
  pi_seq=seq(from=Lmin,to=Lmax,length=2*H+1)
  
  for(h in 2:(2*H)){
    pi=pi_seq[h]
    if(option==1){
      res=Altm(sign(Y-pi),abs(Y-pi),r=truer,type="logistic",start="linear")## recommend
    }else if(option==2){
      res=Altm(sign(Y-pi),abs(Y-pi),r=truer,type="hinge",start="linear")## recommend
    }
    B_fitted[[h]]=res$fitted
  }
  B_fitted=array(unlist(B_fitted),dim=c(dim(Y),2*H-1));
  res=list();
  res$fitted=B_fitted
  res$est=1/2*(apply(sign(B_fitted),1:2,sum)/(2*H)+1)*(Lmax-Lmin)+Lmin
  return(res)
}




Altm=function(Ybar,W,r,type=c("logistic","hinge"),start="linear"){
  result=list()
  d=dim(Ybar)
  if(start=="linear"){
    ini=fit_continuous(Ybar,r)
    A1 = ini$U[[1]];
    scale=matrix(0,nrow=r,ncol=r)
    diag(scale)=ini$lambda
    A2 = ini$U[[2]]%*%scale;
  }else{
    A1 = matrix(runif(d[1]*r,-1,1),nrow = d[1],ncol = r);
    A2 = matrix(runif(d[2]*r,-1,1),nrow = d[2],ncol = r);
  }
  obj=costm(A1,A2,Ybar,W,type)
  binary_obj=binaryloss(Ybar,W,A1%*%t(A2))
  
  error=1;iter=1;
  
  while((error>0.01)&(iter<20)){
    
    
    optimization=optim(c(A2),function(x) costm(A1,matrix(x,ncol=r),Ybar,W,type),function(x) c(gradientm(A1,matrix(x,ncol=r),2,Ybar,W,type)),method="BFGS")
    A2=matrix(optimization$par,ncol=r)
    
    optimization=optim(c(A1),function(x)costm(matrix(x,ncol=r),A2,Ybar,W,type),function(x)gradientm(matrix(x,ncol=r),A2,1,Ybar,W,type),method="BFGS")
    A1=matrix(optimization$par,ncol=r)
    
    obj=c(obj,costm(A1,A2,Ybar,W,type))
    binary_obj=c(binary_obj,binaryloss(Ybar,W,A1%*%t(A2)))
    iter=iter+1
    error=(obj[iter-1]-obj[iter])
    
  }
  result$binary_obj=binary_obj;
  result$obj=obj;
  result$iter=iter;
  result$error=error;
  result$fitted=A1%*%t(A2); ## exact low-rank
  return(result)
}



fit_continuous=function(data,r){
  index=which(is.na(data)==TRUE)
  data[index]=mean(data,na.rm=TRUE)
  original_data=data
  
  if(length(dim(data))>=3){
    sink(tempfile())
    decomp=tryCatch(cp(as.tensor(original_data),r),error=function(c)"degeneracy")
    sink()
    if(inherits(decomp,"character")==TRUE){
      U=list();
      U[[1]]=matrix(0,ncol=r,nrow=dim(data)[1])
      U[[2]]=matrix(0,ncol=r,nrow=dim(data)[2])
      U[[3]]=matrix(0,ncol=r,nrow=dim(data)[3])
      return(list("est"=array(NA,dim=dim(data)),"U"=U,"lambda"=rep(0,r),"info"="degeneracy"))
    }
    res0=1
    res=0
    thresh=10^(-3)
    error=NULL
    
    while((res0-res)>thresh){
      res0=likelihood(original_data,decomp$est@data)
      sink(tempfile())
      decomp=cp(as.tensor(data),r)
      sink()
      res=likelihood(original_data,decomp$est@data)
      data[index]=decomp$est@data[index]
      error=c(error,res)
    }
    
    return(list("est"=decomp$est@data,"U"=decomp$U,"lambda"=decomp$lambda))
  }else if(length(dim(data))==2){
    PQ=svd(data)
    if(r==1){
      est=cbind(PQ$u[,1:r])%*%diag(as.matrix(PQ$d[1:r]))%*%t(cbind(PQ$v[,1:r]))
      U = list(PQ$u[,1:r],PQ$v[,1:r])
      lambda = PQ$d[1:r]
    }else{
      est=PQ$u[,1:r]%*%diag(PQ$d[1:r])%*%t(PQ$v[,1:r])
      U = list(PQ$u[,1:r],PQ$v[,1:r])
      lambda = PQ$d[1:r]
    }
    res0=1
    res=0
    thresh=10^(-3)
    error=NULL
    
    while((res0-res)>thresh){
      res0=likelihood(original_data,est)
      if(r==1){
        est=cbind(PQ$u[,1:r])%*%diag(as.matrix(PQ$d[1:r]))%*%t(cbind(PQ$v[,1:r]))
        U = list(PQ$u[,1:r],PQ$v[,1:r])
        lambda = PQ$d[1:r]
      }else{
        est=PQ$u[,1:r]%*%diag(PQ$d[1:r])%*%t(PQ$v[,1:r])
        U = list(PQ$u[,1:r],PQ$v[,1:r])
        lambda = PQ$d[1:r]
      }
      res=likelihood(original_data,est)
      data[index]=est[index]
      error=c(error,res)
    }
    return(list("est" = est,"U" = U,"lambda"= lambda))
  }
}

