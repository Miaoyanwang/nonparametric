library(pracma)
library(rTensor)
library(quadprog)
library(Matrix)
min.thresh=10^(-3)
max.thresh=10^10
## test
set.seed(1)
d=50
a=randortho(d)[,1]
b=randortho(d)[,1]
c=randortho(d)[,1]

##### simulate graphon models
a=seq(from=0,to=1,length=d)
b=seq(from=0,to=1,length=d)
#c=seq(from=0,to=1,length=d)
signal=graphon_to_tensor(a,b,0,type=10)
signal=signal[,,1]
#plot(svd(signal)$d)
#signal=sigmoid(a%o%b,a=1)
#signal=sigmoid(a%o%b%o%c,a=50)
Y=signal+array(rnorm(length(signal),0,0*max(abs(signal))),dim=dim(signal))
truer=2
hist(Y)
missing=array(rbinom(length(signal),1,0),dim=dim(signal))
Y[missing==1]=NA

Lmin=min(Y,na.rm=T)
Lmax=max(Y,na.rm=T)

set.seed(1)
res=nonMAR(Y,truer,Lmin=min(signal),Lmax=max(signal),H=10,rho=.001)
plot(res$est,signal)
abline(0,1)
plot(res$est[missing==0],signal[missing==0])
plot(res$est[missing==1],signal[missing==1])
mean(abs(res$est[missing==1]-signal[missing==1]))
abline(0,1)


### matrix
est2=fit_continuous(Y,truer)
plot(est2,signal)
plot(est2[missing==0],signal[missing==0])
plot(est2[missing==1],signal[missing==1])
mean(abs(est2[missing==1]-signal[missing==1]))
abline(0,1)


#################### Our method ####################
graphon_to_tensor=function(a,b,c,type){
    d1=length(a);d2=length(b);d3=length(c)
    M=array(0,dim=c(d1,d2,d3))
    if(type==10){
        for(i in 1:d1){
            for(j in 1:d2){
                for(k in 1:d3){
                M[i,j,k]=log(1+0.5*max(a[i],b[j],c[k]))
                }
            }
        }
    }
    if(type==9){
        for(i in 1:d1){
            for(j in 1:d2){
                for(k in 1:d3){
                    M[i,j,k]=exp(-0.5*(min(a[i],b[j],c[k]))) ##
            }
            }
        }
    }
    if(type==6){
        for(i in 1:d1){
            for(j in 1:d2){
                for(k in 1:d3){
                    M[i,j,k]=abs(a[i]-b[j]) ## full rank
                }
            }
        }
    }
    if(type==7){
        for(i in 1:d1){
            for(j in 1:d2){
                for(k in 1:d3){
                    M[i,j,k]=1/(1+exp(-max(a[i],b[j],c[k])^2-min(a[i],b[j],c[k])^4))
                }
            }
        }
    }
    if(type==8){
        for(i in 1:d1){
            for(j in 1:d2){
                for(k in 1:d3){
                    M[i,j,k]=exp(-max(a[i],b[j],c[k])^(3/4))
                }
            }
        }
    }
    return(M)
}


####################  our method ####################
nonMAR=function(Y,truer,H=5,Lmin,Lmax,rho=0.1,lambda=10^(-3)){
    B_fitted=result=list();
    pi_seq=seq(from=Lmin,to=Lmax,length=2*H+1)
    for(h in 2:(2*H)){
        pi=pi_seq[h]
        res=ADMM(sign(Y-pi),abs(Y-pi),r = truer,rho=rho,lambda=lambda)
        result[[h]]=res
        B_fitted[[h]]=res$fitted
    }
    B_fitted=array(unlist(B_fitted),dim=c(dim(Y),2*H-1));
    res=list();
    res$result=result;
    res$fitted=B_fitted
    res$est=1/2*(apply(sign(B_fitted),1:length(dim(Y)),mean)+1+2/H)*(Lmax-Lmin)+Lmin
    return(res)
}

### ADMM for classification
ADMM=function(Ybar,W,r,rho=0.1,lambda=10^(-3)){
  result=list();
  
  Lambda=array(0,dim=dim(Ybar))
  PQ=0; iter=0; obj =residual=error=max.thresh;
  
  rho_list=NULL;
  while(((iter < 10)|(error > 10^-3))){

    PQ_prev=PQ
    
    ### update B
    if((rho+lambda)!=0){
    res=SVM_offset(Ybar,W,OffsetC=(2*rho*PQ-Lambda)/(2*(lambda+rho)),cost=1/(2*(rho+lambda)))
    }else if((rho+lambda)==0){
         res=SVM_offset(Ybar,W,OffsetC=array(0,dim=dim(W)),cost=max.thresh)
    }
    obj=c(obj,res$hinge) ## minimize objective
    B=res$coef
    
    ## Update PQ
    if(rho==0){PQ=B}
    else{
    if(length(dim(Ybar))>2){
    PQ=cp(as.tensor(B+1/(2*rho)*Lambda),r)
        PQ=PQ$est@data
    }else if(length(dim(Ybar))==2){
        PQ=svd(B+1/(2*rho)*Lambda)
        if(r==1){
            PQ=cbind(PQ$u[,1:r])%*%diag(as.matrix(PQ$d[1:r]))%*%t(cbind(PQ$v[,1:r]))
        }else{
            PQ=PQ$u[,1:r]%*%diag(PQ$d[1:r])%*%t(PQ$v[,1:r])
        }
    }
    }
    
    
    residual=c(residual,sqrt(sum((B-PQ)^2)))
    
    ## update Lambda
    Lambda=Lambda+2*rho*(B-PQ)
    
    ## geometric step size
    if(iter>=10){
        rho=rho*1.1;
        lambda=lambda*1.1;
    }
    

    rho_list=c(rho_list,rho)
    iter=iter+1;

    error=abs(-residual[iter+1]+residual[iter])
    if(iter>=200) break
  }
  
  result$obj=obj[-1];
  result$iter=iter;
  result$error=error;
  result$fitted=PQ; ## exact low-rank
  result$B=B; ## approximate low-rank from SVM
  result$residual=residual[-1];result$rho=rho_list;
  result$alpha=res$res$solution
  return(result)
}

SVM_offset=function(Ybar,W,OffsetC,cost=1){
  n=length(Ybar)
  missing=which(is.na(Ybar)==T)
  nonmissing=setdiff(1:n,missing)
  
  m=length(Ybar[nonmissing])
  dvec = 1-c(Ybar[nonmissing]*OffsetC[nonmissing])
  Dmat = diag(1,m)
  Amat = cbind(c(Ybar[nonmissing]),diag(1,m),-diag(1,m))
  bvec = c(rep(0,1+m),-c(cost*W[nonmissing]))
  res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
  
  ## calculate coefficient
  coef=OffsetC
  coef[nonmissing]=coef[nonmissing]+res$solution*Ybar[nonmissing]

  return(list("res"=res,"coef"=coef,"hinge"=objective(coef[nonmissing],Y[nonmissing],W[nonmissing])))
}

hinge = function(y) ifelse(1-y>0,1-y,0)

objective=function(yfit,Y,W){
    return(sum(hinge(Y*yfit)*W))
}

likelihood = function(data,theta){
    index=which(is.na(data)==F & is.na(theta)==F)
   return(sqrt(sum((data[index]-theta[index])^2)))
}

################################### normalize each column of X to be unit-1 ###################################
normalize=function(X){
    d=dim(X)[2]
    for(i in 1:d){
        X[,i]=X[,i]/sqrt(sum(X[,i]^2))
    }
    return(X)
}

##################### construct CP tensor using factor matrices X, Y, Z ###################################
tensorize=function(X,Y,Z){
    r=dim(X)[2]
    tensor=0
    if(is.matrix(X)==0){
        tensor=X%o%Y%o%Z
        return(tensor)
    }
    
    for(i in 1:r){
        tensor=tensor+X[,i]%o%Y[,i]%o%Z[,i]
    }
    return(tensor)
}

fit_continuous=function(data,r){
    original_data=data
    index=which(is.na(data)==T)
    data[index]=mean(data,na.rm=T)
    if(length(dim(data))>=3){
    decomp=cp(as.tensor(data),r)
    res0=1
    res=0
    thresh=10^(-3)
    error=NULL
    
    while((res0-res)>thresh){
    res0=likelihood(original_data,decomp$est@data)
    decomp=cp(as.tensor(data),r)
    res=likelihood(original_data,decomp$est@data)
    data[index]=decomp$est@data[index]
    error=c(error,res)
    }
    return(decomp$est@data)
    }else if(length(dim(data))==2){
        PQ=svd(data)
        if(r==1){
            decomp=cbind(PQ$u[,1:r])%*%diag(as.matrix(PQ$d[1:r]))%*%t(cbind(PQ$v[,1:r]))
        }else{
            decomp=PQ$u[,1:r]%*%diag(PQ$d[1:r])%*%t(PQ$v[,1:r])
        }
        res0=1
        res=0
        thresh=10^(-3)
        error=NULL
        
        while((res0-res)>thresh){
            res0=likelihood(original_data,decomp)
            if(r==1){
                decomp=cbind(PQ$u[,1:r])%*%diag(as.matrix(PQ$d[1:r]))%*%t(cbind(PQ$v[,1:r]))
            }else{
                decomp=PQ$u[,1:r]%*%diag(PQ$d[1:r])%*%t(PQ$v[,1:r])
            }
            res=likelihood(original_data,decomp)
            data[index]=decomp[index]
            error=c(error,res)
        }
        return(decomp)
}
}
