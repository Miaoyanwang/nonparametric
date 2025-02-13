#SMMK
library(pracma)
library(rTensor)
library(quadprog)
library(clusterSim)


#################### Regular Lasso ####################
Lasso=function(xvec,y,xnew,lambda){
    fit <- glmnet(xvec,y,alpha = 1, family = "binomial",type.measure = "mse",lambda=0.1)
    coef = coef(fit)
    length(which(abs(coef)>0))
    lassoB = as.vector(coef)[-1]
    prob = predict(fit,newx = xnew,s = lambda,type = "response")
    res=list();
    res$B_est=lassoB; res$prob=prob
    return(res)
}

## visualization of the active region
visualization=function(B_est,percent=0.8){
d=dim(B_est)[1];H=dim(B_est)[3]
ave=matrix(0,nrow=d,ncol=d)
for(h in 1:(H-1)){
    m=B_est[,,h]
    m[abs(m)<=quantile(abs(m),percent)]=0
    ave=ave+m
}
ave[abs(ave)<=quantile(abs(ave),percent)]=0 ## display
return(ave)
}

#################### Our method ####################
NonparaMatrix=function(X,y,X_new=NULL,truer,support,H=11,option="unsym"){
    pred=NULL;d=dim(X[[1]])[1];
    B_est=array(0,dim=c(d,d,H-1));
    sparse=d-support;
    for(h in 1:(H-1)){
        c = ADMM(X,y,r = truer,srow = sparse,scol =sparse,lambda=0,rho.ini=0.1,p=1/H*h,option=option)
        if(is.null(X_new)) pred=rbind(pred,sign(c$fitted))
        else pred=rbind(pred,unlist(lapply(X_new,function(x)sign(c$predict(x)))))
        B_est[,,h]=c$B
    }
    res=list();
    res$B_est=B_est; res$prob=prob_est(pred);
    return(res)
}

#################### CNN ####################
CNN=function(X,y,X_new,plot.figure=FALSE){
    n=length(X); n1=length(X_new);d=dim(X[[1]])[1];y_recode=(y+1)/2;
    model <- keras_model_sequential() %>%
    layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu",
    input_shape = c(d,d,1)) %>%
    layer_max_pooling_2d(pool_size = c(5,5)) %>%
    
    layer_flatten() %>%
    layer_dense(units = 100, activation = "relu") %>%
    layer_dense(units = 1, activation = "sigmoid")
    
    model %>% compile(
    optimizer = "adam",
    loss = "binary_crossentropy",
    metrics = "accuracy"
    )
    
    xvec=xvec_new=NULL
    for(i in 1:n){
        xvec=rbind(xvec,c(X[[i]]))
    }
    for(i in 1:n1){
        xvec_new=rbind(xvec_new,c(X_new[[i]]))
    }
    mean=apply(xvec,2,mean)
    sd=apply(xvec,2,sd)
    sd[sd==0]=1
    xvec=sweep(xvec, 2, mean, FUN="-")
    xvec=sweep(xvec,2,sd,FUN="/")
    xvec_new=sweep(xvec_new, 2, mean, FUN="-")
    xvec_new=sweep(xvec_new,2,sd,FUN="/")
    
    x_array= array(xvec,dim=c(n,d,d,1))
    x_array_new= array(xvec_new,dim=c(n1,d,d,1))
 
    hold=sample(1:n,round(n/5),replace=F)
    y_test=y_recode[hold];x_test=array(x_array[hold,,,],dim=c(length(hold),d,d,1))
    y_train=y_recode[-hold];x_train=array(x_array[-hold,,,],dim=c(n-length(hold),d,d,1))
    history <- model %>%
    fit(
    x = x_train, y =y_train,
    epochs = 10,
    validation_data = unname(list(x_test,y_test)),
    )
    if(plot.figure==TRUE) plot(history)
    result=list();
    result$prob=predict_proba(model,x_array_new);
    result$class=2*predict_classes(model,x_array_new)-1;
    result$history=history;
    result$acc=output=rev(history$metric$val_accuracy)[1]
    return(result)
}

link=function(type=c("gap","linear","logistic","sqrt","others"),beta=1){
    if(type=="gap") return(function(x) 0.8*(x>=0)+0.2*(x<0))
    if(type=="logistic") return(function(x) sigmoid(x))
    if(type=="linear") return(function(x) pmin(pmax(0,x/2+1/2),1))
    if(type=="sqrt") return(function(x)  pmin(pmax(0,sign(x)*x^2/2+1/2),1))
    else return(function(x) pmin(pmax(0,sign(x)*abs(x)^beta/2+1/2),1))
}


Makepositive = function(mat){
  h = eigen(mat,symmetric = T)
  nmat = (h$vectors)%*%diag(pmax(h$values,10^-4),nrow=nrow(mat))%*%t(h$vectors)
  return(nmat)
}


Karray = function(X,kernel = function(X1,X2) t(X1)%*%X2, type="col"){
  N = length(X);
  if(type=="row"){
    Xt=list();
    for(i in 1:N){
      Xt[[i]]=t(X[[i]])
    }
    X=Xt
  }
  
  m= nrow(X[[1]]); n = ncol(X[[1]]);
  K = array(dim = c(N,N,n,n))
  for(i in 1:N){
    for(j in 1:N){
      K[i,j, , ] = kernel(X[[i]],X[[j]])
    }
  }
  return(as.tensor(K))
}


expkernel = function(Y,Z){
  n = ncol(Y)
  A = matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      A[i,j] = t(Y[,i]-Z[,j])%*%(Y[,i]-Z[,j])
    }
  }
  return(exp(-A))
}

polykernel = function(Y,Z,deg = 3){
  n = ncol(Y)
  return((t(Y)%*%Z+matrix(1,n,n))^deg)
}


linearkernel = function(X1,X2) t(X1)%*%X2

constkernel = function(X1,X2) return(matrix(0,nrow=ncol(X1),ncol=ncol(X1)))

SMMK_con = function(X,y,r,kernel_row = c("linear","poly","exp","const"),kernel_col = c("linear","poly","exp","const"), cost = 10, rep = 1, p = .5,sparse=0){
  result = list()
  
  # Default is linear kernel.
  kernel_row <- match.arg(kernel_row)
  if (kernel_row == "linear") {
    kernel_row = linearkernel
  }else if(kernel_row == "poly"){
    kernel_row = polykernel
  }else if(kernel_row =="exp"){
    kernel_row = expkernel
  }else if(kernel_row == "const"){
    kernel_row = constkernel
  }
  
  kernel_col <- match.arg(kernel_col)
  if (kernel_col == "linear") {
    kernel_col = linearkernel
  }else if(kernel_col == "poly"){
    kernel_col = polykernel
  }else if(kernel_col =="exp"){
    kernel_col = expkernel
  }else if(kernel_col == "const"){
    kernel_col = constkernel
  }
  
  d1 = nrow(X[[1]]); d2 = ncol(X[[1]]); n = length(X)
  K_row = Karray(X,kernel_row,type="row")
  K_col = Karray(X,kernel_col,type="col")
  compareobj = 10^10
  
  
  for(nsim in 1:rep){
    error = 10; iter = 0; obj = 10^10
    # initialize P_row,P_col
    P_row = randortho(d1)[,1:r,drop = F]; P_col = randortho(d2)[,1:r,drop = F]  
    
    
    
    while((iter < 20)&(error >10^-3)){
      # update C
      W_row = P_row%*%t(P_row); W_col = P_col%*%t(P_col)
      Dmat=matrix(unfold(K_row,c(1,2),c(3,4))@data%*%c(W_row)+unfold(K_col,c(1,2),c(3,4))@data%*%c(W_col),nrow=n,ncol=n)
      
      dvec = rep(1,n)
      Dmat = Makepositive((y%*%t(y))*Dmat)
      Amat = cbind(y,diag(1,n),-diag(1,n))
      bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
      res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      alpha=res$solution
      
      CPh_row=ttl(K_row,list(t(as.matrix(y*alpha)),t(P_row)),ms=c(1,3))
      CPh_col=ttl(K_col,list(t(as.matrix(y*alpha)),t(P_col)),ms=c(1,3))
      
      CC_row=ttl(CPh_row,list(t(as.matrix(y*alpha)),t(P_row)),ms=c(2,4))
      CC_col=ttl(CPh_col,list(t(as.matrix(y*alpha)),t(P_col)),ms=c(2,4))
      
      CC_row=as.matrix(CC_row@data[1,1,,])
      CC_col=as.matrix(CC_col@data[1,1,,])
      
      factor_row=unfold(ttm(CPh_row,sqrtm(Makepositive(CC_row))$Binv,3),2,c(1,3,4))@data
      factor_col=unfold(ttm(CPh_col,sqrtm(Makepositive(CC_col))$Binv,3),2,c(1,3,4))@data
      Dmat=factor_row%*%t(factor_row)+factor_col%*%t(factor_col)
      
      
      dvec = rep(1,n)
      Dmat = Makepositive((y%*%t(y))*Dmat)
      Amat = cbind(y,diag(1,n),-diag(1,n))
      bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
      res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      alpha=res$solution
      obj=c(obj,-res$value)
      iter = iter+1
      error = abs(-obj[iter+1]+obj[iter])/obj[iter]
      
      # P formula
      P_row=ttm(CPh_row,t(as.matrix(alpha*y)),2)[1,1,,]@data
      P_col=ttm(CPh_col,t(as.matrix(alpha*y)),2)[1,1,,]@data
      
      P_row=matrix(P_row,nrow=r)
      P_col=matrix(P_col,nrow=r)
      
      P_row = svd(P_row)$v
      P_col = svd(P_col)$v
      
      #### sparse model
      
      B=0
      for(i in 1:n){
        B=B+alpha[i]*y[i]*P_row%*%t(P_row)%*%X[[i]]+alpha[i]*y[i]*X[[i]]%*%P_col%*%t(P_col)
      }
      if(sparse>=1){
        B=sparse_matrix(B,r,sparse,sparse)
        P_row=svd(B)$u[,1:r]
        P_col=svd(B)$v[,1:r]
      }
      
    }
    if(compareobj>obj[iter+1]){
      P_row_optimum=P_row; P_col_optimum=P_col;
      obj_optimum=obj;
      compareobj=obj[iter+1]
    }
  }
  
  
  
  P_row= P_row_optimum; P_col= P_col_optimum;
  W_row = P_row%*%t(P_row); W_col = P_col%*%t(P_col);
  
  Dmat=K=matrix(unfold(K_row,c(1,2),c(3,4))@data%*%c(W_row)+unfold(K_col,c(1,2),c(3,4))@data%*%c(W_col),nrow=n,ncol=n)
  
  dvec = rep(1,n)
  Dmat = Makepositive((y%*%t(y))*Dmat)
  Amat = cbind(y,diag(1,n),-diag(1,n))
  bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
  res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
  alpha=res$solution
  
  slope = function(Xnew){
    newK = rep(0,n)
    for( i in 1:n){
      newK[i] = sum(W_row*kernel_row(t(Xnew),t(X[[i]])))+
        sum(W_col*kernel_col(Xnew,X[[i]]))
    }
    
    return(sum(alpha*y*newK))
  }
  
  # intercept part estimation (update b)
  yfit=K%*%(alpha*y) ## faster than lapply
  
  B=0;
  for(i in 1:n){
    B=B+alpha[i]*y[i]*P_row%*%t(P_row)%*%X[[i]]+alpha[i]*y[i]*X[[i]]%*%P_col%*%t(P_col)
  }
  if(sparse>=1){
    B=sparse_matrix(B,r,sparse,sparse);
    slope = function(Xnew) sum(B*Xnew)
    yfit = unlist(lapply(X,slope))
  }
  
  positive=min(yfit[y==1])
  negative=max(yfit[y==-1])
  if ((1-positive)<(-1-negative)) {
    intercept = -(positive+negative)/2
  }else{
    gridb0 = seq(from = -1-negative,to = 1-positive,length = 100)
    intercept = gridb0[which.min(sapply(gridb0,function(b) objective(b,yfit,y,p = p)))]
  }
  compareobj = obj[iter+1]
  predictor = function(Xnew) sign(slope(Xnew)+intercept)
  
  result$alpha = alpha
  result$slope = slope; result$predict = predictor
  result$intercept = intercept;
  result$P_row = P_row; result$P_col = P_col;
  result$obj = obj[-1]; result$iter = iter;
  result$error = error;
  result$fitted=yfit+intercept; ## add fitted value as a criterium to select cost
  result$B=B
  return(result)
}

hinge = function(y) ifelse(1-y>0,1-y,0)

objective=function(b,yfit,y,p){
  return(sum(hinge(y*(yfit+b))[y==1]*(1-p))+sum(hinge(y*(yfit+b))[y==-1]*p))
}
gradient=function(W,x,yfit,y,p){
  yfit=W%*%x+yfit;
  gradient1=matrix(-W[which((yfit<1)&(y==1)),],ncol=length(x))
  gradient2=matrix(W[which(((-yfit)<1)&(y==-1)),],ncol=length(x))
  
  return(apply(gradient1,2,sum)*(1-p)+apply(gradient2,2,sum)*p)
}

sparse_matrix=function(M,r,srow,scol,option="sym"){
  srow_target=srow;scol_target=scol;
  ## scheme 1: first low rank, then sparse
  #res=svd(M);step=0;
  #M_low=res$u[,1:r]%*%diag(res$d[1:r],nrow=r)%*%t(res$v[,1:r])
  #while((srow>0)&(scol>0)){
      #row_norm=diag(M_low%*%t(M_low))
    #M_low[sort(row_norm,index=T)$ix[step+1],]=0
    #srow=srow-1
    #col_norm=diag(t(M_low)%*%M_low)
    #M_low[,sort(col_norm,index=T)$ix[step+1]]=0
    #scol=scol-1; step=step+1
    #}
  #M_low_output=M_low; value1=sum((M-M_low)^2)
  
  ## scheme 2: first sparse, then low rank
  srow_index=scol_index=NULL
  M_low=M;srow=srow_target;scol=scol_target;step=0;
  while((srow>0)&(scol>0)){
    row_norm=diag(M_low%*%t(M_low))
    M_low[sort(row_norm,index=T)$ix[step+1],]=0
    srow_index=c(srow_index,sort(row_norm,index=T)$ix[step+1])
    srow=srow-1
    col_norm=diag(t(M_low)%*%M_low)
    if(option=="sym"){
      M_low[,sort(row_norm,index=T)$ix[step+1]]=0
    }else M_low[,sort(col_norm,index=T)$ix[step+1]]=0
    scol_index=c(scol_index,sort(col_norm,index=T)$ix[step+1])
    scol=scol-1; step=step+1;
  }
  res=svd(M_low) ## some zero entries become nonzero due to numerical precision
  M_low=res$u[,1:r]%*%diag(res$d[1:r],nrow=r)%*%t(res$v[,1:r])
  M_low[srow_index,]=0 ## reset sparse entries to zero
  M_low[,scol_index]=0 #
  if(option=="sym") diag(M_low)=0
  value2=sum((M-M_low)^2)
  # if(value2<=value1) M_low_output=M_low;
  M_low_output=M_low;
  return(M_low_output)
}

### alternative method for simutanously low-rank+sparse model

ADMM=function(X,y,Covariate=NULL,r,srow,scol,rho.ini=0.1,p=0.5,lambda=0,option="unsym"){
  result=list();
  n=length(X);
  average_scale=mean(unlist(lapply(X,function(x)norm(x,"F"))))
  rho=rho.ini/sqrt(n)*average_scale^2;  ## increase with X scale
  ##rho=rho.ini/sqrt(n)
  
  Lambda=matrix(0,nrow=nrow(X[[1]]),ncol=ncol(X[[1]]))
  obj=residual=NULL
  PQ=0; iter=0; obj = 10^10;error=10;residual=10^10
  rho_list=NULL;
  
  while(((iter < 10)|(error > 10^-4))){

    PQ_prev=PQ
    ### update B abd c
    res=SVM_offset(X,y,Covariate,OffsetC=(2*rho*PQ-Lambda)/(2*(lambda+rho)),p=p,cost=1/(2*(rho+lambda)))
    obj=c(obj,res$hinge) ## minimize objective
    B=res$coef
    
    ## Update PQ
    PQ=sparse_matrix(B+1/(2*rho)*Lambda,r,srow,scol,option)
    
    residual=c(residual,norm(B-PQ,"F")) ## primal residual
    #residual_dual=c(residual_dual,rho*norm(PQ-PQ_prev,"F")) ## dual residual
    
    ## update Lambda
    Lambda=Lambda+2*rho*(B-PQ)
    
    ## geometric step size
    if(iter>=10) rho=rho*1.1

    rho_list=c(rho_list,rho)
    gc();
    
    
    iter=iter+1;

    ##error=abs(-residual[iter+1]+residual[iter])
    error=residual[iter+1]/(r*(nrow(X[[1]])+ncol(X[[1]])-srow-scol));
    if(iter>=200) break
  }
  ## take PQ (exact low-rank + sparse) as the output
  slope = function(Xnew) sum(Xnew*PQ)
  
  if(length(Covariate)>=1){
    W=cbind(1,matrix(unlist(Covariate),nrow=n,byrow=TRUE))
  } else {
    W=as.matrix(rep(1,n))}
  
   yfit = unlist(lapply(X,function(x) sum(PQ*x)))
   b=rep(0,dim(W)[2])
   intercept_opt=optim(b,function(x)objective(W%*%x,yfit,y,p=p),function(x) gradient(W,x,yfit,y,p=p),method="BFGS")
   intercept=intercept_opt$par
   yfit=yfit+intercept
  
  
  predictor = function(Xnew,Covariate=NULL){
    if(length(Covariate)>=1){
      W=cbind(1,matrix(unlist(Covariate),nrow=1,byrow=TRUE))
      sign(sum(Xnew*PQ)+intercept)
    }
    else sign(sum(Xnew*PQ)+intercept)
  }
  
  result$slope=slope;result$predict=predictor;
  
  result$intercept=intercept;
  
  result$P_row=svd(PQ)$u[,1:r]; result$P_col=svd(PQ)$v[,1:r];
  result$obj=obj[-1];
  result$obj_exact=intercept_opt$value
  result$iter=iter;
  result$error=error;
  result$fitted=yfit;
  result$B=PQ;
  result$residual=residual[-1];result$B_appr=res$coef;result$rho=rho_list;
  
  return(result)
}

SVM_offset=function(X,y,Covariate,OffsetC,p=0.5,cost=1){
  offset=unlist(lapply(X,function(x)sum(x*OffsetC)))
  n=length(X)
  Dmat=K=matrix(unlist(X),nrow=n,byrow=TRUE)%*%t(matrix(unlist(X),nrow=n,byrow=TRUE))
  dvec = 1-y*offset
  Dmat = Makepositive((y%*%t(y))*Dmat)
  Amat = cbind(y,diag(1,n),-diag(1,n))
  
  bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
  res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
  
  ## calculate coefficient
  B=matrix((y*res$solution)%*%matrix(unlist(X),nrow=n,byrow=TRUE),nrow=nrow(X[[1]]))
  coef=B+OffsetC
  
  ### calculate hinge loss
  yfit=K%*%(res$solution*y)+offset
  #positive=min(yfit[y==1])
  #negative=max(yfit[y==-1])
  #if ((1-positive)<(-1-negative)) {
  #    intercept = -(positive+negative)/2
  #}else{
  #        gridb0 = seq(from = -1-negative,to = 1-positive,length = 100)
  #    intercept = gridb0[which.min(sapply(gridb0,function(b) objective(b,yfit,y,p = p)))]
  #    }
  
  if(length(Covariate[[1]])>=1){
    W=cbind(1,matrix(unlist(Covariate),nrow=n,byrow=TRUE))
  } else {
    W=as.matrix(rep(1,n))}
  b=rep(0,dim(W)[2])
  intercept=optim(b,function(x)objective(W%*%x,yfit,y,p=p),function(x) gradient(W,x,yfit,y,p=p),method="BFGS")$par
  
  yfit=yfit+W%*%intercept
  return(list("res"=res,"coef"=coef,"fitted"=yfit,"hinge"=objective(0,yfit,y,p=p),"intercept"=intercept))
}

### estimate probability from sequence of classifications.
## Input: cum a H-by-n matrix. each row is the classification result for h = (1,...H)/(H+1)
prob_est=function(cum,option=3){
    H=dim(cum)[1]
    if(option==1){
        cum_sum=apply(cum,2,cumsum)
        prob=apply(cum_sum,2,which.max)/(H+1)
        return(prob) ## use maximum cumulative probability
    }else if(option==2){
        d=dim(cum)[2]
        prob=rep(0,d)
        for(i in 1:d){
            vector=cum[,i]
            vector=c(1,vector,-1)
            index1=max(which(vector==1))-1
            index2=max(which(rev(vector)==-1))-1
            prob[i]=(index1+H-index2)/(2*(H+1)) ## use proposal in Biometrika (2008) Wang, Shen & Liu.
        }
        return(prob)
    }else if(option==3){
        prob=apply(cum,2,sum)/(2*(H+1))+1/2
        return(prob)
    }
}



################################### Option 1 ###################################

SMM = function(X,y,r,kernel_row = c("linear","poly","exp","const"),kernel_col = c("linear","poly","exp","const"),cost = 10, rep = 1, p = .5){
  result = list()
  
  # Default is linear kernel.
  kernel_row <- match.arg(kernel_row)
  if (kernel_row == "linear") {
    kernel_row = linearkernel
  }else if(kernel_row == "poly"){
    kernel_row = polykernel
  }else if(kernel_row =="exp"){
    kernel_row = expkernel
  }else if(kernel_row == "const"){
    kernel_row = constkernel
  }
  
  kernel_col <- match.arg(kernel_col)
  if (kernel_col == "linear") {
    kernel_col = linearkernel
  }else if(kernel_col == "poly"){
    kernel_col = polykernel
  }else if(kernel_col =="exp"){
    kernel_col = expkernel
  }else if(kernel_col == "const"){
    kernel_col = constkernel
  }
  
  d1 = nrow(X[[1]]); d2 = ncol(X[[1]]); n = length(X)
  K = Karray(X,kernel_row,type="row")
  #K_col = Karray(X,kernel_col,type="col")
  compareobj = 10^10
  
  
  for(nsim in 1:rep){
    error = 10; iter = 0; obj = 10^10
    # initialize P_row,P_col
    P =  randortho(d1)[,1:r,drop = F]
    
    
    
    while((iter < 20)&(error >10^-3)){
      # update C
      W = P%*%t(P);# W_col = P_col%*%t(P_col)
      Dmat=matrix(unfold(K,c(1,2),c(3,4))@data%*%c(W),nrow=n,ncol=n)
      
      dvec = rep(1,n)
      Dmat = Makepositive((y%*%t(y))*Dmat)
      Amat = cbind(y,diag(1,n),-diag(1,n))
      bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
      res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      alpha=res$solution
      
      CPh=ttl(K,list(t(as.matrix(y*alpha)),t(P)),ms=c(1,3))
      #CPh_col=ttl(K_col,list(t(as.matrix(y*alpha)),t(P_col)),ms=c(1,3))
      
      CC=ttl(CPh,list(t(as.matrix(y*alpha)),t(P)),ms=c(2,4))
      #CC_col=ttl(CPh_col,list(t(as.matrix(y*alpha)),t(P_col)),ms=c(2,4))
      
      CC=as.matrix(CC@data[1,1,,])
      #CC_col=as.matrix(CC_col@data[1,1,,])
      
      factors=unfold(ttm(CPh,sqrtm(Makepositive(CC))$Binv,3),2,c(1,3,4))@data
      #factor_col=unfold(ttm(CPh_col,sqrtm(Makepositive(CC_col))$Binv,3),2,c(1,3,4))@data
      Dmat=factors%*%t(factors)#+factor_col%*%t(factor_col)
      
      
      dvec = rep(1,n)
      Dmat = Makepositive((y%*%t(y))*Dmat)
      Amat = cbind(y,diag(1,n),-diag(1,n))
      bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
      res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      alpha=res$solution
      obj=c(obj,-res$value)
      iter = iter+1
      error = abs(-obj[iter+1]+obj[iter])/obj[iter]
      
      # P formula
      P=ttm(CPh,t(as.matrix(alpha*y)),2)[1,1,,]@data
      #P_col=ttm(CPh_col,t(as.matrix(alpha*y)),2)[1,1,,]@data
      
      P=matrix(P,nrow=r)
      #P_col=matrix(P_col,nrow=r)
      
      P = svd(P)$v
      #P_col = svd(P_col)$v
      
      #### sparse model
      
      B=0
      A = 0
      for(i in 1:n){
        B=B+alpha[i]*y[i]*P%*%t(P)%*%X[[i]]#+alpha[i]*y[i]*X[[i]]%*%P_col%*%t(P_col)
        A = A+alpha[i]*y[i]*X[[i]]
      }
      
      
    }
    if(compareobj>obj[iter+1]){
      P_optimum=P; #P_col_optimum=P_col;
      obj_optimum=obj;
      compareobj=obj[iter+1]
    }
  }
  
  
  
  P= P_optimum;# P_col= P_col_optimum;
  W = P%*%t(P);# W_col = P_col%*%t(P_col);
  
  Dmat= K = matrix(unfold(K,c(1,2),c(3,4))@data%*%c(W),nrow=n,ncol=n)
  
  dvec = rep(1,n)
  Dmat = Makepositive((y%*%t(y))*Dmat)
  Amat = cbind(y,diag(1,n),-diag(1,n))
  bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
  res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
  alpha=res$solution
  
  slope = function(Xnew){
    newK = rep(0,n)
    for( i in 1:n){
      newK[i] = sum(W*kernel_row(t(Xnew),t(X[[i]])))
    }
    
    return(sum(alpha*y*newK))
  }
  
  # intercept part estimation (update b)
  yfit=K%*%(alpha*y) ## faster than lapply
  
  B=0;
  for(i in 1:n){
    B=B+alpha[i]*y[i]*P%*%t(P)%*%X[[i]]#+alpha[i]*y[i]*X[[i]]%*%P_col%*%t(P_col)
  }
  
  
  positive=min(yfit[y==1])
  negative=max(yfit[y==-1])
  if ((1-positive)<(-1-negative)) {
    intercept = -(positive+negative)/2
  }else{
    gridb0 = seq(from = -1-negative,to = 1-positive,length = 100)
    intercept = gridb0[which.min(sapply(gridb0,function(b) objective(b,yfit,y,p = p)))]
  }
  compareobj = obj[iter+1]
  predictor = function(Xnew) sign(slope(Xnew)+intercept)
  
  result$alpha = alpha
  result$slope = slope; result$predict = predictor
  result$intercept = intercept;
  result$P_row = P; result$P_col = P;
  result$obj = obj[-1]; result$iter = iter;
  result$error = error;
  result$fitted=yfit+intercept; ## add fitted value as a criterium to select cost
  result$B=B
  return(result)
}



################################### Option 1 ###################################

SMMK = function(X,y,r,kernel_row = c("linear","poly","exp","const"),kernel_col = c("linear","poly","exp","const"), cost = 10, rep = 1, p = .5){
  result = list()
  
  # Default is linear kernel.
  kernel_row <- match.arg(kernel_row)
  if (kernel_row == "linear") {
    kernel_row = linearkernel
  }else if(kernel_row == "poly"){
    kernel_row = polykernel
  }else if(kernel_row =="exp"){
    kernel_row = expkernel
  }else if(kernel_row == "const"){
    kernel_row = constkernel
  }
  
  kernel_col <- match.arg(kernel_col)
  if (kernel_col == "linear") {
    kernel_col = linearkernel
  }else if(kernel_col == "poly"){
    kernel_col = polykernel
  }else if(kernel_col =="exp"){
    kernel_col = expkernel
  }else if(kernel_col == "const"){
    kernel_col = constkernel
  }
  
  d1 = nrow(X[[1]]); d2 = ncol(X[[1]]); n = length(X)
  K_row = Karray(X,kernel_row,type="row")
  K_col = Karray(X,kernel_col,type="col")
  compareobj = 10^10
  
  # Choose non zero columns and rows
  
  
  
  for(nsim in 1:rep){
    error = 10; iter = 0; obj = 10^10
    # initialize P_row,P_col
    if (d1 == d2) {
      P_row <- P_col <- randortho(d1)[,1:r,drop = F]
    }else{
      P_row = randortho(d1)[,1:r,drop = F]; P_col = randortho(d2)[,1:r,drop = F]  
    }
    
    
    while((iter < 20)&(error >10^-3)){
      # update C
      W_row = P_row%*%t(P_row); W_col = P_col%*%t(P_col)
      Dmat=matrix(unfold(K_row,c(1,2),c(3,4))@data%*%c(W_row)+unfold(K_col,c(1,2),c(3,4))@data%*%c(W_col),nrow=n,ncol=n)
      
      dvec = rep(1,n)
      Dmat = Makepositive((y%*%t(y))*Dmat)
      Amat = cbind(y,diag(1,n),-diag(1,n))
      bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
      res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      alpha=res$solution
      
      CPh_row=ttl(K_row,list(t(as.matrix(y*alpha)),t(P_row)),ms=c(1,3))
      CPh_col=ttl(K_col,list(t(as.matrix(y*alpha)),t(P_col)),ms=c(1,3))
      
      CC_row=ttl(CPh_row,list(t(as.matrix(y*alpha)),t(P_row)),ms=c(2,4))
      CC_col=ttl(CPh_col,list(t(as.matrix(y*alpha)),t(P_col)),ms=c(2,4))
      
      CC_row=as.matrix(CC_row@data[1,1,,])
      CC_col=as.matrix(CC_col@data[1,1,,])
      
      factor_row=unfold(ttm(CPh_row,sqrtm(Makepositive(CC_row))$Binv,3),2,c(1,3,4))@data
      factor_col=unfold(ttm(CPh_col,sqrtm(Makepositive(CC_col))$Binv,3),2,c(1,3,4))@data
      Dmat=factor_row%*%t(factor_row)+factor_col%*%t(factor_col)
      
      
      dvec = rep(1,n)
      Dmat = Makepositive((y%*%t(y))*Dmat)
      Amat = cbind(y,diag(1,n),-diag(1,n))
      bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
      res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      alpha=res$solution
      obj=c(obj,-res$value)
      iter = iter+1
      error = abs(-obj[iter+1]+obj[iter])/obj[iter]
      
      # P formula
      P_row=ttm(CPh_row,t(as.matrix(alpha*y)),2)[1,1,,]@data
      P_col=ttm(CPh_col,t(as.matrix(alpha*y)),2)[1,1,,]@data
      
      P_row=matrix(P_row,nrow=r)
      P_col=matrix(P_col,nrow=r)
      
      P_row = svd(P_row)$v
      P_col = svd(P_col)$v
      
      
      
      B=0
      for(i in 1:n){
        B=B+alpha[i]*y[i]*P_row%*%t(P_row)%*%X[[i]]+alpha[i]*y[i]*X[[i]]%*%P_col%*%t(P_col)
      }
      
    }
    if(compareobj>obj[iter+1]){
      P_row_optimum=P_row; P_col_optimum=P_col;
      obj_optimum=obj;
      compareobj=obj[iter+1]
    }
  }
  
  
  
  P_row= P_row_optimum; P_col= P_col_optimum;
  W_row = P_row%*%t(P_row); W_col = P_col%*%t(P_col);
  
  Dmat=K=matrix(unfold(K_row,c(1,2),c(3,4))@data%*%c(W_row)+unfold(K_col,c(1,2),c(3,4))@data%*%c(W_col),nrow=n,ncol=n)
  
  dvec = rep(1,n)
  Dmat = Makepositive((y%*%t(y))*Dmat)
  Amat = cbind(y,diag(1,n),-diag(1,n))
  bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
  res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
  alpha=res$solution
  
  slope = function(Xnew){
    newK = rep(0,n)
    for( i in 1:n){
      newK[i] = sum(W_row*kernel_row(t(Xnew),t(X[[i]])))+
        sum(W_col*kernel_col(Xnew,X[[i]]))
    }
    
    return(sum(alpha*y*newK))
  }
  
  # intercept part estimation (update b)
  yfit=K%*%(alpha*y) ## faster than lapply
  
  B=0;
  for(i in 1:n){
    B=B+alpha[i]*y[i]*P_row%*%t(P_row)%*%X[[i]]+alpha[i]*y[i]*X[[i]]%*%P_col%*%t(P_col)
  }
  
  
  positive=min(yfit[y==1])
  negative=max(yfit[y==-1])
  if ((1-positive)<(-1-negative)) {
    intercept = -(positive+negative)/2
  }else{
    gridb0 = seq(from = -1-negative,to = 1-positive,length = 100)
    intercept = gridb0[which.min(sapply(gridb0,function(b) objective(b,yfit,y,p = p)))]
  }
  compareobj = obj[iter+1]
  predictor = function(Xnew) sign(slope(Xnew)+intercept)
  
  result$alpha = alpha
  result$slope = slope; result$predict = predictor
  result$intercept = intercept;
  result$P_row = P_row; result$P_col = P_col;
  result$obj = obj[-1]; result$iter = iter;
  result$error = error;
  result$fitted=yfit+intercept; ## add fitted value as a criterium to select cost
  result$B=B
  return(result)
}


################################### Combination 1 ###################################
SMMK_sparse = function(X,y,r,kernel_row = c("linear","poly","exp","const"),kernel_col = c("linear","poly","exp","const"),option = c("approximate","exact"), cost = 10, rep = 1, p = .5,sparse=0){
  result = list()
  
  option <- match.arg(option)
  d1 = nrow(X[[1]]); d2 = ncol(X[[1]]); n = length(X)
  if(sparse>0){
    
    res = SMMK(X,y,r,kernel_row,kernel_col,cost,rep,p)
    initB = res$B
    row_o = order(diag(initB%*%t(initB)),decreasing = T)[1:(d1-sparse)]
    col_o = order(diag(t(initB)%*%initB),decreasing = T)[1:(d1-sparse)]
  }else{
    row_o = 1:d1; col_o = 1:d2
  }
  
  
  
  X_sp = lapply(X,function(x) x[row_o,col_o,drop = F])
  d1sp = nrow(X_sp[[1]]); d2sp = ncol(X_sp[[1]]); n = length(X_sp)
  if(option == "exact"){
    res = SMM(X_sp,y,r,kernel_row,kernel_col,cost,rep,p)
  }else{
    res = SMMK(X_sp,y,r,kernel_row,kernel_col,cost,rep,p)
  }
  
  
  B = matrix(0,nrow = d1,ncol = d2)
  B[row_o,col_o] = res$B
  
  slope = function(Xnew) res$slope(Xnew[row_o,col_o])
  intercept = res$intercept
  
  predictor = function(Xnew) sign(slope(Xnew)+intercept)
  
  result$alpha = res$alpha
  result$slope = slope; result$predict = predictor
  result$intercept = intercept;
  result$P_row = res$P_row; result$P_col = res$P_col;
  result$obj = res$obj; result$iter = res$iter;
  result$error = res$error;
  result$fitted=res$fitted;
  result$B=B
  return(result)
}


convpb = function(temp,H){ 
  temp=c(1,temp,-1)
  index1=max(which(temp==1))-1
  index2=max(which(rev(temp)==-1))-1
  return((index1+H-index2)/(2*(H+1)))
}

