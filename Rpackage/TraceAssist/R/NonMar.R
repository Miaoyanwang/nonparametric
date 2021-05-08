

f=function(x,area,sign=1,rev=1){
  x=rev*x+(1-rev)*(1-x)
  if(area==1) return(sign*6*x^2)
  if(area==2) return(sign*3*sqrt(x))
  if(area==3) return(sign*(5-5*exp(-x)))
  if(area==4) return(sign*3*(exp(x)-1))
  if(area==5) return(sign*2*log(5*x+1))
  if(area==6) return(sign*4*sin(x))
  if(area==7) return(sign*3*tan(x))
  if(area==8) return(sign*8*x^3)
  if(area==9) return(sign*10*x^4)
  if(area==10) return(sign*5*x^(3/2))
}


############## functional matrix decomposition #########################

function_matrix_denoise=function(res,X,node_ID1,node_ID2,group){
  n=length(X)
  d=dim(X[[1]])[1]
  order=sort(res$prob,index=T,decreasing=F)$ix
  pb_sort=sort(res$prob,decreasing=F)
  new_X=array(0,dim=c(d,d,n))
  for(i in 1:n){
    new_X[,,i]=X[order][[i]]
  }
  if((length(node_ID1)==1)&(length(node_ID2)==1)){
    effect=new_X[node_ID1,node_ID2,]
    sd=0
  }else if((length(node_ID1)==1)&(length(node_ID2)>1)){
    effect=apply(new_X[node_ID1,node_ID2,],2,mean)
    sd=apply(new_X[node_ID1,node_ID2,],2,sd)
  }else if((length(node_ID1)>1)&(length(node_ID2)==1)){
    effect=apply(new_X[node_ID1,node_ID2,],1,mean)
    sd=apply(new_X[node_ID1,node_ID2,],1,sd)
  }else{
    effect=apply(new_X[node_ID1,node_ID2,],3,mean)
    sd=apply(new_X[node_ID1,node_ID2,],3,sd)
  }
  data=data.frame(effect=effect,prob=pb_sort,CI_lower=effect+sd,CI_upper=effect-sd,group=group)


  cut=unique(pb_sort); ave=se=NULL
  cut=c(cut,1)
  for(i in 1:(length(cut)-1)){
    vector=effect[(pb_sort>=cut[i])&(pb_sort<cut[i+1])]
    ave=c(ave,mean(vector))
    if(length(vector)==1) se=c(se,0)
    else se=c(se,sd(vector))
  }
  summary=data.frame(effect=ave,se=se,prob=cut[-length(cut)],group=group)


  return(list("data"=data,"summary"=summary))
}



#################### Alternative methods : Regualr Lasso, CNN ############################
#################### Regular Lasso #######################################################


#' Logistic probability model via penalized maximum likelihood
#'
#' Fit a logistic probability model based on Lasso penalty
#' @param xvec An input matrix. Each row is a vectorized predictor.
#' @param y Binary response variable.
#' @param xnew New predictors in the test data. Organized as a matrix with each row being a data point.
#' @param lambda The regularization penalty.
#' @return The returned object is a list of components.
#' @return \code{B_est} - The estimated coefficient vector of linear predictor.
#' @return \code{prob} - The predicted probabilities for the test data.
#' @usage Lasso(xvec,y,xnew,lambda)
#' @export
#' @importFrom glmnet "glmnet"
#' @importFrom stats "predict"

Lasso=function(xvec,y,xnew,lambda){
  fit <- glmnet(xvec,y,alpha = 1, family = "binomial",type.measure = "mse",lambda=lambda)
  coef = coef(fit)
  length(which(abs(coef)>0))
  lassoB = as.vector(coef)[-1]
  prob = predict(fit,newx = xnew,s = lambda,type = "response")
  res=list();
  res$B_est=lassoB; res$prob=prob
  return(res)
}


#################### CNN ##################################################################

#' Convolutional Neural Network (CNN) with two hidden layers
#'
#' Implement a CNN with two hidden layers and ReLU activation.
#' @param X A list of matrix-valued predictors.
#' @param y Binary response variable.
#' @param X_new A list of new matrices in the test data.
#' @param plot.figure Option for plotting trajectory of accuracy over epochs.
#' @return The returned object is a list of components.
#' @return \code{prob} - The predicted probabilities for the test data.
#' @return \code{class} - The estimated binary response for the test data.
#' @return \code{history} - The trajectory of classification accuracy over epochs.
#' @return \code{acc} - The classification accuracy on test data.
#' @usage CNN(X,y,X_new,plot.figure = FALSE)
#' @import keras
#' @export

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



################################# Some basic functions #####################################

link=function(type=c("gap","linear","logistic","sqrt","others"),beta=1){
  if(type=="gap") return(function(x) 0.6*(x>=0)-0.6*(x<0))
  if(type=="logistic") return(function(x) (1-exp(-x))/(1+exp(-x)))
  if(type=="linear") return(function(x) pmin(pmax(0,x/2+1/2),1))
  if(type=="sqrt") return(function(x)  pmin(pmax(0,sign(x)*x^2/2+1/2),1))
  else return(function(x) pmin(pmax(0,sign(x)*abs(x)^beta/2+1/2),1))
}


Makepositive = function(mat){
  h = eigen(mat,symmetric = T)
  nmat = (h$vectors)%*%diag(pmax(h$values,10^-4),nrow=nrow(mat))%*%t(h$vectors)
  return(nmat)
}

hinge = function(y) ifelse(1-y>0,1-y,0)

objective=function(b,yfit,ybar,Weight,type=c("logistic","hinge","binary")){
  if(type=="hinge")  return(sum(hinge(ybar*(yfit+b))*Weight))
  if(type=="logistic") return(sum(log2(1+exp(-ybar*(yfit+b)))*Weight))
  if(type=="binary") return(sum(abs(sign(ybar)-sign(yfit+b))*Weight))
}


########## Nonparametric Trace Regression via Sign Series Representation ################
################################### ASSIST  #############################################



#' Aggregation of structured sign series for trace regression (ASSIST)
#'
#' Main function for fitting the nonparametric trace regression. The algorithm uses a learning reduction approach to estimate the nonparametric trace regression via ASSIST.
#' @param X A list of matrix-valued predictors.
#' @param y A vector of response variables.
#' @param X_new A list of new matrices in the test data. \code{X_new = NULL} returns fitted values in the training data.
#' @param r The rank of sign representable function to be fitted.
#' @param sparse_r The number of zero rows in coefficient matrix.
#' @param sparse_c The number of zero columns in coefficient matrix.
#' @param H Resoution parameter that controls the number of classifiers to aggregate.
#' @param lambda Lagrangian multiplier.
#' @param rho.ini Initial step size.
#' @param min Minimum value of the response variables
#' @param max Maximum value of the response variables.
#' @return The returned object is a list of components.
#' @return \code{B_est} - An array that collects a series of coefficient matrices for the classifiers used in the algorithm.
#' @return \code{fitted} - The predicted responses in the test data.
#' @return \code{sign_fitted} - A matrix that collects a series of predicted signs for the classifiers used in the algorithm.
#' @usage TraceAssist(X,y,X_new=NULL,r,sparse_r,sparse_c,H=10,lambda=0,rho.ini=0.1,min,max)
#' @references C. Lee, L. Li, H. Zhang, and M. Wang (2021). Nonparametric Trace Regression via Sign Series Representation. \emph{arXiv preprint arXiv:2105.01783}.
#' @examples
#' ######### Generate matrices in the training data ################
#' X = list()
#' for(i in 1:10){
#'  X[[i]] = matrix(runif(4,-1,1),nrow = 2,ncol = 2)
#' }
#'
#' ######### Generate coefficient matrix ###########################
#' B = runif(2,-1,1)%*%t(runif(2,-1,1))
#'
#' ######### Generate response variables ###########################
#' y = NULL;signal = NULL
#' for(i in 1:10){
#'  signal = c(signal,sum(X[[i]]*B))
#'  y = c(y,sum(X[[i]]*B)+rnorm(1,sd = 0.1))
#' }
#'
#'
#' ######### Run ASSIST ############################################
#' res =TraceAssist(X,y,r = 1,sparse_r = 0,sparse_c = 0,min = min(y),max = max(y))
#' mean(abs(res$fitted-signal))
#'
#'
#' ######### Generate new matrices in the test data ################
#' X_new = list()
#' for(i in 1:10){
#'   X_new[[i]] = matrix(runif(4,-1,1),nrow = 2,ncol = 2)
#' }
#'
#' ######### Generate response variables from X_new ################
#' y_new = NULL
#' for(i in 1:10){
#'   y_new = c(y_new,sum(X_new[[i]]*B))
#' }
#'
#' ######### Run ASSIST #############################################
#' res =TraceAssist(X,y,X_new,r = 1,sparse_r = 0,sparse_c = 0,min = min(y),max = max(y))
#' mean(abs(res$fitted-y_new))
#'
#' @export
#' @importFrom quadprog "solve.QP"
#' @importFrom stats "optim"

TraceAssist=function(X,y,X_new=NULL,r,sparse_r,sparse_c,H=10,lambda=0,rho.ini=0.1,min,max){

  pred=NULL;d1=dim(X[[1]])[1];d2=dim(X[[1]])[2]

  B_est=array(0,dim=c(d1,d2,2*H+1));
  pi_seq=seq(from=min,to=max,length=(2*H+1))

  for(h in 2:(2*H)){
    pi=pi_seq[h]
    c = ADMM(X,sign(y-pi),abs(y-pi),r = r,srow = sparse_r,scol =sparse_c,lambda=lambda,rho.ini=rho.ini)
    if(is.null(X_new)) pred=rbind(pred,sign(c$fitted))
    else pred=rbind(pred,unlist(lapply(X_new,function(x)sign(c$predict(x)))))
    B_est[,,h]=c$B
  }
  res=list();
  res$B_est=B_est; res$fitted=min+prob_est(pred)*(max-min); res$sign_fitted=pred;
  return(res)
}


#' ADMM algorithm for weighted classification
#'
#' Implement an ADMM algorithm to optimize the weigthed classificiation loss.
#' @param X A list of matrix-valued predictors.
#' @param ybar A vector of  shifted response variables.
#' @param Weight Classification weight.
#' @param Covariate Additional covariates including intercept. \code{Covariate = NULL} indicates no covariates.
#' @param r The rank of coefficient matrix to be fitted.
#' @param srow The number of zero rows in coefficient matrix.
#' @param scol The number of zero columns in coefficient matrix.
#' @param lambda Lagrangian multiplier. Default is zero.
#' @param rho.ini Initial step size. Default is 1.
#' @return The returned object is a list of components.
#' @return \code{intercept} - The estimated intercept of the classifier.
#' @return \code{P_row} - The left-singular vectors of the coefficient matrix.
#' @return \code{P_col} - The right-singular vectors of the coefficient matrix.
#' @return \code{obj} - Trajectory of weighted classification loss values over iterations.
#' @return \code{iter} - The number of iterations.
#' @return \code{fitted} - A vector of fitted reponse from estimated classifier.
#' @return \code{B} - The estimated coefficient matrix of the classifier.
#' @usage ADMM(X,ybar,Weight,Covariate=NULL,r,srow,scol,lambda=0,rho.ini=1)
#' @references C. Lee, L. Li, H. Zhang, and M. Wang (2021). Nonparametric Trace Regression via Sign Series Representation. \emph{arXiv preprint arXiv:2105.01783}.
#' @examples
#' #### Generate matrix predictors  ##########
#' X = list()
#' for(i in 1:10){
#'  X[[i]] = matrix(runif(4,-1,1),nrow = 2,ncol = 2)
#' }
#'
#' #### Generate coefficient matrix #########
#' B = runif(2,-1,1)%*%t(runif(2,-1,1))
#'
#' #### Generate response variables #########
#' y = NULL
#' for(i in 1:10){
#'  y = c(y,sign(sum(X[[i]]*B)+rnorm(1,sd = 0.1)))
#' }
#'
#' #### classification with equal weights #########
#' res = ADMM(X,y,rep(1,10),r = 1,srow = 0,scol = 0)
#'
#' ### Misclassification rate on training data ######
#' mean(sign(res$fitted)-y)
#' @export
#' @importFrom quadprog "solve.QP"
#' @importFrom stats "optim"


ADMM=function(X,ybar,Weight,Covariate=NULL,r,srow,scol,lambda=0,rho.ini=1){
  ### normalization
  average_scale=mean(unlist(lapply(X,function(x)norm(x,"F"))))
  X=lapply(X,function(x)x/average_scale)
  ####

  result=list();
  n=length(X);
  rho=rho.ini/n;  ## decrease with sample size

  Lambda=matrix(0,nrow=nrow(X[[1]]),ncol=ncol(X[[1]]))

  PQ=0;
  #PQ=B.ini;

  iter=0; obj = 10^10;error=10;residual=10^10
  rho_list=NULL;
  while(((iter < 10)|(error > 10^-4))){

    PQ_prev=PQ

    ### update B
    res=SVM_offset(X,ybar,Weight,Covariate,OffsetC=(2*rho*PQ-Lambda)/(2*(lambda+rho)),cost=1/(2*(rho+lambda)))
    obj=c(obj,res$hinge) ## minimize objective
    B=res$coef

    ## Update PQ
    PQ=sparse_matrix(B+1/(2*rho)*Lambda,r,srow,scol)

    residual=c(residual,norm(B-PQ,"F")) ## primal residual
    #residual_dual=c(residual_dual,rho*norm(PQ-PQ_prev,"F")) ## dual residual

    ## update Lambda
    Lambda=Lambda+2*rho*(B-PQ)

    ## geometric step size
    if(iter>=10) {
      rho=rho*1.1
      lambda=lambda*1.1
    }

    rho_list=c(rho_list,rho)
    iter=iter+1;

    ##error=abs(-residual[iter+1]+residual[iter])
    error=residual[iter+1]/(r*(nrow(X[[1]])+ncol(X[[1]])-srow-scol));
    if(iter>=50) break
  }


  ### find intercept based on exact decomposition PQ (not B!)
  if(length(Covariate)>=1){
    W=cbind(1,matrix(unlist(Covariate),nrow=n,byrow=TRUE))
  } else {
    W=as.matrix(rep(1,n))}

  yfit = unlist(lapply(X,function(x) sum(PQ*x)))
  b=rep(0,dim(W)[2])
  intercept_opt=optim(b,function(x)objective(W%*%x,yfit,ybar,Weight=Weight,type="hinge"),method="BFGS")
  intercept=intercept_opt$par
  yfit=yfit+intercept


  PQ=PQ*average_scale
  ## take PQ (exact low-rank + sparse) as the output
  predictor = function(Xnew,Covariate=NULL){
    if(length(Covariate)>=1){
      W=cbind(1,matrix(unlist(Covariate),nrow=1,byrow=TRUE))
      sign(sum(Xnew*PQ)+intercept)
    }
    else sign(sum(Xnew*PQ)+intercept)
  }

  result$predict=predictor;## output function

  result$intercept=intercept;
  result$P_row=svd(PQ)$u[,1:r]; result$P_col=svd(PQ)$v[,1:r];
  result$obj=obj[-1];
  result$iter=iter;
  result$fitted=yfit;
  result$B=PQ;
  return(result)
}


##### Some inner functions for ASSIST algorithm ####################################################
SVM_offset=function(X,ybar,Weight,Covariate,OffsetC,cost=1){

  offset=unlist(lapply(X,function(x)sum(x*OffsetC)))


  n=length(X)
  Dmat=K=matrix(unlist(X),nrow=n,byrow=TRUE)%*%t(matrix(unlist(X),nrow=n,byrow=TRUE))
  dvec = 1-ybar*offset
  Dmat = Makepositive((ybar%*%t(ybar))*Dmat)
  Amat = cbind(ybar,diag(1,n),-diag(1,n))
  Weight[Weight==0]=min(Weight[Weight!=0]) ## avoid degenerate
  bvec = c(rep(0,1+n),-c(cost*Weight))


  res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)

  ## calculate coefficient
  coef=OffsetC
  coef=coef+matrix((ybar*res$solution)%*%matrix(unlist(X),nrow=n,byrow=TRUE),nrow=nrow(X[[1]]))

  ### calculate hinge loss

  yfit=K%*%(res$solution*ybar)+offset

  if(length(Covariate[[1]])>=1){
    W=cbind(1,matrix(unlist(Covariate),nrow=n,byrow=TRUE))
  } else {
    W=as.matrix(rep(1,n))}
  b=rep(0,dim(W)[2])

  ## calculate intercept
  intercept=optim(b,function(x)objective(W%*%x,yfit,ybar,Weight,type="hinge"),method="BFGS")$par

  yfit=yfit+W%*%intercept
  return(list("res"=res,"coef"=coef,"fitted"=yfit,"hinge"=objective(0,yfit,ybar,Weight,type="hinge"),"intercept"=intercept))
}

#### sparse matrix decomposition
sparse_matrix=function(M,r,srow,scol){
  srow_target=srow;scol_target=scol;

  ## scheme : first sparse, then low rank
  srow_index=scol_index=NULL
  M_low=M;srow=srow_target;scol=scol_target;step=0;
  while((srow>0)&(scol>0)){
    row_norm=diag(M_low%*%t(M_low))
    M_low[sort(row_norm,index=T)$ix[step+1],]=0
    srow_index=c(srow_index,sort(row_norm,index=T)$ix[step+1])
    srow=srow-1
    col_norm=diag(t(M_low)%*%M_low)
    M_low[,sort(col_norm,index=T)$ix[step+1]]=0
    scol_index=c(scol_index,sort(col_norm,index=T)$ix[step+1])
    scol=scol-1; step=step+1;
  }
  res=svd(M_low) ## some zero entries become nonzero due to numerical precision
  M_low=res$u[,1:r]%*%diag(res$d[1:r],nrow=r)%*%t(res$v[,1:r])
  M_low[srow_index,]=0 ## reset sparse entries to zero
  M_low[,scol_index]=0 #
  value2=sum((M-M_low)^2)
  # if(value2<=value1) M_low_output=M_low;
  M_low_output=M_low;
  return(M_low_output)
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


