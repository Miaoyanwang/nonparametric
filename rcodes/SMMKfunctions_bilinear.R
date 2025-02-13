#SMMK
library(pracma)
library(quadprog)
library(psych)
library(MASS)
library(tictoc)
library(rTensor)
# Make sure Q matrix is positive definite
Makepositive = function(mat){
  h = eigen(mat,symmetric = T)
  nmat = (h$vectors)%*%diag(pmax(h$values,10^-4),nrow=nrow(mat))%*%t(h$vectors)
  return(nmat)
}


# save kernel values
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


## Some kernels (Expkernel does not work)
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

# Main function (not available for weighed loss yet).
## input rank r should be a vector, r=(r_row,r_col).
smmk_new = function(X,y,r,kernel_row = c("linear","poly","exp","const"),kernel_col = c("linear","poly","exp","const"),cost = 10,rep = 1,p = 0.5){
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
  
  ################ Initilization and prestore all pairwise kernels ##############
  ##set.seed(1)
  d_row = nrow(X[[1]]); d_col = ncol(X[[1]]); N = length(X)
  K_row = Karray(X,kernel_row,type="row")
  K_col = Karray(X,kernel_col,type="col")
  r_row=r[1]; r_col=r[2];
  
  compareobj=10^10
  for(nsim in 1:rep){
  error = 10; iter = 0; obj=NULL
  
  P_row = randortho(d_row)[,1:r_row,drop = F]
  P_col = randortho(d_col)[,1:r_col,drop = F]
  
  while((iter < 20)&(error >10^-4)){
      
      ######################### Step 1. implicitly update core #########################
      res=update_core(X,y,P_row,P_col,K_row,K_col,cost=cost,p=p,intercept=F)
      obj = c(obj,-res$value)
      alpha = res$alpha
        
        
     ######################### Step 2. update row projection #########################
     
     ############### auxiliary quantities. Notation follows from 0821.pdf ###############
    
     W_col=P_col%*%t(P_col)
     W_row=P_row%*%t(P_row)
     
     
     ## speed up
     col_sum=matrix(unfold(K_col,c(1,2),c(3,4))@data%*%c(W_col),nrow=N,ncol=N)
     M=sum(W_col)*K_row+col_sum%o%(matrix(1,nrow=d_row,ncol=d_row))
     A=ttl(M,list(t(as.matrix(alpha*y)),t(P_row)),ms=c(2,4))
     B=ttl(A,list(t(as.matrix(alpha*y)),t(P_row)),ms=c(1,3))
     
     B=as.matrix(B@data[1,1,,]);
     
     ################ kernel used in dual problem  ###############
     temp=unfold(ttm(A,sqrtm(Makepositive(B))$Binv,4),1,c(2,3,4))@data
     K=temp%*%t(temp)
      
      ################ kernel used in dual problem  ###############
      
      dvec = rep(1,length(X))
      Dmat = Makepositive((y%*%t(y))*K)
      Amat = cbind(y,diag(1,N),-diag(1,N))
      bvec = c(rep(0,1+N),ifelse(y==1,-cost*(1-p),-cost*p))
      res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      beta = res$solution
      obj = c(obj,-res$value)
     
     ################ update row projection ################
 P_rownew=ttl(M,list(t(as.matrix(beta*y)),t(as.matrix(alpha*y)),t(P_row)),ms=c(1,2,4))[1,1,,]@data ## speed up
    
     P_row=svd(P_rownew)$u
     W_row = P_row%*%t(P_row)
    
     ########### Step 3. implicite update Core  #########################
     res=update_core(X,y,P_row,P_col,K_row,K_col,cost=cost,p=p,intercept=F)
     obj = c(obj,-res$value)
     alpha = res$alpha
     
     ########### Step 4. update column projection ######################
     
     ############### auxiliary quantities. Notation follows from 0821.pdf ##########
     row_sum=matrix(unfold(K_row,c(1,2),c(3,4))@data%*%c(W_row),nrow=N,ncol=N)
     M=sum(W_row)*K_col+row_sum%o%(matrix(1,nrow=d_col,ncol=d_col))
     A=ttl(M,list(t(as.matrix(alpha*y)),t(P_col)),ms=c(2,4))
     B=ttl(A,list(t(as.matrix(alpha*y)),t(P_col)),ms=c(1,3))
     
     
     B=as.matrix(B@data[1,1,,]);
     
     ################ kernel used in dual problem  ###############
     temp=unfold(ttm(A,sqrtm(Makepositive(B))$Binv,4),1,c(2,3,4))@data
     K=temp%*%t(temp)

     
     dvec = rep(1,length(X))
     Dmat = Makepositive((y%*%t(y))*K)
     Amat = cbind(y,diag(1,N),-diag(1,N))
     bvec = c(rep(0,1+N),ifelse(y==1,-cost*(1-p),-cost*p))
     res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
     gamma = res$solution
     obj = c(obj,-res$value)
    
     ############### update column projection  ###############
     P_colnew=ttl(M,list(t(as.matrix(gamma*y)),t(as.matrix(alpha*y)),t(P_col)),ms=c(1,2,4))[1,1,,]@data ## speed up
     
     
     P_col=svd(P_colnew)$u
     W_col = P_col%*%t(P_col)
     
     iter=iter+1;
     
     if(iter==1) error=10
     else error=abs(obj[iter*4]-obj[(iter-1)*4])/obj[iter*4] ## each loop involves 4 alternating updates
    } ## stop here for single initilization, when rep = 1
  
    ## take optimum over multiple initilizations
    if(compareobj>obj[iter*4]){
        P_row_optimum=P_row; P_col_optimum=P_col;
        W_row_optimum=W_row; W_col_optimum=W_col;
        obj_optimum=obj;
        compareobj=obj[iter*4]
     }
}
  
  ## Outside the loop. Estimate the intercept
    res=update_core(X,y,P_row_optimum,P_col_optimum,K_row,K_col,cost=cost,p=p,intercept=TRUE)
    alpha=res$alpha
    b0hat=res$intercept
    obj_optimum=c(obj_optimum,-res$value)
    
    
   
    slope = function(nX){
        ## pairwise kernel between new X and all training X
        newK=rep(0,N)
        for(i in 1:N){
                newK[i] = sum(W_col_optimum)*sum(W_row_optimum*kernel_row(t(nX),t(X[[i]])))+sum(W_row_optimum)*sum(W_col_optimum*kernel_col(nX,X[[i]]))
        }
        sp=sum(alpha*y*newK)
        return(sp)
    }
    predictor = function(nX) sign(slope(nX)+b0hat)
   
   
   result$slope = slope; result$predict = predictor; result$error = error; result$obj = obj_optimum; result$iter = iter; result$P_col=P_col_optimum;result$P_row = P_row_optimum;result$alpha=alpha; result$intercept = b0hat; result$fitted=res$fitted
   
   return(result)
}


update_core=function(X,y,P_row,P_col,K_row,K_col,cost=10,p=0.5,intercept=TRUE){
    
    N=length(X)
    W_row = P_row%*%t(P_row)
    W_col = P_col%*%t(P_col)
 K=matrix(sum(W_col)*unfold(K_row,c(1,2),c(3,4))@data%*%c(W_row)+sum(W_row)*unfold(K_col,c(1,2),c(3,4))@data%*%c(W_col),nrow=N,ncol=N)

 

    ## update core based on the dual problem
    dvec = rep(1,length(X))
    Dmat = Makepositive((y%*%t(y))*K)
    Amat = cbind(y,diag(1,N),-diag(1,N))
    bvec = c(rep(0,1+N),ifelse(y==1,-cost*(1-p),-cost*p))
    res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
    alpha = res$solution
    value = res$value
    
    #penalty = t(alpha*y)%*%K%*%(alpha*y) ## equal to squared norm of C
    
    b0=yfit=NULL ## estimate the intercept when needed
    if(intercept==TRUE){
    yfit=K%*%(alpha*y)
    negative=min(yfit)
    positive=max(yfit)
    gridb0 = seq(from = -1-negative,to = 1-positive,length = 100)
    b0=gridb0[which.min(sapply(gridb0,function(b) sum(pmax(1-yfit*y-b*y,0))))]
    }
    
    ## return core, objective in the dual (negative primal obj), intercept
    return(list("alpha"=alpha,"value"=value,"intercept"=b0,"fitted"=yfit+b0))
}
