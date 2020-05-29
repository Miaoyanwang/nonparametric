## May 29, 2020. by Miaoyan. SDR. How does the choice of rank in SMM affect SDR?
library(MixMatrix)
library(rTensor)
source("../SMMfunctions.R")

N=100; d1=4; d2=4; r1=2; r2=2; Nsim=10
error_L=error_R=array(0,dim=c(4,Nsim))
precision = 0.05
set.seed(1)
for(nsim in 1:Nsim){
for(r in 1:4){
set.seed(nsim)
X = rmatrixnorm(N,mean=matrix(0,nrow=d1,ncol=d2),U=diag(d1),V=diag(d2))
X = lapply(seq(N), function(i) rbind(X[,,i])) ## convert to a list
noise =rnorm(N,0,0)
y = sign(unlist(lapply(X,function(x)truefun(x)$val))+noise)

## weighted SMM
normalize=normalize_X(X)
X_norm=normalize$X
Bhat=array(0,dim=c(d1,d2,(1/precision-1)))

for(i in 1:(1/precision-1)){
   fit=smm(X_norm,y,r,cost=1,p = i*precision)
   Bhat[,,i]=normalize$L%*%fit$B%*%normalize$R
}

res=tucker(as.tensor(Bhat),ranks=c(r1,r2,r1*r2))
R_est=res$U[[1]][,1:r1]
L_est=res$U[[2]][,1:r2]

error_R[r,nsim]=norm(R_est%*%t(R_est)-diag(c(1,1,0,0)),"F")
error_L[r,nsim]=norm(L_est%*%t(L_est)-diag(c(1,1,0,0)),"F")
}
}


### plot
R_true=L_true=diag(c(1,1,0,0))
levelplot(R_true, col.regions = gray(100:0/100),main="True column SDR")
levelplot(R_est%*%t(R_est), col.regions = gray(100:0/100),main="Estimated column SDR (r=3)")
levelplot(L_est%*%t(L_est), col.regions = gray(100:0/100),main="Estimated row SDR (r=3)")

################ auxiliary functions ################
truefun=function(X){
    coef=matrix(0,nrow=nrow(X),ncol=ncol(X))
    coef[1,1]=2
    coef[1,2]=2
    coef[2,1]=2
    
    val=sum(X*coef)
    
    return(list("val"=val))
}

normalize_X=function(X){
    ## input X: a list of matrices
    N=length(X)
    d1=nrow(X[[1]])
    d2=ncol(X[[1]])
    meanX_est=Reduce('+',X)/N
    X=lapply(X,function(x)(x-meanX_est))
    
    X_trans=lapply(X,function(x)(t(x)))
    matrix_R=matrix(unlist(X_trans),ncol=d2,byrow=TRUE)
    R=sqrtm(cov(matrix_R))$Binv
    
    matrix_L=matrix(unlist(X),ncol=d1,byrow=TRUE)
    L=sqrtm(cov(matrix_L))$Binv
    norm_X=lapply(X, function(x) L%*%x%*%R)
    return(list("X"=norm_X,"R"=R,"L"=L))
}
