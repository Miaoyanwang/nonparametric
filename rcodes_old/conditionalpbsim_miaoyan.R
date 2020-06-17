############ simulation ###################
source("SMMfunctions.R")
library(ggplot2)
library(fields)
library(mvtnorm)

##################################################################
################# auxiliary functions  #########################
############################################################
### change the input format for posterior functions --> allow vector-type inputs
est=function(x.1,x.2,X,y){
return(apply(cbind(x.1,x.2),1,function(entry) posterior(X,y,cost=100,precision=0.1,test=entry)))
}
est_radial=function(x.1,x.2,x,y){
return(apply(cbind(x.1,x.2),1,function(entry) posterior2(data.frame(x=x, y=y),precision=0.1,test=data.frame(rbind(entry)))))
}

############## Simulation for probability estimation; added by Miaoyan on 04/10/2020 ######
### specify the true class probability over 2d space
prob=function(x1,x2){
pb=NULL
for(i in 1:length(x1)){
    ##pb=c(pb,0.2*(2*x1[i]+x2[i]>0)+0.8*(2*x1[i]+x2[i]<=0)) # linear prob.
pb=c(pb,dmvnorm(c(x1[i],x2[i]))/ dmvnorm(c(0,0))) ## Gaussian prob.
}
return(pb)
}


################# Classification ###########
N=70
x <- matrix(rnorm(N*2), ncol = 2)
y=rbinom(N,1,prob(x[,1],x[,2]))*2-1 ## generate class labels based on probability
dat <- data.frame(x=x, y=ifelse(y==-1,0,y))
ggplot(data = dat,aes(x.1,x.2,colour = y))+geom_point()

X = lapply(seq_len(nrow(x)),function(i) x[i,,drop = F])
fit = svm(X,y,cost = 100)
bhat = fit$B; b0hat = fit$b0
fit1 = svm(X,y,cost = 100,p = 0.01)
bhat1 = fit1$B; b0hat1 = fit1$b0
fit2 = svm(X,y,cost = 100,p = 0.99)
bhat2 = fit2$B; b0hat2 = fit2$b0

ggplot(data = dat,aes(x.1,x.2,colour = as.factor(y)))+geom_point()+labs(colour = "classification")+
geom_abline(intercept = -(b0hat)/bhat[2],slope = -bhat[1]/bhat[2],color = "blue")+
geom_abline(intercept = -(b0hat1)/bhat1[2],slope = -bhat1[1]/bhat1[2],color = "red",lty = 2)+
geom_abline(intercept = -(b0hat2)/bhat2[2],slope = -bhat2[1]/bhat2[2],color = "orange",lty = 2)

par(mfrow=c(3,1))
###################### Ground truth  ######################
x1_seq = x2_seq=seq(-2,2,length.out=20)
ygrid_true=outer(x1_seq,x2_seq,prob)
image.plot(x1_seq,x2_seq,ygrid_true,main="Ground truth probability",col  = gray((10:5)/10),xlab="x1",ylab="x2")

##################### Linear kernel #########################
ygrid=outer(x1_seq,x2_seq,function(x.1,x.2)est(x.1,x.2,X,y))
image.plot(x1_seq,x2_seq,ygrid,main="Estimation (Linear Kernel)",col= gray((10:5)/10),xlab="x1",ylab="x2")
points(dat$x.1,dat$x.2,col = as.factor(dat$y),pch=16)
abline(a= -(b0hat)/bhat[2],b= -bhat[1]/bhat[2],col = "orange",lwd=5)
abline(a= -(b0hat1)/bhat1[2],b= -bhat1[1]/bhat1[2],col = "darkblue",lwd=5)
abline(a= -(b0hat2)/bhat2[2],b= -bhat2[1]/bhat2[2],col = "darkgreen",lwd=5)

################### Gaussian kernel ################
ygrid_radial=outer(x1_seq,x2_seq,function(x.1,x.2)est_radial(x.1,x.2,x,y))
image.plot(x1_seq,x2_seq,ygrid_radial,main="Estimation (Gaussian Kernel)",col= gray((10:5)/10),xlab="x1",ylab="x2")
points(dat$x.1,dat$x.2,col = as.factor(dat$y),pch=16)
####################### end of simulation ####################


## What to do next:

## Classification:
## 1. (Method) Think about the implementation for non-linear low-rank matrix kernels. See my note.  
## 2. (Theory) Find relavent literature about convergence rate of misclassification error. Extend to linear low-rank matrix kernels.  

## Probability Estimation:
## 3. (Clarification) Are linear kernels (with varing C) rich enough to approximate most functions?

## 04/21: test kernel methods for matrix-valued predictors.

source("SMMfunctions.R")
set.seed(1818)
m = 5; n = 4; r = 2; N = 20; b0 = 0.1
result = gendat(m,n,r,N,b0)
X = result$X; y = result$y; dat = result$dat
B = result$B


result2 = smm(X,y,3)
library(psych)

kernel_tr=function(X,H,type="row"){
    m=nrow(X[[1]]);n=ncol(X[[1]]); N=length(X)
    output=matrix(nrow=N,ncol=N)
    
    H=H%*%solve(t(H)%*%H)%*%t(H) ## H is n-by-r or m-by-r
    
    if(type=="row"){
        for(i in 1:N){
            for(j in 1:N){
                #output[i,j]=tr(H%*%(X[[i]]%*%t(X[[j]])))
        output[i,j]=tr(H%*%expm(-t(X[[i]]-X[[j]])%*%(X[[i]]-X[[j]])))
            }
        }
    }
        if(type=="col"){
            for(i in 1:N){
                for(j in 1:N){
                    # output[i,j]=tr((t(X[[i]])%*%X[[j]])%*%H)
                     output[i,j]=tr(expm(-t(X[[i]]-X[[j]])%*%(X[[i]]-X[[j]]))%*%H)

                }
            }
        }
        h = eigen(output)
        output= (h$vectors)%*%diag(pmax(h$values,eps))%*%t(h$vectors)
        return(output)
}

svm_kernel=function(y,cost=10,kernel){
    N=length(y)
    dvec=rep(1,N)
    Amat = cbind(y,diag(1,N),-diag(1,N))
    bvec = c(rep(0,1+N),ifelse(y==1,-cost*.5,-cost*.5))
    
    val=NULL
    ### data generation
    set.seed(1818)
    V = randortho(n)[,1:r]
    
    for(nsim in 1:10){
    Dmat=y%*%t(y)*kernel_tr(X,V,type="col")
    alpha = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
    val=c(val,alpha$value)
    
    UTU=matrix(0,nrow=n,ncol=n)
    for(i in 1:N){
        for(j in 1:N){
            #UTU=UTU+alpha$solution[i]*alpha$solution[j]*y[i]*y[j]*t(X[[i]])%*%X[[j]]
            UTU=UTU+alpha$solution[i]*alpha$solution[j]*y[i]*y[j]*expm(-t(X[[i]]-X[[j]])%*%(X[[i]]-X[[j]]))
        }
    }
    UTU=solve(t(V)%*%V)%*%t(V)%*%UTU%*%V%*%solve(t(V)%*%V)
    
    UTH=array(0,dim=c(r,n,N))
    for(j in 1:N){
        for(i in 1:N){
            ##UTH[,,j]=UTH[,,j]+alpha$solution[i]*y[i]*solve(t(V)%*%V)%*%t(V)%*%t(X[[i]])%*%X[[j]]
            
            UTH[,,j]=UTH[,,j]+alpha$solution[i]*y[i]*solve(t(V)%*%V)%*%t(V)%*%expm(-t(X[[i]]-X[[j]])%*%(X[[i]]-X[[j]]))
        }
    }
    
    kernel_new=matrix(0,nrow=N,ncol=N)
    for(i in 1:N){
        for(j in 1:N){
            kernel_new[i,j]=tr(solve(UTU)%*%UTH[,,j]%*%t(UTH[,,i]))
        }
    }
    h = eigen(kernel_new,symmetric=T)
    kernel_new= (h$vectors)%*%diag(pmax(h$values,eps))%*%t(h$vectors)
    
    
    Dmat=y%*%t(y)*kernel_new
    beta = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
    val=c(val,beta$value)
    
    V_new=matrix(0,nrow=r,ncol=n)
    for(i in 1:N){
        V_new=V_new+beta$solution[i]*y[i]*solve(UTU)%*%UTH[,,i]
    }
    #V=t(V_new)
    
    V=svd(t(V_new))$u
    }
    
    
    
    
}
svm_kernel=function(y,cost=10,kernel){
    N=length(y)
    dvec=rep(1,N)
    Amat = cbind(y,diag(1,N),-diag(1,N))
    bvec = c(rep(0,1+N),ifelse(y==1,-cost*.5,-cost*.5))
    
    val=NULL
    U = randortho(m)[,1:r]

    for(nsim in 1:10){
    Dmat=y%*%t(y)*kernel_tr(X,U,type="row")
    alpha = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
    beta=y*alpha$solution
    val=c(val,alpha$value) ## maximize
    
    U_new=matrix(0,nrow=m,ncol=m)
    for(i in 1:N){
        for(j in 1:N){
            U_new=U_new+beta[i]*beta[j]*(X[[i]])%*%t(X[[j]])
        }
    }
    U_new=(1/C)*U_new
    
    
    for(i in 1:N){
        for(j in 1:N){
            U_new=U_new-(alpha$solution[i]!=0)*y[i]*(X[[i]])%*%t(X[[j]])*beta[j]
        }
    }
    
    Fun=function(U){tr(t(U)%*%(U_new%*%U))}
    df=function(U){U_new%*%U}
    U_est=optStiefel(Fun,df,Vinit=U,maxIters=10,maxLineSearchIters=10)
    U=U_est
    #step=svd(U_new)$d[1]
    #test=c(test,step)
    #U=svd(U-2/step*U_new%*%U)$u[,1:r]
    }
    
    

svm_kernel=function(y,cost=10,kernel){
    N=length(y)
    dvec=rep(1,N)
    Amat = cbind(y,diag(1,N),-diag(1,N))
    bvec = c(rep(0,1+N),ifelse(y==1,-cost*.5,-cost*.5))
    val=NULL
    set.seed(0)
    U = randortho(m)[,1:r] ##initilization
    V = randortho(n)[,1:r]
    
    for(nsim in 1:10){
    ##update U fixing V
    Dmat=y%*%t(y)*kernel_tr(X,V,type="col")
    alpha_new = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
    val=c(val,alpha_new$value)
    
    ##update V fixing U
    Dmat=y%*%t(y)*kernel_tr(X,U,type="row")
    alpha = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
    ##val=c(val,alpha$value)
    
    ## explicitly update U and V
    U_new=matrix(0,nrow=m,ncol=r)
    V_new=matrix(0,nrow=n,ncol=r)
    for(i in 1:N){
        for(j in 1:N){
    U_new=U_new+y[i]*y[j]*alpha$solution[i]*alpha_new$solution[j]*(X[[i]])%*%t(X[[j]])%*%U
    V_new=V_new+y[i]*y[j]*alpha$solution[i]*alpha_new$solution[j]*t(X[[i]])%*%X[[j]]%*%V
        }
    }
    U=svd(U_new)$u[,1:r]
    V=svd(V_new)$u[,1:r]
    
    ## get cost value
    
    }


