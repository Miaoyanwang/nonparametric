## nonparametric estimation of probability matrix, Aug 31, 2020
source("SMMKfunctions_con.R")
library("tensorordinal")
library(gridExtra)
set.seed(1)
m = 11; n = 11;b0 = 0;
## simulate X
X=list();index=0;
for(i in 1:m){
    for(j in 1:n){
        index=index+1;
        X[[index]]=matrix(0,nrow=m,ncol=n)
        X[[index]][i,j]=1
    }
}
## simulate probability matrix
#B1=matrix(runif(m*r,-1,1),nrow = m)
#B2=matrix(runif(m*r,0,1),nrow = m)
#signal=5*B1%*%t(B2)
#levelplot(signal)

signal=matrix(2,nrow=n,ncol=m)
signal[1:9,1:9]=-2
signal[1:5,1:5]=2
levelplot(signal)


###simulate probability based on logistic model
prob=exp(signal)/(1+exp(signal))
hist(prob)
## true step probability
H=10
prob_step=array(0,dim=c(n,m,H-1))
for(h in 1:(H-1)){
    prob_step[,,h]=1*(prob<h/H)
}


## simulate Y (1 replicate)
nrep=1
Y=array(0,dim=c(m,n,nrep))
for(rep in 1:nrep){
    Y[,,rep]=matrix(rbinom(m*n,1,prob),nrow=m,ncol=n,byrow=F)*2-1
}
levelplot(apply(Y,c(1,2),mean))
#levelplot(prob)

y_est=array(0,dim=c(n,m,nrep,H-1))
for(h in 1:(H-1)){
    for(rep in 1:nrep){
        con=SMMK_con(X,c(t(Y[,,rep])),r=1,kernel_row="linear",kernel_col="linear",cost = 1, rep = 2, p =1-h/H)
        y_est[,,rep,h]=matrix(con$fitted,nrow=m,ncol=n,byrow=T)
    }
}


est=prob_est(sign(y_est)[,,1,])

cum=array(0,dim=c(n,m,(H-1)))
for(h in 1:(H-1)){
    for(i in 1:n){
        for(j in 1:m){
            cum[i,j,h]=sum(sign(y_est)[i,j,1,1:h]==-1)/H
        }
    }
}


### approach 2: parametric model
data=array(0,dim=c(m,n,2))
data[,,1]=Y
data[,,2]=Y
res=fit_ordinal(data[,,1:2],c(4,4,2),omega=0,alpha=10)
est2=estimation(res$theta,res$omega,type="mode")
p = matrix(theta_to_p(res$theta,res$omega)[,1],nrow=n,ncol=m)
levelplot(1-p,main="parametric estimation", at=seq(0, 1, length.out=100),col.regions = gray(100:0/100))

## plot
p1=levelplot(prob,main="true probability",at=seq(0, 1, length.out=100),col.regions = gray(100:0/100),cex=2)
p2=levelplot((1+Y[,,1])*0.5,main="binary observation",at=seq(0, 1, length.out=100),col.regions = gray(100:0/100))
p3=levelplot(est,main="nonparametric estimation",at=seq(0, 1, length.out=100),col.regions = gray(100:0/100))
p4=levelplot(1-p,main="parametric (logistic-model) estimation", at=seq(0, 1, length.out=100),col.regions = gray(100:0/100))
grid.arrange(p1, p2, p3,p4,ncol=4)

levelplot(sign(prob_step),main="intermediate estimation: step function", at=seq(0, 1, length.out=100),col.regions = gray(100:0/100))
levelplot(cum,main="intermediate estimation: cum function",at=seq(0, 1, length.out=100),col.regions = gray(100:0/100))



#### help function
prob_est=function(est){
    n=dim(est)[1];m=dim(est)[2];H=dim(est)[3]+1
    prob_est=matrix(0,nrow=n,ncol=m)
    for(i in 1:n){
        for(j in 1:m){
            prob_est[i,j]=which.max(cumsum(est[i,j,]))/H
    }
    }
        return(prob_est)
}


