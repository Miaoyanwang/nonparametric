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


## questions to investigate:
## 1. extend simulation to low-rank matrix kernel.
## 2. Impact of regularization parameter, C (cost). Does C allow to vary for different pi? --> check Biometrika Vol 95 No 1 p149-167.
## 3. Is the linear kernel rich enough to approximate most functions if C is allowed to vary?
## 4. (Theory) what family of functions can be approximated by linear kernels? How about Gaussian kernels? Similar questions for matrix kernels. What function family can be well approxiamte by rank-r matrix kernels? How about nonlinear rank-r matrix kernel?


