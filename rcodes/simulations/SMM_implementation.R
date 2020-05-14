source("SMMfunctions.R")

###############################################
## CV
###############################################
### simulation 1
### data generation
library(kableExtra)
set.seed(1818)
m = 5; n = 4; r = 3; N = 100; b0 = 0.1
result = gendat(m,n,r,N,b0)
X = result$X; y = result$y; dat = result$dat
B = result$B


### CV
ind = sample(100,100)

cvresult = matrix(nrow = 2, ncol = 5)
objresult = matrix(nrow = 2,ncol = 5)
for (i in 1:5) {
  testind = ind[((20*(i-1)+1):(20*i))]

  trX = X[-testind]
  result = svm(trX,y[-testind])
  cvresult[1,i] = length(which(unlist(lapply(X[testind],result$predict))==y[testind]))/20
  objresult[1,i] = result$obj
  result = smm(trX,y[-testind],3)
  cvresult[2,i] = length(which(unlist(lapply(X[testind],result$predict))==y[testind]))/20
  objresult[2,i] = objv(result$B,result$b0,trX,y[-testind])
  print(paste(i,"-th is done"))
}
cv = cbind(cvresult,apply(cvresult,1,mean))
colnames(cv) = c("1st","2nd","3rd","4th","5th","average")
rownames(cv) = c("SVM","SMM")
kable(cv,"latex")
cvobj = cbind(objresult,apply(objresult,1,mean))
colnames(cvobj) = c("1st","2nd","3rd","4th","5th","average")
rownames(cvobj) = c("SVM","SMM")
kable(cvobj,"latex")




###############################################
## pca plotting
###############################################

## data generation

result1 = svm(X,y)
result2 = smm(X,y,3)
x = dat[,-1]


library(ggplot2)
coord = eigen(cov(x))$vectors
fcoef = t(coord)%*%t(x)
rownames(fcoef) = paste("PC",1:(m*n),sep = "")
ndat = as.data.frame(t(rbind(y,fcoef)))

xgrid1 = seq(min(fcoef[1,]),max(fcoef[1,]),length = 100)
xgrid2 = seq(min(fcoef[2,]),max(fcoef[2,]),length = 100)
xgrid3 = seq(min(fcoef[3,]),max(fcoef[3,]),length = 100)




#### TRUE pc1 vs pc2 and pc1 vs pc3
xgrid12 = expand.grid(X1 = xgrid1,X2 = xgrid2)
xgrid13 = expand.grid(X1 = xgrid1,X2 = xgrid3)


gridX12 = lapply(1:nrow(xgrid12),
               function(i) matrix(coord[,c(1,2)]%*%t(xgrid12[i,]),nrow = m,ncol = n))
gridX13 = lapply(1:nrow(xgrid13),
                 function(i) matrix(coord[,c(1,3)]%*%t(xgrid13[i,]),nrow = m,ncol = n))


ygrid12 = unlist(lapply(gridX12,function(x) sign(sum(B*x)+b0)))
ygrid13 = unlist(lapply(gridX13,function(x) sign(sum(B*x)+b0)))

totgrid12 = as.data.frame(cbind(xgrid12,ygrid12))
totgrid13 = as.data.frame(cbind(xgrid13,ygrid13))

library(ggpubr)
gt12 = ggplot(data = totgrid12,aes(x = X1,y = X2,colour = as.factor(ygrid12)))+geom_point(size = .1)+
  labs(colour = "y",x ="PC1",y = "PC2")+
  geom_point(data =ndat,aes(x = PC1,y = PC2,colour = as.factor(y)))
gt12

gt13 = ggplot(data = totgrid13,aes(x = X1,y = X2,colour = as.factor(ygrid13)))+geom_point(size = .1)+
  labs(colour = "y",x ="PC1",y = "PC3")+
  geom_point(data =ndat,aes(x = PC1,y = PC3,colour = as.factor(y)))
gt13



## SVM method

Bhat1 = result1$B; b0hat1 = result1$b0
ygrid12 = unlist(lapply(gridX12,function(x) sign(sum(Bhat1*x)+b0hat1)))
ygrid13 = unlist(lapply(gridX13,function(x) sign(sum(Bhat1*x)+b0hat1)))

totgrid12 = as.data.frame(cbind(xgrid12,ygrid12))
totgrid13 = as.data.frame(cbind(xgrid13,ygrid13))

library(ggpubr)
g12 = ggplot(data = totgrid12,aes(x = X1,y = X2,colour = as.factor(ygrid12)))+geom_point(size = .1)+
  labs(colour = "y",x ="PC1",y = "PC2")+
  geom_point(data =ndat,aes(x = PC1,y = PC2,colour = as.factor(y)))
g12

g13 = ggplot(data = totgrid13,aes(x = X1,y = X2,colour = as.factor(ygrid13)))+geom_point(size = .1)+
  labs(colour = "y",x ="PC1",y = "PC3")+
  geom_point(data =ndat,aes(x = PC1,y = PC3,colour = as.factor(y)))
g13



#### Linear SMM estimation pc1 vs pc2 and pc1 vs pc3

Bhat2 = result2$B; b0hat2 = result2$b0
ygrid12 = unlist(lapply(gridX12,function(x) sign(sum(Bhat2*x)+b0hat2)))
ygrid13 = unlist(lapply(gridX13,function(x) sign(sum(Bhat2*x)+b0hat2)))

totgrid12 = as.data.frame(cbind(xgrid12,ygrid12))
totgrid13 = as.data.frame(cbind(xgrid13,ygrid13))


gm12 = ggplot(data = totgrid12,aes(x = X1,y = X2,colour = as.factor(ygrid12)))+geom_point(size = .1)+
  labs(colour = "y",x ="PC1",y = "PC2")+
  geom_point(data =ndat,aes(x = PC1,y = PC2,colour = as.factor(y)))
gm12

gm13 = ggplot(data = totgrid13,aes(x = X1,y = X2,colour = as.factor(ygrid13)))+geom_point(size = .1)+
  labs(colour = "y",x ="PC1",y = "PC3")+
  geom_point(data =ndat,aes(x = PC1,y = PC3,colour = as.factor(y)))
gm13

ggarrange(gt12,g12, gm12,gt13,g13, gm13,
          labels = c("A", "B","C","D","E","F"),
          ncol = 3, nrow = 2)



###############simulation 2
############## new data set data generation

m = 5; n = 4; r = 3
u1 = matrix(runif(m*r,-1,1),nrow = m,ncol = r); v1 = matrix(runif(n*r,-1,1),nrow = n,ncol = r)
u2 = matrix(runif(m*r,-1,1),nrow = m,ncol = r); v2 = matrix(runif(n*r,-1,1),nrow = n,ncol = r)
X = list()
for (i in 1:50) {
  X[[i]] = u1%*%t(v1)+matrix(rnorm(m*n,0,2),nrow = m,ncol = n)
  X[[i+50]] = u2%*%t(v2)+matrix(rnorm(m*n,0,2),nrow = m,ncol = n)
}
y = c(rep(1,50),rep(-1,50))
x = matrix(nrow =100,ncol = m*n)
for(i in 1:100){
  x[i,] = as.vector(X[[i]])
}
dat = data.frame(y = factor(y), x)

############cv##########################

ind = sample(100,100)

cvresult = matrix(nrow = 2, ncol = 5)
objresult = matrix(nrow = 2,ncol = 5)
for (i in 1:5) {
  testind = ind[((20*(i-1)+1):(20*i))]

  trX = X[-testind]
  result = svm(trX,y[-testind],minimization = TRUE)
  cvresult[1,i] = length(which(unlist(lapply(X[testind],result$predict))==y[testind]))/20
  objresult[1,i] = result$obj
  result = smm(trX,y[-testind],r)
  cvresult[2,i] = length(which(unlist(lapply(X[testind],result$predict))==y[testind]))/20
  objresult[2,i] = objv(result$B,result$b0,trX,y[-testind])
  print(paste(i,"th is done"))
}

cv = cbind(cvresult,apply(cvresult,1,mean))
colnames(cv) = c("1st","2nd","3rd","4th","5th","average")
rownames(cv) = c("SVM","SMM")
kable(cv,"latex")
cvobj = cbind(objresult,apply(objresult,1,mean))
colnames(cvobj) = c("1st","2nd","3rd","4th","5th","average")
rownames(cvobj) = c("SVM","SMM")
kable(cvobj,"latex")

############### plotting #####################3


result1 = svm(X,y)
result2 = smm(X,y,3)
x = dat[,-1]


library(ggplot2)
coord = eigen(cov(x))$vectors
fcoef = t(coord)%*%t(x)
rownames(fcoef) = paste("PC",1:(m*n),sep = "")
ndat = as.data.frame(t(rbind(y,fcoef)))

xgrid1 = seq(min(fcoef[1,]),max(fcoef[1,]),length = 100)
xgrid2 = seq(min(fcoef[2,]),max(fcoef[2,]),length = 100)
xgrid3 = seq(min(fcoef[3,]),max(fcoef[3,]),length = 100)




#### TRUE pc1 vs pc2 and pc1 vs pc3
xgrid12 = expand.grid(X1 = xgrid1,X2 = xgrid2)
xgrid13 = expand.grid(X1 = xgrid1,X2 = xgrid3)


gridX12 = lapply(1:nrow(xgrid12),
                 function(i) matrix(coord[,c(1,2)]%*%t(xgrid12[i,]),nrow = m,ncol = n))
gridX13 = lapply(1:nrow(xgrid13),
                 function(i) matrix(coord[,c(1,3)]%*%t(xgrid13[i,]),nrow = m,ncol = n))


## SVM method

Bhat1 = result1$B; b0hat1 = result1$b0
ygrid12 = unlist(lapply(gridX12,function(x) sign(sum(Bhat1*x)+b0hat1)))
ygrid13 = unlist(lapply(gridX13,function(x) sign(sum(Bhat1*x)+b0hat1)))

totgrid12 = as.data.frame(cbind(xgrid12,ygrid12))
totgrid13 = as.data.frame(cbind(xgrid13,ygrid13))

library(ggpubr)
g12 = ggplot(data = totgrid12,aes(x = X1,y = X2,colour = as.factor(ygrid12)))+geom_point(size = .1)+
  labs(colour = "y",x ="PC1",y = "PC2")+
  geom_point(data =ndat,aes(x = PC1,y = PC2,colour = as.factor(y)))
g12

g13 = ggplot(data = totgrid13,aes(x = X1,y = X2,colour = as.factor(ygrid13)))+geom_point(size = .1)+
  labs(colour = "y",x ="PC1",y = "PC3")+
  geom_point(data =ndat,aes(x = PC1,y = PC3,colour = as.factor(y)))
g13



#### Linear SMM estimation pc1 vs pc2 and pc1 vs pc3

Bhat2 = result2$B; b0hat2 = result2$b0
ygrid12 = unlist(lapply(gridX12,function(x) sign(sum(Bhat2*x)+b0hat2)))
ygrid13 = unlist(lapply(gridX13,function(x) sign(sum(Bhat2*x)+b0hat2)))

totgrid12 = as.data.frame(cbind(xgrid12,ygrid12))
totgrid13 = as.data.frame(cbind(xgrid13,ygrid13))


gm12 = ggplot(data = totgrid12,aes(x = X1,y = X2,colour = as.factor(ygrid12)))+geom_point(size = .1)+
  labs(colour = "y",x ="PC1",y = "PC2")+
  geom_point(data =ndat,aes(x = PC1,y = PC2,colour = as.factor(y)))
gm12

gm13 = ggplot(data = totgrid13,aes(x = X1,y = X2,colour = as.factor(ygrid13)))+geom_point(size = .1)+
  labs(colour = "y",x ="PC1",y = "PC3")+
  geom_point(data =ndat,aes(x = PC1,y = PC3,colour = as.factor(y)))
gm13

ggarrange(g12, gm12,g13, gm13,
          labels = c("A", "B","C","D"),
          ncol = 2, nrow = 2)






