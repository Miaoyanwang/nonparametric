source("SMMfunctions.R")

### data generation
set.seed(1818)
m = 4; n = 3; r = 2; N = 50; b0 = 0.1
result = gendat(m,n,r,N,b0)
X = result$X; y = result$y; dat = result$dat
B = result$B

ind = sample(100,100)

cvresult = matrix(nrow = 2, ncol = 5)
objresult = training=matrix(nrow = 2,ncol = 5)
for (i in 1:5) {
  testind = ind[((20*(i-1)+1):(20*i))]

  trX = X[-testind]
  result1 = svm(trX,y[-testind])
  cvresult[1,i] = length(which(unlist(lapply(X[testind],result1$predict))==y[testind]))/20
  objresult[1,i] = result1$obj
  #norm(result1$B/norm(result1$B,"F")-B/norm(B,"F"),"F")
  #training[1,i] = length(which(unlist(lapply(X[-testind],result1$predict))==y[-testind]))/80
  
  result = smm(trX,y[-testind],r=r,rep=10)
  cvresult[2,i] = length(which(unlist(lapply(X[testind],result$predict))==y[testind]))/20
  objresult[2,i] = objv(result$B,result$b0,trX,y[-testind])
  #norm(result$B/norm(result$B,"F")-B/norm(B,"F"),"F")
  #training[2,i] = length(which(unlist(lapply(X[-testind],result$predict))==y[-testind]))/80
  

}
cvresult;apply(cvresult,1,mean)
objresult


## data generation
u1 = matrix(runif(40,-1,1),nrow = 20,ncol = 2); v1 = matrix(runif(30,-1,1),nrow = 15,ncol = 2)
u2 = matrix(runif(40,-1,1),nrow = 20,ncol = 2); v2 = matrix(runif(30,-1,1),nrow = 15,ncol = 2)
X = list()
for (i in 1:50) {
  X[[i]] = u1%*%t(v1)+matrix(rnorm(300,0,4),nrow = 20,ncol = 15)
  X[[i+50]] = u2%*%t(v2)+matrix(rnorm(300,0,4),nrow = 20,ncol = 15)
}
y = c(rep(1,50),rep(-1,50))
x = matrix(nrow =100,ncol = 300)
for(i in 1:100){
  x[i,] = as.vector(X[[i]])
}
dat = data.frame(y = factor(y), x)



ind = sample(100,100)

cvresult = matrix(nrow = 2, ncol = 5)
objresult = matrix(nrow = 2,ncol = 5)
for (i in 1:5) {
  testind = ind[((20*(i-1)+1):(20*i))]

  trX = X[-testind]
  result = svm(trX,y[-testind])
  cvresult[1,i] = length(which(unlist(lapply(X[testind],result$predict))==y[testind]))/20
  objresult[1,i] = result$obj
  result = smm(trX,y[-testind],10)
  cvresult[2,i] = length(which(unlist(lapply(X[testind],result$predict))==y[testind]))/20
  objresult[2,i] = objv(result$B,result$b0,trX,y[-testind])
}
cvresult;apply(cvresult,1,mean)
objresult;apply(objresult,1,mean)




# SVM
fit = e1071::svm(factor(y) ~ ., data = dat, scale = FALSE, kernel = "linear", cost = 10)
fit$coefs
fit$rho
Bvec = t(fit$coefs)%*%fit$SV
hatB = matrix(Bvec,nrow = m,ncol = n)
#training result
length(which(y == predict(fit,dat[,-1])))/N


# SMM
result = smm(X,y,r,10);result
#traing result
length(which(unlist(lapply(X,result$predict))==y))/N


#obj value using svm
objv(hatB,fit$rho,X,y)
#obj value using smm
objv(result$B,result$b0,X,y)
#obj value at true B,b0
objv(B,b0,X,y)




###################040120##############################################
source("SMMfunctions.R")
set.seed(1818)
m = 10; n = 8; r = 1; N = 100; b0 = 0.1
result = gendat(m,n,r,N,b0)
X = result$X; y = result$y; dat = result$dat
B = result$B
k = 50
svmres = svm(X[1:k],y[1:k])
smmres = smm(X[1:k],y[1:k],5);smmres
subspace(smmres$B,B)
subspace(svmres$B,B)

lim = c(min(B/b0),max(B/b0))

par(mfrow = c(1,2))
plot(B/b0,svmres$B/abs(svmres$b0),xlab = "True B",ylab = "Estimated B",main = paste("SVM when N =", k),
     xlim = lim, ylim = lim)
abline(0,1,col = "red")
plot(B/b0,smmres$B/abs(smmres$b0),xlab = "True B",ylab = "Estimated B",main = paste("SMM when N =", k),
     xlim = lim, ylim = lim)
abline(0,1,col = "red")


#### consistency test
k = 50
par(mfrow = c(1,1))
svmres = svm(X[1:k],y[1:k])
fsmmres = smm(X[1:k],y[1:k],8)
plot(fsmmres$B/abs(fsmmres$b0),svmres$B/abs(svmres$b0),xlab = "SMM B",ylab = "SVM B",main = paste("When N =", k),
     xlim = c(min(fsmmres$B/abs(fsmmres$b0)),max(fsmmres$B/abs(fsmmres$b0))),
     ylim = c(min(svmres$B/abs(svmres$b0)),max(svmres$B/abs(svmres$b0))))
abline(0,1,col = "red")




#### Stability test
k = 100
smmres = smm(X[1:k],y[1:k],5,rep=1);smmres

lim = c(min(B/b0),max(B/b0))

par(mfrow = c(1,2))
plot(B/sqrt(sum(B^2)),smmres$B/sqrt(sum(smmres$B^2)),xlab = "True B",ylab = "Estimated B",main = paste("SMM when N =", k),
     xlim = lim, ylim = lim)
abline(0,1,col = "red")


smmres = smm(X[1:k],y[1:k],5,rep=1);smmres
plot(B/sqrt(sum(B^2)),smmres$B/sqrt(sum(smmres$B^2)),xlab = "True B",ylab = "Estimated B",main = paste("SMM when N =", k),
     xlim = lim, ylim = lim)
abline(0,1,col = "red")



