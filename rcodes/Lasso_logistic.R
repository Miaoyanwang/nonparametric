load(file = "bbnet68_spatial_orientation.RData")
bdat = read.table(file = "brain_vsplot.txt",head = T)
cost = 1


X = list()
for(n in 1:212){
  X[[n]] = A[,,n]
}
y = ifelse(VSPLOT=="low",-1,+1)

n = 212;d = 68

set.seed(1818)
# l1 = split(sample(1:57,57),as.factor(1:5))
# l2 = split(sample(58:114,57),as.factor(1:5))
l1 = split(sample(1:(n/2),n/2),as.factor(1:5))
l2 = split(sample((n/2+1):n,n/2),as.factor(1:5))
cindex = list()
for (k in 1:5) {
  cindex[[k]] = c(l1[[k]],l2[[k]])
}




x = matrix(nrow = n,ncol = d*(d-1)/2)
for(i in 1:n){
  x[i,] =  t(X[[i]])[lower.tri(t(X[[i]]))]
}


library(glmnet)

cv.out <- cv.glmnet(x,y,alpha = 1, family = "binomial",type.measure = "mse")
plot(cv.out)
lambda_lse <- cv.out$lambda.1se

coef = coef(cv.out,s = lambda_lse)
# the number of non zero coefficients
length(which(abs(coef)>0))

lassoB = matrix(0,nrow = d,ncol = d)
s = 0
for(i in 2:d){
  for(j in 1:(i-1)){
    s = s+1
    lassoB[j,i]= lassoB[i,j] = coef[s+1]
  }
}

nonzeroind = which(abs(lassoB)>0,arr.ind = T)
length(unique(nonzeroind[,1])) #non zero columns or rows => sparsity = 68-nonzero
ord = sort(unique(nonzeroind[,1])) 
length(which(svd(lassoB[ord,ord])$d>0.00001)) # rank
library('plot.matrix')
par(mar=c(6, 4, 4, 4) + 0.1)
plot(lassoB,col = hcl.colors(100, "YlOrRd", rev = T),main = "Submatrix of B")

# estimation result
prob = predict(cv.out,newx = x,s = lambda_lse,type = "response")
sum(log(prob)[which(y==1)])+sum(log(1-prob)[which(y== -1)])
library(ggplot2)
ggplot(bdat,aes(x = VSPLOT_TC,y = prob))+geom_point()+geom_hline(yintercept = 0.5,col = "red")+ylim(0,1)+
  labs(x= "VSPLOT",y = "Estimated probability")


# cross validation

cvresultlasso = matrix(nrow = 2,ncol = 5)
temp = NULL
for(k in 1:5){
  test_index = cindex[[k]]
  train_index = setdiff(1:n,test_index)
  x_train = x[train_index,]
  y_train = y[train_index]
  x_test = x[test_index,]
  y_test = y[test_index]
  cv.out <- cv.glmnet(x_train,y_train,alpha = 1, family = "binomial",type.measure = "mse")
  lambda_lse <- cv.out$lambda.1se
  prob = predict(cv.out,newx = x_test,s = lambda_lse,type = "response")
  #cvresultlasso[3,k] = sum(((y_test+1)/2-prob)^2)
  #mean(abs(y_test-ifelse(prob>0.5,1,-1))/2)
  #temp = c(temp,mean(abs(ifelse(prob>0.5,1,-1)-y_test))/2 )
  
  cvresultlasso[1,k] = (sum(log(prob)[which(y_test==1)])+sum(log(1-prob)[which(y_test== -1)]))
  
  prob = predict(cv.out,newx = x_train,s = lambda_lse,type = "response")
  #cvresultlasso[4,k] =  sum(((y_train+1)/2-prob)^2)
  cvresultlasso[2,k] = (sum(log(prob)[which(y_train==1)])+sum(log(1-prob)[which(y_train== -1)]))
}
apply(cvresultlasso,1,mean)



