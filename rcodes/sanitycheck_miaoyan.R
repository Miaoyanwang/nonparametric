source("SMMKfunctions_miaoyan.R")
set.seed(1)
N = 10
x = rbind(cbind(rnorm(N,1,.5),rnorm(N,1,.5),rnorm(N,-1,.5),rnorm(N,-1,.5)),
cbind(rnorm(N,sd = .5),rnorm(N,sd = .5),rnorm(N,sd = .5),rnorm(N,sd = .5)))
X =  lapply(seq_len(nrow(x)),function(i) matrix(x[i,,drop = F],2,2))
y = c(rep(1,N),rep(-1,N))

set.seed(1)
presult = smmk_new(X,y,c(1,1),kernel_row="linear",kernel_col="linear",rep=2)
plot(presult$obj)

source("SMMKfunctions.R")
sresult = smmk(X,y,kernel="linear",r=1)

sresult = smmk(X,y,2)


lresult = smmk(X,y,1)
presult = smmk(X,y,1,polykernel)
eresult = smmk(X,y,1,expkernel)




source("SMMfunctions.R")
set.seed(1818)
m = 3; n = 3; r = 1; N = 50; b0 = 0.1
result = gendat(m,n,r,N,b0)
X = result$X; y = result$y; dat = result$dat
B = result$B
svmres = svm(X,y)
smmres = smm(X,y,1);
smmres
resid=svd(smmres$B)$u[,1]%*%t(svd(smmres$B)$u[,1])-svd(B)$u[,1]%*%t(svd(B)$u[,1])
sum(resid^2)

subspace(smmres$B,B)
subspace(svmres$B,B)

source("SMMKfunctions.R")
smmkres= smmk(X,y,r=1,kernel = "linear",cost = 10,rep = 5,p = .5)
smmkres$V

smmres = smm(X,y,1);

source("SMMKfunctions_miaoyan.R")
set.seed(1)
res = smmk_new(X,y,kernel_row="linear",kernel_col="linear",r=c(1,1),cost=2);
plot(res$obj)

## works
set.seed(1)
res = smmk_new(X,y,kernel_row="linear",kernel_col="const",r=c(1,1),cost=1);
plot(res$obj)

#Best1=0
#for(i in 1:N){
#Best1=Best1+res$alpha[i]*(res$P_row)%*%t(res$P_row)%*%X[[i]]
#}

## residual
resid=res$P_row%*%t(res$P_row)-svd(B)$u[,1]%*%t(svd(B)$u[,1])
sum(resid^2)

