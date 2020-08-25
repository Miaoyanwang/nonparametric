## generate data with linear kernel, Aug 25, 2020
source("SMMfunctions.R")
set.seed(1818)
m = 20; n = 20; r = 1; N = 100; b0 = 0.1
result = gendat(m,n,r,N,b0) ## please try simulaiton with real X (from HCP data) as features, and simuluate binary y based on our model. Find suitable hyper-parameters.
X = X_t=result$X; y = result$y; dat = result$dat; B = result$B
y_true=rep(0,N)
for(i in 1:N){
    X[[i]]=X[[i]]/sqrt(sum(X[[i]]^2)) ## rescale features to have norm 1
    y_true[i]=sum(B*X[[i]])
    X_t[[i]]=t(X[[i]])
}

source("SMMKfunctions_con.R")
con=SMMK_con(X,y,r=1,kernel_row="linear",kernel_col="linear",cost = 1, rep = 1, p = .5) ## try different cost, m, n, N, r
sum(sign(con$fitted)!=y) ## diagnostic
plot(y_true,con$fitted) ## diagnostic
plot(svd(B)$u[,1], con$P_row) ## estimation error
plot(svd(B)$v[,1], con$P_col) ## estimation error

