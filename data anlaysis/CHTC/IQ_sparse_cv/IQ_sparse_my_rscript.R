library(pracma)
library(rTensor)
library(quadprog)
source("SMMKfunctions_con.R")

args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
    BATCH <- as.numeric(args[1])
} else {
    stop()
}

load(file = "brain_binary_IQ.RData")
hinge = function(y) ifelse(1-y>0,1-y,0)
X = list()
for(n in 1:114){
    X[[n]] = A[,,n]
}
y = ifelse(FSIQ>120,1,-1)

set.seed(1818)
l1 = split(sample(1:57,57),as.factor(1:5))
l2 = split(sample(58:114,57),as.factor(1:5))
cindex = list()
for (k in 1:5) {
    cindex[[k]] = c(l1[[k]],l2[[k]])
}

cvresult = matrix(nrow = 0, ncol = 7)
colnames(cvresult) = c("trainzloss","trainhloss","zloss","hloss","rank","sparse","k")

d=68
r = (1:d)[BATCH]
for(sparse in 1:(d-r+1)){
for(k in 1:5){
  cv=c(rep(0,4),r,sparse,k)

  test_index = cindex[[k]]
  train_index = setdiff(1:114,test_index)
  train_X = X[train_index]; train_y =  y[train_index]
  #con=SMMK_con(train_X,train_y,r,kernel_row="linear",kernel_col="linear",cost=1, rep = 10, p = .5,sparse=sparse-1)

con=ADMM(train_X,train_y,r=r,srow=sparse-1,scol=sparse-1,rho.ini=1,p=0.5)
  
  #0-1 loss
  predict = unlist(lapply(X,con$predict))
  cv[1]= mean(predict[train_index]!=y[train_index])
  cv[3] = mean(predict[test_index]!=y[test_index])
  
  #hinge loss
  hloss = hinge(unlist(lapply(X,function(x) con$slope(x)+con$intercept))*y)
  cv[2] = sum(hloss[train_index]);
  cv[4]= sum(hloss[test_index]);
  
  cvresult=rbind(cvresult,cv)
  
  print(paste("cv",r,"-",sparse,"-",k," is done",sep = ""))
  # get error for each CV
}
}

save(cvresult,file = paste("IQ_CV_",r,".RData",sep=""))

########  after obtain the results from server #### ####
h=NULL
for(r in 1:68){
   load(paste("IQ_CV_",r,".RData",sep = ""))
   h=rbind(h,cvresult)
}

m=matrix(0,nrow=68,ncol=68)
for(r in 1:68){
    for(s in 1:68){
        m[r,s]=mean(h[(h[,5]==r)&(h[,6]==s),4])
        }
}

fields::image.plot(m, add = F, col = hcl.colors(12, "YlOrRd", rev = TRUE),xlab="rank/68",ylab="sparsity proportion",main="hinge loss in test set (IQ)")
#contour(m,add=TRUE,levels = pretty(c(0,3), 6), labcex = 1)
dev.copy(pdf,"hinge_test.pdf")

##### sequence of classifications. Try different combination of (r,srow,scol). How to choose the best one? should we use the one that gives best classification with p=0.5? Or something else, e.g. K-L divergence based on prob(y_grid) vs. y?
## sepcial case (r,srow,scol)=(1,0,0) reduces to earlier low-rank version
## try (r,srow,scol)=(2,60,60)?
y_grid=NULL
for(h in 1:20){
    con=ADMM(X,y,r=2,srow=60,scol=60,rho.ini=1,p=0.05*h)
    y_grid=rbind(y_grid,c(sign(con$fitted)))
}
image(y_grid)
plot(FSIQ,prob_est(y_grid,option=1))


