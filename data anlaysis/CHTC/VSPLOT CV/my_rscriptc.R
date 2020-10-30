library(pracma)
library(rTensor)
library(quadprog)
args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  BATCH <- as.numeric(args[1])
} else {
  stop()
}
source("SMMKfunctions_con.R")
load(file = "bbnet68_spatial_orientation.RData")
cost = 1

#  save the combinations of rank and sparsity 
indexset = matrix(nrow = 2346,ncol = 2)
s = 0
for(i in 1:68){
  for(j in 1:(68-i+1)){
    s = s+1
    indexset[s,1] = i
    indexset[s,2] = j
  }
}


# Assign rank and sparsity according to BATCH (job number)
r = indexset[BATCH,1]
sparse = indexset[BATCH,2]


hinge = function(y) ifelse(1-y>0,1-y,0)


X = list()
for(n in 1:212){
  X[[n]] = A[,,n]
}
y = ifelse(VSPLOT=="low",-1,+1)





# Define test and training index set for cross validation 
set.seed(1818)
l1 = split(sample(1:106,106),as.factor(1:5))
l2 = split(sample(107:212,106),as.factor(1:5))
cindex = list()
for (k in 1:5) {
  cindex[[k]] = c(l1[[k]],l2[[k]])
}

# save prediced y labels from full dataset given pi = 0.02,...,0.98
predy =   matrix(nrow = 49, ncol = 212)
# save perdiced y_test labels from learning training dataset given pi = 0.02,....,0.98
predy1 = list();
# save perdiced y_training labels from learning training dataset given pi = 0.02,....,0.98
train_predy1 = list(); 

for(k in 1:5){
  predy1[[k]] = matrix(nrow = 49,ncol =length(cindex[[k]]))
}
for(k in 1:5){
  train_predy1[[k]] = matrix(nrow = 49,ncol =212-length(cindex[[k]]))
}






for(i in 1:49){
  p = i*0.02
  # get fitted value from full datset
  b = ADMM(X,y,r = r,srow = sparse-1,scol = sparse-1, p = p,option ="sym")
  predy[i,] = sign(b$fitted)
  
  for(k in 1:5){
    
    # define test and training datset
    test_index = cindex[[k]]
    train_index = setdiff(1:212,test_index)
    train_X = X[train_index]; train_y =  y[train_index]
    test_X = X[test_index]; test_y = y[test_index]
    
    # learn from traning dataset
    con2=ADMM(train_X,train_y,r=r,srow=sparse-1,scol=sparse-1,p=p,option="sym") 
    
    # save predicted y_test given p = i*0.02 => can change intercept as follows
    # yfit = unlist(lapply(test_X,function(x) con2$slope(x)))
    # positive=min(yfit[y==1])
    # negative=max(yfit[y==-1])
    # if ((1-positive)<(-1-negative)) {
    #   intercept = -(positive+negative)/2
    # }else{
    #   gridb0 = seq(from = -1-negative,to = 1-positive,length = 100)
    #   intercept = gridb0[which.min(sapply(gridb0,function(b) objective(b,yfit,y,p = p)))]
    # }
    # predy1[[k]][i,] = unlist(lapply(test_X,function(x) con2$slope(x)+intercept))
    
    predy1[[k]][i,] = unlist(lapply(test_X,function(x) con2$predict(x)))
    
    # save predicted y_training given p = i*0.02
    train_predy1[[k]][i,] = sign(con2$fitted)
    
    print(paste("cv",p,"-",r,"-",sparse,"-",k," is done",sep = ""))
  }
}



save(predy,predy1,train_predy1,file = paste("brainfitcv_admm_",r,"_",sparse,".RData",sep = ""))


