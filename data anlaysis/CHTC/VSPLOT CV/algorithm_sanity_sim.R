#############################algorithm sanity check ###################################################################
X = list()
n = 100; d = 10
for(i in 1:n){
  X[[i]] = matrix(nrow =d,ncol = d)
  for(j in 1:d){
    for(k in j:d){
      X[[i]][j,k] = rbinom(1,1,.5)
      X[[i]][k,j] =  X[[i]][j,k]
    }
  }
  diag(X[[i]]) = 0
}


truer = 3
rU = matrix(runif(d*truer,-1,1),nrow = d)
rV = matrix(runif(d*truer,-1,1),nrow = d)

B = rU%*%t(rV)
sparse = 3
ord = sort(sample(1:10,sparse))

B[ord,] = 0;  B[,ord] = 0
B = 2*(B+t(B))
ystar = unlist(lapply(X,function(x) sum(B*x)))
length(which(sign(ystar)==1))

## logistic model
pb = sigmoid(ystar)
y = ifelse(sapply(pb,function(x) rbinom(1,1,x))==1,1,-1)
length(which(y==1))

save(B,pb,y,X,file = "sim.RData")
load("sim.RData")
plot(pb,y)
source("SMMKfunctions_con.R")

svd(B)$d # real rank is actually = 5









########### classification

r = 5;sparse = 3
a = SMMK_sparse(X,y,5,sparse =3,option = "exact")
b = SMMK_sparse(X,y,5,sparse = 3)
c = ADMM(X,y,r = 5,srow = 3,scol = 3)



par(mfrow = c(1,4))
plot(c(B)/sqrt(sum(B^2)),c(a$B)/sqrt(sum(a$B^2)),xlab = "true B",ylab = "esitmated B (SMMK opt1)")+abline(0,1,col = "red")
plot(c(B)/sqrt(sum(B^2)),c(b$B)/sqrt(sum(b$B^2)),xlab = "true B",ylab = "esitmated B (SMMK opt2)")+abline(0,1,col = "red")
plot(c(B)/sqrt(sum(B^2)),c(c$B)/sqrt(sum(c$B^2)),xlab = "true B",ylab = "esitmated B (ADMM)")+abline(0,1,col = "red")
plot(c(B)/sqrt(sum(B^2)),c(lassoB)/sum(lassoB^2),xlab = "true B",ylab = "esitmated B (Lasso)")+abline(0,1,col = "red")

logistic = function(x) exp(x)/(1+exp(x))
loss = function(pb,y) sum(log(pb)[which(y==1)])+sum(log(1-pb)[which(y== -1)])
loss(logistic(a$fitted),y)
loss(logistic(b$fitted),y)
loss(logistic(c$fitted),y)

subspace(B,a$B)
subspace(B,b$B)
subspace(B,c$B)





####### probability estimation
#################### cross validation ###############################
# I made 3 different comparison, SMMK opt1, SMMK opt2, ADMM. We can use ADMM only 

length(which(y==1))
set.seed(1818)
l1 = split(sample(which(y==1)),as.factor(1:5))
l2 = split(sample(which(y==-1)),as.factor(1:5))
cindex = list()
for (k in 1:5) {
  cindex[[k]] = c(l1[[k]],l2[[k]])
}

predy = list()
predy[[1]] =   matrix(nrow = 49, ncol = 100) #SMMK opt1
predy[[2]] =   matrix(nrow = 49, ncol = 100) #SMMK opt2
predy[[3]] =   matrix(nrow = 49, ncol = 100) #ADMM

predy1 = list(); predy2 = list(); predy3 = list()
train_predy1 = list(); train_predy2 = list(); train_predy3 = list()
for(k in 1:5){
  predy1[[k]] = matrix(nrow = 49,ncol =length(cindex[[k]]))
  predy2[[k]] = matrix(nrow = 49,ncol =length(cindex[[k]]))
  predy3[[k]] = matrix(nrow = 49,ncol =length(cindex[[k]]))
  
}
for(k in 1:5){
  train_predy1[[k]] = matrix(nrow = 49,ncol =100-length(cindex[[k]]))
  train_predy2[[k]] = matrix(nrow = 49,ncol =100-length(cindex[[k]]))
  train_predy3[[k]] = matrix(nrow = 49,ncol =100-length(cindex[[k]]))
}

# rank =1,3,5,7
for(r in c(1,3,5,7)){
  for(j in 0:(10-r)){
    sparse = j
    for(i in 1:49){
      p = i*0.02
      
      a = SMMK_sparse(X,y,r = r,sparse =sparse,option = "exact",p = p)
      b = SMMK_sparse(X,y,r = r,sparse = sparse,p = p)
      c = ADMM(X,y,r = r,srow = sparse,scol = sparse, p = p,option ="symmetric")
      
      predy[[1]][i,] = sign(a$fitted)
      predy[[2]][i,] = sign(b$fitted)
      predy[[3]][i,] = sign(c$fitted)
      
      for(k in 1:5){
        
        
        test_index = cindex[[k]]
        train_index = setdiff(1:100,test_index)
        train_X = X[train_index]; train_y =  y[train_index]
        test_X = X[test_index]; test_y = y[test_index]
        
        con1=SMMK_sparse(train_X,train_y,r = r,kernel_row="linear",kernel_col="linear",option = "exact",cost=1, rep = 1, p = p,sparse=sparse)
        con2=SMMK_sparse(train_X,train_y,r = r,kernel_row="linear",kernel_col="linear",cost=1, rep = 1, p = p,sparse=sparse)
        con3=ADMM(train_X,train_y,r = r,srow = sparse,scol = sparse,p = p)    
        predy1[[k]][i,] = unlist(lapply(test_X,function(x) con1$predict(x)))
        train_predy1[[k]][i,] = sign(con1$fitted)
        
        predy2[[k]][i,] = unlist(lapply(test_X,function(x) con2$predict(x)))
        train_predy2[[k]][i,] = sign(con2$fitted)
        
        predy3[[k]][i,] = unlist(lapply(test_X,function(x) con3$predict(x)))
        train_predy3[[k]][i,] = sign(con3$fitted)
        
        print(paste("cv",p,"-",r,"-",sparse,"-",k," is done",sep = ""))
        # get error for each CV
      }
    }
    save(predy,predy1,predy2,predy3,train_predy1,train_predy2,train_predy3,
         file = paste("rk",r,"sparse",sparse,".RData",sep= "")) 
  }
  
}



##### The result  #####################################################
prob = list()
prob1 = list();prob2 = list(); prob3 = list()
train_prob1 = list(); train_prob2 = list(); train_prob3 = list()
cvresult1 = matrix(nrow = 2,ncol = 5);cvresult2 = matrix(nrow = 2,ncol = 5);cvresult3 = matrix(nrow = 2,ncol = 5)

cvmat = matrix(nrow = 28*3,ncol = 5)
s = 0
for(i in c(1,3,5,7)){
  for(j in 0:(10-i)){
    s = s+1
    rank = i; sparse = j
    cvmat[s,c(1,2,5)]= c(rank,sparse,1)
    cvmat[s+28,c(1,2,5)] = c(rank,sparse,2)
    cvmat[s+56,c(1,2,5)] = c(rank,sparse,3)
    load(paste("rk",rank,"sparse",sparse,".RData",sep= ""))
    for(k in 1:5){
      testy = y[cindex[[k]]]
      trainy = y[setdiff(1:100,cindex[[k]])]
      
      prob1[[k]] = pmax(apply(predy1[[k]],2,function(x) convpb(x,49)),0.005)
      prob2[[k]] = pmax(apply(predy2[[k]],2,function(x) convpb(x,49)),0.005)
      prob3[[k]] = pmax(apply(predy3[[k]],2,function(x) convpb(x,49)),0.005)
      
      train_prob1[[k]] = apply(train_predy1[[k]],2,function(x) convpb(x,49))
      train_prob2[[k]] = apply(train_predy2[[k]],2,function(x) convpb(x,49))
      train_prob3[[k]] = apply(train_predy3[[k]],2,function(x) convpb(x,49))
      
      cvresult1[1,k] = sum(log(prob1[[k]])[which(testy==1)])+sum(log(1-prob1[[k]])[which(testy== -1)])
      cvresult1[2,k] = sum(log(train_prob1[[k]])[which(trainy==1)])+sum(log(1-train_prob1[[k]])[which(trainy== -1)])
      cvresult2[1,k] = sum(log(prob2[[k]])[which(testy==1)])+sum(log(1-prob2[[k]])[which(testy== -1)])
      cvresult2[2,k] = sum(log(train_prob2[[k]])[which(trainy==1)])+sum(log(1-train_prob2[[k]])[which(trainy== -1)])
      cvresult3[1,k] = sum(log(prob3[[k]])[which(testy==1)])+sum(log(1-prob3[[k]])[which(testy== -1)])
      cvresult3[2,k] = sum(log(train_prob3[[k]])[which(trainy==1)])+sum(log(1-train_prob3[[k]])[which(trainy== -1)])
    }
    cvmat[s,3:4] = apply(cvresult1,1,mean)
    cvmat[s+28,3:4] = apply(cvresult2,1,mean)
    cvmat[s+56,3:4] = apply(cvresult3,1,mean)
    
  }
}

cvmat = as.data.frame(cvmat)
colnames(cvmat) = c("rank","sparse","test","train","type")
cvmat$rank = as.factor(cvmat$rank)
cvmat$type = as.factor(cvmat$type)

cvmat[which.max(cvmat$test),]
library(ggplot2)
g1 = ggplot(data = cvmat[cvmat$type !=1,], aes(x = sparse,y = test,col = type))+geom_point()+geom_line()+geom_hline(yintercept = 20*log(0.5),linetype = "dotted")+
  geom_hline(yintercept = -8.095777,linetype ="dotted",col = "red")+
  facet_wrap(~ rank,ncol =4)+scale_color_discrete(name = "Type" ,labels = c("SMMK","ADMM"))+labs(y = "log-likelihood on test")
g1
g2 = ggplot(data = cvmat[cvmat$type !=1,], aes(x = sparse,y = train,col = type))+geom_point()+geom_line()+geom_hline(yintercept = 80*log(0.5),linetype = "dotted")+
  geom_hline(yintercept = -24.814460,linetype = "dotted",col = "red")+
  facet_wrap(~ rank,ncol = 4)+scale_color_discrete(name = "Type" ,labels = c("SMMK", "SMMK","ADMM"))+labs(y = "log-likelihood on test")
g2



cvmat[cvmat$type==1,][which.max(cvmat[cvmat$type==1,]$test),] 
cvmat[cvmat$type==2,][which.max(cvmat[cvmat$type==2,]$test),] #(5,2,-6.95,-6.12)
cvmat[cvmat$type==3,][which.max(cvmat[cvmat$type==3,]$test),] #(5,4,-6.12,-8.04)
