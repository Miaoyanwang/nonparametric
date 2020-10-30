## Getting cvresult
indexset = matrix(nrow = 2346,ncol = 2)
s = 0
for(i in 1:68){
  for(j in 1:(68-i+1)){
    s = s+1
    indexset[s,1] = i
    indexset[s,2] = j
  }
}

# getting probability from set of labels given pi = 0.02,...,0.98
convpb = function(temp,H){ 
  temp=c(1,temp,-1)
  index1=max(which(temp==1))-1
  index2=max(which(rev(temp)==-1))-1
  return((index1+H-index2)/(2*(H+1)))
}





# combination from rank 1 to 3
cvmat = matrix(nrow = 201,ncol = 4)
s = 0
for(BATCH in 1:201){
  s = s+1
  r = indexset[BATCH,1];sparse = indexset[BATCH,2]

  load(paste("brainfitcv_admm_",r,"_",sparse,".RData",sep = ""))

  cvmat[s,1] = r; cvmat[s,2] = sparse; 
  prob1 = list();
  train_prob1 = list(); 
  cvresult1 = matrix(nrow = 2,ncol = 5);cvresult2 = matrix(nrow = 2,ncol = 5)
  for(k in 1:5){
    testy = y[cindex[[k]]]
    trainy = y[setdiff(1:212,cindex[[k]])]
    # to avoid pb = 1 or 0 I use pmin, pmax with threshold 0.02 and 0.98
    prob1[[k]] = pmin(pmax(apply(predy1[[k]],2,function(x) convpb(x,49)),0.02),0.98)
    train_prob1[[k]] = pmin(pmax(apply(train_predy1[[k]],2,function(x) convpb(x,49)),0.02),0.98)
    cvresult1[1,k] = sum(log(prob1[[k]])[which(testy==1)])+sum(log(1-prob1[[k]])[which(testy== -1)])
    cvresult1[2,k] = sum(log(train_prob1[[k]])[which(trainy==1)])+sum(log(1-train_prob1[[k]])[which(trainy== -1)])
  }
  
  cvmat[s,3:4] = apply(cvresult1,1,mean)
}


cvmat = as.data.frame(cvmat)
colnames(cvmat) = c("rank","sparse","test","train")
cvmat$rank = as.factor(cvmat$rank)
cvmat[which.max(cvmat[,4]),]

library(ggplot2)
# -27.10 is lasso cv performance
g1 = ggplot(data = cvmat,aes(x = sparse,y = test))+geom_point()+geom_line()+geom_hline(yintercept = -27.10557,linetype = "dashed",col = "red")+
  geom_hline(yintercept= (212*log(0.5))/5,linetype = "dashed")+facet_wrap(~ rank)
g1

# -81.4257 is lasso cv performance
g2 = ggplot(data = cvmat,aes(x = sparse,y = train))+geom_point()+geom_line()+geom_hline(yintercept = -81.42570,linetype = "dashed",col = "red")+
  geom_hline(yintercept= (4*212*log(0.5))/5,linetype = "dashed")+facet_wrap(~ rank)
g2
