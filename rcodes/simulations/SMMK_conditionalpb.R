source("SMMKfunctions.R")


##### Simulation #########################################
library(mvtnorm)
S = diag(c(1,2))
x = rmvnorm(40, mean = rep(0, nrow(S)), sigma = S)
a = apply(x,1,function(x) dmvnorm(x, mean = rep(0, 2), sigma = S, log = FALSE))
y = sign(a-mean(a))
X = lapply(seq_len(nrow(x)),function(i) matrix(x[i,,drop = F],2,1))

dat = as.data.frame(cbind(x,y))
colnames(dat) = c("x1","x2","y")
ggplot(data =dat, aes(x = x1, y = x2, colour = as.factor(y)))+geom_point()

x1range = range(x[,1])
x1_seq = seq(x1range[1],x1range[2],length.out = 100)
x2range = range(x[,2])
x2_seq = seq(x2range[1],x2range[2],length.out = 100)
xgrid = expand.grid(x1_seq,x2_seq)

lprd = smmk(X,y,1,"linear")$predict
pprd = smmk(X,y,1,"poly")$predict
eprd = smmk(X,y,1,"exp")$predict

ylgrid = apply(xgrid,1,function(x) lprd(matrix(x,nrow = 2)))
ypgrid = apply(xgrid,1,function(x) pprd(makesym(matrix(x,nrow = 2))))
yegrid = apply(xgrid,1,function(x) eprd(makesym(matrix(x,nrow = 2))))
ytgrid = apply(xgrid,1,
               function(x) sign(dmvnorm(x,mean = rep(0,2),sigm = S,log = F)-mean(a)))

lgrid = as.data.frame(cbind(xgrid,ylgrid))
colnames(lgrid) = c("x1","x2","y")
plt1 = ggplot(data = lgrid, aes(x = x1,y=x2,colour = as.factor(y)))+geom_point(size =0.1)+
  geom_point(data =dat,aes(x = x1,y = x2,colour = as.factor(y)))

pgrid = as.data.frame(cbind(xgrid,ypgrid))
colnames(pgrid) = c("x1","x2","y")
plt2 = ggplot(data = pgrid, aes(x = x1,y=x2,colour = as.factor(y)))+geom_point(size =0.1)+
  geom_point(data =dat,aes(x = x1,y = x2,colour = as.factor(y)))

egrid = as.data.frame(cbind(xgrid,yegrid))
colnames(egrid) = c("x1","x2","y")
plt3 = ggplot(data = egrid, aes(x = x1,y=x2,colour = as.factor(y)))+geom_point(size =0.1)+
  geom_point(data =dat,aes(x = x1,y = x2,colour = as.factor(y)))

tgrid = as.data.frame(cbind(xgrid,ytgrid))
colnames(tgrid) = c("x1","x2","y")
plt4 = ggplot(data = tgrid, aes(x = x1,y=x2,colour = as.factor(y)))+geom_point(size =0.1)+
  geom_point(data =dat,aes(x = x1,y = x2,colour = as.factor(y)))


###################################################################################
## conditional pb on poly
fit = list()
precision = 0.1
for(i in 1:(1/precision-1)){
  fit[[i]] = smmk(X,y,1,"poly",p = i*precision)$predict
}

estimate = function(x1,x2){
  pb = NULL
  for(i in 1:length(x1)){
    a = 1:(1/precision-1)
    for(j in 1:(1/precision-1)){
      a[j] = fit[[j]](makesym(matrix(c(x1[i],x2[i]),nrow = 2)))
    }
    if (all(a==1)) {
      cpb = 1
    }else if(all(a==-1)){
      cpb = 0
    }else{
      cpb = (max(which(a==1))+min(which(a==-1)))/(2/precision)
    }
    pb = c(pb,cpb)
  }
  return(pb)
}

ypgrid_est = estimate(xgrid[,1],xgrid[,2])

pgrid_est = as.data.frame(cbind(xgrid[,1],xgrid[,2],ypgrid_est))
colnames(pgrid_est) = c("x1","x2","pb")
plt5 = ggplot(data = pgrid_est, aes(x = x1,y=x2,colour = pb))+geom_point()

###################################################################################
## conditional pb on exp
fit = list()
precision = 0.1
for(i in 1:(1/precision-1)){
  fit[[i]] = smmk(X,y,1,"exp",p = i*precision)$predict
}

estimate = function(x1,x2){
  pb = NULL
  for(i in 1:length(x1)){
    a = 1:(1/precision-1)
    for(j in 1:(1/precision-1)){
      a[j] = fit[[j]](makesym(matrix(c(x1[i],x2[i]),nrow = 2)))
    }
    if (all(a==1)) {
      cpb = 1
    }else if(all(a==-1)){
      cpb = 0
    }else{
      cpb = (max(which(a==1))+min(which(a==-1)))/(2/precision)
    }
    pb = c(pb,cpb)
  }
  return(pb)
}

yegrid_est = estimate(xgrid[,1],xgrid[,2])

egrid_est = as.data.frame(cbind(xgrid[,1],xgrid[,2],yegrid_est))
colnames(egrid_est) = c("x1","x2","pb")
plt6 = ggplot(data = egrid_est, aes(x = x1,y=x2,colour = pb))+geom_point()

ggarrange(plt4, plt1,plt2,plt3,labels = c("A", "B","C","D"),ncol = 2, nrow = 2)

ggarrange(plt5,plt6,labels = c("A","B"),ncol = 2)


