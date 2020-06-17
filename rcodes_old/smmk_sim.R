######################### conditional probability ####################
# P(x1|y=1)~ N((1,1),I), P(x2|y=1)~N((-1,-1),I)
# P(x1|y=-1)~ N((0,0),I), P(x2|y=1)~N((0,0),I)
N = 20
x = rbind(cbind(rnorm(N,1),rnorm(N,1),rnorm(N,-1),rnorm(N,-1)),
          cbind(rnorm(N),rnorm(N),rnorm(N),rnorm(N)))
X =  lapply(seq_len(nrow(x)),function(i) matrix(x[i,,drop = F],2,2))
y = c(rep(1,20),rep(-1,20))

### Realization of the data
dat = as.data.frame(cbind(x[,1]+x[,2],x[,3]+x[,4],y))
colnames(dat) = c("z1","z2","y")
ggplot(data =dat, aes(x = z1, y = z2, colour = y))+geom_point()

sym=function(x){
    d=dim(x)
    res=matrix(0,nrow=sum(d),ncol=sum(d))
    res[(1:d[1]),(d[1]+1):sum(d)]=x
    res[(d[1]+1):sum(d),(1:d[1])]=t(x)
    return(res)
}
sec=function(X,B){
    d=dim(X[[1]])
    B1=B[(1:d[1]),(d[1]+1):sum(d)]
    B2=B[(d[1]+1):sum(d),(1:d[1])]
    return(list(B1,t(B2)))
}


lresult = SMM(X,y,r)
presult = SMM(X,y,1,polykernel)

set.seed(1818)
m = 3; n = 4; r = 2; N = 50; b0 = 0.1
result = gendat(m,n,r,N,b0)
X = result$X; y = result$y; dat = result$dat
B = result$B
X_new=lapply(X,function(x)sym(x))
result=smm(X_new,y,r,drow=m)
sec(X,result$B)[[1]]/norm(sec(X,result$B)[[1]],"F")
sec(X,result$B)[[2]]/norm(sec(X,result$B)[[2]],"F")
B/norm(B,"F")


lresult = smm(X,y,r,rep=20)
lresult$B/norm(lresult$B,"F")
B/norm(B,"F")

lpredictor = lresult$predict
ppredictor = presult$predict

length(which(unlist(lapply(X,lpredictor))==y))/length(y)
length(which(unlist(lapply(X,ppredictor))==y))/length(y)


x1_seq = seq(min(x[,1]),max(x[,1]),length.out = 20)
x2_seq = seq(min(x[,2]),max(x[,2]),length.out = 20)
x3_seq = seq(min(x[,3]),max(x[,3]),length.out = 20)
x4_seq = seq(min(x[,4]),max(x[,4]),length.out = 20)
xgrid = expand.grid(x1_seq,x2_seq,x3_seq,x4_seq)

ylgrid = apply(xgrid,1,function(x) lpredictor(matrix(x,nrow = 2)))
ypgrid = apply(xgrid,1,function(x) ppredictor(matrix(x,nrow = 2)))

lgrid = as.data.frame(cbind(xgrid[,1]+xgrid[,2],xgrid[,3]+xgrid[,4],ylgrid))
colnames(lgrid) = c("z1","z2","y")
plt1 = ggplot(data = lgrid, aes(x = z1,y=z2,colour = as.factor(y)))+geom_point(size =0.01)+
  geom_point(data =dat,aes(x = z1,y = z2,colour = as.factor(y)))

pgrid = as.data.frame(cbind(xgrid[,1]+xgrid[,2],xgrid[,3]+xgrid[,4],ypgrid))
colnames(pgrid) = c("z1","z2","y")
plt2 = ggplot(data = pgrid, aes(x = z1,y=z2,colour = as.factor(y)))+geom_point(size =0.01)+
  geom_point(data =dat,aes(x = z1,y = z2,colour = as.factor(y)))

library(ggpubr)
ggarrange(plt1, plt2,labels = c("A", "B"),ncol = 2, nrow = 1)

