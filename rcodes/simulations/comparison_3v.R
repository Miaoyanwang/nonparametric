library(ggplot2)
N = 100
x = rbind(cbind(rnorm(N,1,.5),rnorm(N,1,.5),rnorm(N,-1,.5),rnorm(N,-1,.5)),
          cbind(rnorm(N,sd = .5),rnorm(N,sd = .5),rnorm(N,sd = .5),rnorm(N,sd = .5)))
X =  lapply(seq_len(nrow(x)),function(i) matrix(x[i,,drop = F],2,2))
y = c(rep(1,N),rep(-1,N))

### Realization of the data
dat = as.data.frame(cbind(x[,1]+x[,2],x[,3]+x[,4],y))
colnames(dat) = c("z1","z2","y")
ggplot(data =dat, aes(x = z1, y = z2, colour = y))+geom_point()

makesym = function(mat){
  m = nrow(mat); n = ncol(mat)
  nmat = rbind(cbind(matrix(0,n,n),t(mat)),cbind(mat,matrix(0,m,m)))
  return(nmat)
}

nX = lapply(X,makesym)

r1 = smm(X,y,1,rep = 1); r1$B

objv(r1$B,r1$b0,X,y,cost = 10,prob =.5)

r1$obj
r2 = smm(nX,y,2,rep = 1); r2$B

objv(r2$B,r2$b0,nX,y,prob =.5)
B = makesym(r1$B)
B = matrix(0,4,4); B[1:2,3:4] = r2$B[1:2,3:4]; B[3:4,1:2] = r2$B[3:4,1:2]


min(sapply(seq(-10,10,by = .01),function(x) objv(B,x,nX,y,prob=.5)))

Bsum = t(r2$B[1:2,3:4])+r2$B[3:4,1:2]
min(sapply(seq(-10,10,by = .01),function(x) objv(Bsum,x,X,y,prob=.5)))


B = makesym(Bsum/2)
min(sapply(seq(-10,10,by = .01),function(x) objv(B,x,nX,y,prob=.5)))

r2 = smm(nX,y,2,rep = 50); r2$B
objv(r2$B,r2$b0,nX,y,cost = 10,prob =.5)
B = matrix(0,4,4); B[1:2,3:4] = r2$B[1:2,3:4]; B[3:4,1:2] = r2$B[3:4,1:2]

min(sapply(seq(-10,10,by = .01),function(x) objv(B,x,nX,y,prob=.5)))

Bsum = t(r2$B[1:2,3:4])+r2$B[3:4,1:2]
min(sapply(seq(-10,10,by = .01),function(x) objv(Bsum,x,X,y,prob=.5)))
B = makesym(Bsum/2)
min(sapply(seq(-10,10,by = .01),function(x) objv(B,x,nX,y,prob=.5)))



###################################################################################
N  = 25
x = rbind(cbind(rnorm(N,1,.5),rnorm(N,1,.5),rnorm(N,-1,.5),rnorm(N,-1,.5)),
          cbind(rnorm(N,sd = .5),rnorm(N,sd = .5),rnorm(N,sd = .5),rnorm(N,sd = .5)))
X =  lapply(seq_len(nrow(x)),function(i) matrix(x[i,,drop = F],2,2))
y = c(rep(1,N),rep(-1,N))
nX = lapply(X,makesym)





result1 = smm(X,y,1)
result2 = smm(nX,y,2,rep = 10)
result3 = SMMK_con(X,y,1,rep =20)


B = result1$B
BtildeB = result2$B
Pr = result3$P_row; Pc = result3$P_col; alpha = result3$alpha
W_row = Pr%*%t(Pr);W_col= Pc%*%t(Pc)
B1 =0; B2 = 0;
for(i in 1:50){
  B1= B1+ alpha[i]*y[i]*W_row%*%X[[i]]
  B2 = B2+ t(alpha[i]*y[i]*W_col%*%t(X[[i]]))
}
min(sapply(seq(-10,10,by = .01),function(x) objv(B1+B2,x,X,y,prob=.5)))
B1
B2
B1+B2
objv(B,result1$b0,X,y,prob = .5)




predictor1 = result1$predict
predictor2 = result2$predict
predictor3 = result3$classifier

length(which(unlist(lapply(X,predictor1))==y))/length(y)
length(which(unlist(lapply(nX,predictor2))==y))/length(y)
length(which(unlist(lapply(X,predictor3))==y))/length(y)



x1_seq = seq(min(x[,1]),max(x[,1]),length.out = 20)
x2_seq = seq(min(x[,2]),max(x[,2]),length.out = 20)
x3_seq = seq(min(x[,3]),max(x[,3]),length.out = 20)
x4_seq = seq(min(x[,4]),max(x[,4]),length.out = 20)
xgrid = expand.grid(x1_seq,x2_seq,x3_seq,x4_seq)

ygrid1 = apply(xgrid,1,function(x) predictor1(matrix(x,nrow = 2)))
ygrid2 = apply(xgrid,1,function(x) predictor2(makesym(matrix(x,nrow = 2))))
ygrid3 = apply(xgrid,1,function(x) predictor3(matrix(x,nrow = 2)))

#length 160000
length(which(ygrid1-ygrid2!=0))
length(which(ygrid2-ygrid3!=0))
length(which(ygrid1-ygrid3!=0))





grid1 = as.data.frame(cbind(xgrid[,1]+xgrid[,2],xgrid[,3]+xgrid[,4],ygrid1))
grid2 = as.data.frame(cbind(xgrid[,1]+xgrid[,2],xgrid[,3]+xgrid[,4],ygrid2))
grid3 = as.data.frame(cbind(xgrid[,1]+xgrid[,2],xgrid[,3]+xgrid[,4],ygrid3))

colnames(grid1) = colnames(grid2) = colnames(grid3) = c("z1","z2","y")

plt1 = ggplot(data = grid1, aes(x = z1,y=z2,colour = as.factor(y)))+geom_point(size =0.01)+
  geom_point(data =dat,aes(x = z1,y = z2,colour = as.factor(y)))

plt2 = ggplot(data = grid2, aes(x = z1,y=z2,colour = as.factor(y)))+geom_point(size =0.01)+
  geom_point(data =dat,aes(x = z1,y = z2,colour = as.factor(y)))

plt3 = ggplot(data = grid3, aes(x = z1,y=z2,colour = as.factor(y)))+geom_point(size =0.01)+
  geom_point(data =dat,aes(x = z1,y = z2,colour = as.factor(y)))



library(ggpubr)
ggarrange(plt1, plt2,plt3,labels = c("v1", "v1/sym","v2"),ncol = 3, nrow = 1)




