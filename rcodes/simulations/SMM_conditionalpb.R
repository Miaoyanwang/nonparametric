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
colnames(dat) = c("z1","z2","class")
ggplot(data =dat, aes(x = z1, y = z2, colour = y))+geom_point()


##  true distribution of P(y = 1|(sum(x1),sum(x2)))
truecondp = function(x1,x2){
  pb = NULL
  for(i in 1:length(x1)){
    pb = c(pb,1/(1+exp(-sum(c(2,-2)*(c(x1[i],x2[i])-c(1,-1)))/4)))
  }
  return(pb)
}

x1_seq = seq(min(x[,1]),max(x[,1]),length.out = 20)
x2_seq = seq(min(x[,2]),max(x[,2]),length.out = 20)
x3_seq = seq(min(x[,3]),max(x[,3]),length.out = 20)
x4_seq = seq(min(x[,4]),max(x[,4]),length.out = 20)


z1range = range(x1_seq+x2_seq)
z1_seq = seq(z1range[1],z1range[2],length.out = 100)
z2range = range(x3_seq+x4_seq)
z2_seq = seq(z2range[1],z2range[2],length.out = 100)
zgrid = expand.grid(z1_seq,z2_seq)
ygrid_true = truecondp(zgrid[,1],zgrid[,2])
tgrid = as.data.frame(cbind(zgrid[,1],zgrid[,2],ygrid_true))

colnames(tgrid) = c("z1","z2","pb")
plt1 = ggplot(data = tgrid, aes(x = z1,y=z2,colour = pb))+geom_point()





#### Estimated distribution

fit = list()
precision = 0.1
for(i in 1:(1/precision-1)){
  fit[[i]] = smm(X,y,1,100,p = i*precision)$predict
}

estimate = function(x1,x2,x3,x4){
  pb = NULL
  for(i in 1:length(x1)){
    a = 1:(1/precision-1)
    for(j in 1:(1/precision-1)){
      a[j] = fit[[j]](matrix(c(x1[i],x2[i],x3[i],x4[i]),nrow = 2))
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

xgrid = expand.grid(x1_seq,x2_seq,x3_seq,x4_seq)
ygrid_est = estimate(xgrid[,1],xgrid[,2],xgrid[,3],xgrid[,4])

grid = as.data.frame(cbind(xgrid[,1]+xgrid[,2],xgrid[,3]+xgrid[,4],ygrid_est))
colnames(grid) = c("z1","z2","pb")
plt2 = ggplot(data = grid, aes(x = z1,y=z2,colour = pb))+geom_point()

ggarrange(plt1, plt2,labels = c("A", "B"),ncol = 2, nrow = 1)


