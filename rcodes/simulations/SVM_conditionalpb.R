############ simulation ###################
source("SMMfunctions.R")
x <- matrix(rnorm(30*2), ncol = 2)
y <- c(rep(-1,15), rep(1,15))
x[y==1,] <- x[y==1,] + 3/2
dat <- data.frame(x=x, y=ifelse(y==-1,0,y))
ggplot(data = dat,aes(x.1,x.2,colour = y))+geom_point()



X = lapply(seq_len(nrow(x)),function(i) x[i,,drop = F])
fit = svm(X,y,cost = 100)
bhat = fit$B; b0hat = fit$b0
fit1 = svm(X,y,cost = 100,p = 0.01)
bhat1 = fit1$B; b0hat1 = fit1$b0
fit2 = svm(X,y,cost = 100,p = 0.99)
bhat2 = fit2$B; b0hat2 = fit2$b0

ggplot(data = dat,aes(x.1,x.2,colour = as.factor(y)))+geom_point()+labs(colour = "classification")+
  geom_abline(intercept = -(b0hat)/bhat[2],slope = -bhat[1]/bhat[2],color = "blue")+
  geom_abline(intercept = -(b0hat1)/bhat1[2],slope = -bhat1[1]/bhat1[2],color = "red",lty = 2)+
  geom_abline(intercept = -(b0hat2)/bhat2[2],slope = -bhat2[1]/bhat2[2],color = "orange",lty = 2)


#####################################################################
ycondp = vector(length = 30)
for(i in 1:30){
  ycondp[i] = posterior(X,y,cost = 100,x[i,])
  print(paste(i,"th point is done lol"))
}

xgrid = expand.grid(X1 = x[,1], X2 = x[,2])
ygrid = vector(length = 900)
for(i in 1:900){
  ygrid[i] = posterior(X,y,cost = 100,xgrid[i,])
  print(paste(i,"th point is done lol"))
}


condpb = as.data.frame(cbind(xgrid,ygrid))

ggplot(data = condpb,aes(xgrid[,1],xgrid[,2],colour = ygrid))+geom_point(size = .5)+
  labs(colour = "probability",x ="x1",y = "x2")+
  scale_colour_gradientn(colours = c("darkblue","orange","darkgreen"),
                         breaks = breaks, labels = format(breaks))+
  geom_abline(intercept = -(b0hat)/bhat[2],slope = -bhat[1]/bhat[2],color = "orange")+
  geom_abline(intercept = -(b0hat1)/bhat1[2],slope = -bhat1[1]/bhat1[2],color = "darkblue")+
  geom_abline(intercept = -(b0hat2)/bhat2[2],slope = -bhat2[1]/bhat2[2],color = "darkgreen")+
  geom_point(data = dat,aes(x.1,x.2,colour = y),size = 3,shape = 5)


#########radial############# I used e1071 method
rm(svm)
library(e1071)
posterior2 = function(dat,test){
  a = 1:99
  for(i in 1:99){
    classcosts <- table(as.factor(dat$y))  # the weight vector must be named with the classes names
    classcosts[1] <- i# a class -1 mismatch has a terrible cost
    classcosts[2] <- 100 - i   # a class +1 mismatch not so much...

    fit = svm(factor(y) ~ ., data = dat, scale = FALSE, kernel = "radial",
              class.weights = classcosts)
    a[i] = ifelse(predict(fit, test)==1,1,-1)
  }
  if (all(a==1)) {
    return(1)
  }else if(all(a==-1)){
    return(0)
  }else{
    return((max(which(a==1))+min(which(a==-1)))/200)
  }
}



colnames(dat) = c("x.1","x.2","y")

xgrid = expand.grid(X1 = x[,1], X2 = x[,2])
colnames(xgrid) = c("x.1","x.2")
ygrid3 = vector(length = 900)
for(i in 1:900){
  ygrid3[i] = posterior2(dat,xgrid[i,])
  print(paste(i,"th point is done lol"))
}


condpb3 = as.data.frame(cbind(xgrid,ygrid3))

ggplot(data = condpb3,aes(xgrid[,1],xgrid[,2],colour = ygrid3))+geom_point(size = .5)+
  labs(colour = "probability",x ="x1",y = "x2")+
  scale_colour_gradientn(colours = c("darkblue","orange","darkgreen"),
                         breaks = breaks, labels = format(breaks))+
  geom_point(data = dat,aes(x.1,x.2,colour = y),size = 3)










