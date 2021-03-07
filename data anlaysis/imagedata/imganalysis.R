#Data generation
library(imager)

im <- load.image("imgsample2.jpg")
im <-  grayscale(im)
plot(im,axes = F)
imarray = as.array(im)
dim(im)


imglist = list()
s = 0
for (i in 1:4) {
  for(j in 1:5){
    s = s+1
    set.seed(j)
    missing=array(rbinom(length(imarray),1,i/5),dim=dim(imarray))
    missing_img = imarray
    missing_img[missing==1]=NA
    immissing = as.cimg(missing_img)
    plot(immissing,axes = F)
    complete_img1 = complete_img2 = missing_img
    
    missing_v = array(missing_img[,,1,],dim = dim(im)[-3])
    imglist2[[s]] = missing_v
  }
}

timg = array(imarray[,,1,],dim = dim(im)[-3])
save(imglist,timg,file = "imglist.RData")

load(file ="imglist.RData")


## Output check#################################################################

load("imglist2.RData")

missing_img = complete_img1 = complete_img2 = imarray

missing_img[,,1,] = imglist2[[16]]
missing_img  = as.cimg(missing_img)
plot(missing_img,axes = F)


load(file = paste("imgoutput/img_matrix2_",16,"_",10,".RData",sep = ""))
complete_img1[,,1,] = A
complete_img1  = as.cimg(complete_img1)
plot(complete_img1,axes = F)

complete_img2[,,1,] = B
complete_img2  = as.cimg(complete_img2)
plot(complete_img2,axes = F)


## Output graph #######################################################################
index = matrix(nrow = 200,ncol = 2)
s = 0
for(i in 1:20){
  for(j in 1:10){
    s = s+1
    index[s,] = c(i,2*j)
  }
}

imgout = as.data.frame(matrix(nrow = 80,ncol = 5))
names(imgout) = c("Missing","Rank","MAE","sd","Method")
imgout[,1] = rep(rep(c("20%","40%","60%","80%"),each = 10),2)
imgout[,2] = rep(rep(2*(1:10),4),2)
imgout[,5] = rep(c("NonparaM","CPT"),each = 40)


m1 = NULL;m2 = NULL
for(s in 1:200){
  imgnum = index[s,1];r = index[s,2]
  if(imgnum==10&r ==14){
    load(file = paste("imgoutput/img_matrix2_",imgnum,"_",r-2,".RData",sep = ""))
    missing_v = imglist2[[imgnum]]
    naindex = which(is.na(missing_v))
    m1 = c(m1,mean(abs(timg[naindex]-A[naindex])));m2 = c(m2,mean(abs(timg[naindex]-B[naindex])))
  }else if(imgnum==19&r ==18){
    load(file = paste("imgoutput/img_matrix2_",imgnum,"_",r-2,".RData",sep = ""))
    missing_v = imglist2[[imgnum]]
    naindex = which(is.na(missing_v))
    m1 = c(m1,mean(abs(timg[naindex]-A[naindex])));m2 = c(m2,mean(abs(timg[naindex]-B[naindex])))
  }else{
  load(file = paste("imgoutput/img_matrix2_",imgnum,"_",r,".RData",sep = ""))
  missing_v = imglist2[[imgnum]]
  naindex = which(is.na(missing_v))
  m1 = c(m1,mean(abs(timg[naindex]-A[naindex])));m2 = c(m2,mean(abs(timg[naindex]-B[naindex])))
  }
}

M1 = matrix(m1,ncol = 10,byrow = T); M2 = matrix(m2,ncol = 10,byrow = T)
MAE1 = NULL; MAE2 = NULL
SD1 = NULL; SD2 = NULL
for(i in 1:4){
  MAE1 = c(MAE1,apply(M1[(5*i-4):(5*i),],2,mean)); MAE2 = c(MAE2,apply(M2[(5*i-4):(5*i),],2,mean))
  SD1  = c(SD1,apply(M1[(5*i-4):(5*i),],2,sd)); SD2  = c(SD2,apply(M2[(5*i-4):(5*i),],2,sd))
}
imgout[,3] = c(MAE1,MAE2)
imgout[,4] = c(SD1,SD2)

imgout$Missing = as.factor(imgout$Missing)
imgout$Method = as.factor(imgout$Method)

library(ggplot2)
ggplot(imgout[imgout$Missing=="20%",],aes(x = Rank,y = MAE,color = Method))+geom_point()+
  geom_line()+ ylim(range(imgout$MAE))+labs(title = "a. 20% Missing rate")+
  theme(text = element_text(size = 15))+geom_errorbar(aes(ymin=MAE-sd, ymax=MAE+sd), width=.2,position=position_dodge(0.05))   

ggplot(imgout[imgout$Missing=="40%",],aes(x = Rank,y = MAE,color = Method))+geom_point()+
  geom_line()+ ylim(range(imgout$MAE))+labs(title = "b. 40% Missing rate")+
  theme(text = element_text(size = 15))+geom_errorbar(aes(ymin=MAE-sd, ymax=MAE+sd), width=.2,position=position_dodge(0.05))  

ggplot(imgout[imgout$Missing=="60%",],aes(x = Rank,y = MAE,color = Method))+geom_point()+
  geom_line()+ ylim(range(imgout$MAE))+labs(title = "c. 60% Missing rate")+
  theme(text = element_text(size = 15))+geom_errorbar(aes(ymin=MAE-sd, ymax=MAE+sd), width=.2,position=position_dodge(0.05))  

ggplot(imgout[imgout$Missing=="80%",],aes(x = Rank,y = MAE,color = Method))+
  geom_point()+geom_line()+ ylim(range(imgout$MAE))+labs(title = "d. 80% Missing rate")+
  theme(text = element_text(size = 15))+geom_errorbar(aes(ymin=MAE-sd, ymax=MAE+sd), width=.2,position=position_dodge(0.05))    


mean(abs(A-timg))
mean(abs(B-timg))
A
B
MAE
