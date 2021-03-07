args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  BATCH <- as.numeric(args[1])
} else {
  stop()
}
source("signT.R")
load(file = "imglist.RData")




index = matrix(nrow = 200,ncol = 2)
s = 0
for(i in 1:20){
  for(j in 1:10){
    s = s+1
    index[s,] = c(i,2*j)
  }
}

imgnum = index[BATCH,1]
missing_v = imglist[[imgnum]]
r = index[BATCH,2]

Lmin = min(missing_v,na.rm = T)
Lmax = max(missing_v,na.rm = T)

missing_res = SignT(missing_v,r,Lmin=Lmin,Lmax=Lmax,H=10,option=1)
A = missing_res$est

est3 = fit_continuous(missing_v,r)
B = est3$est

rm(A)

naindex = which(is.na(missing_v))
MAE = c(mean(abs(timg[naindex]-A[naindex])),mean(abs(timg[naindex]-B[naindex])))
save(A,B,MAE,file = paste("img_matrix_",imgnum,"_",r,".RData",sep = ""))


