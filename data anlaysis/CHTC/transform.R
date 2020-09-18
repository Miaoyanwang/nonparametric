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


Makepositive = function(mat){
  h = eigen(mat,symmetric = T)
  nmat = (h$vectors)%*%diag(pmax(h$values,10^-4),nrow=nrow(mat))%*%t(h$vectors)
  return(nmat)
}


Karray = function(X,kernel = function(X1,X2) t(X1)%*%X2, type="col"){
  N = length(X);
  if(type=="row"){
    Xt=list();
    for(i in 1:N){
      Xt[[i]]=t(X[[i]])
    }
    X=Xt
  }
  
  m= nrow(X[[1]]); n = ncol(X[[1]]);
  K = array(dim = c(N,N,n,n))
  for(i in 1:N){
    for(j in 1:N){
      K[i,j, , ] = kernel(X[[i]],X[[j]])
    }
  }
  return(as.tensor(K))
}

# K = function(A,B,kernel = function(a,b) sum(a*b),type = c("row","col")){
#   if(nrow(A) != nrow(B)) stop("two matrix sizes are different")
#   if(ncol(A) != ncol(B)) stop("two matrix sizes are different")
#   d1 = nrow(A); d2 = ncol(A)
#   if (type == "row") {
#     Kab = matrix(nrow= d1,ncol = d1)
#     for(i in 1:d1){
#       for(j in 1:d1){
#         Kab[i,j] = kernel(A[i,],B[j,])
#       }
#     }
#   }else{
#     Kab = matrix(nrow= d2,ncol = d2)
#     for(i in 1:d2){
#       for(j in 1:d2){
#         Kab[i,j] = kernel(A[,i],B[,j])
#       }
#     }
#   }
#   return(Kab)
# }

expkernel = function(Y,Z){
  n = ncol(Y)
  A = matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      A[i,j] = t(Y[,i]-Z[,j])%*%(Y[,i]-Z[,j])
    }
  }
  return(exp(-A))
}

polykernel = function(Y,Z,deg = 3){
  n = ncol(Y)
  return((t(Y)%*%Z+matrix(1,n,n))^deg)
}


linearkernel = function(X1,X2) t(X1)%*%X2

constkernel = function(X1,X2) return(matrix(0,nrow=ncol(X1),ncol=ncol(X1)))

SMMK_con = function(X,y,r,kernel_row = c("linear","poly","exp","const"),kernel_col = c("linear","poly","exp","const"), cost = 10, rep = 1, p = .5){
  result = list()
  
  # Default is linear kernel.
  kernel_row <- match.arg(kernel_row)
  if (kernel_row == "linear") {
    kernel_row = linearkernel
  }else if(kernel_row == "poly"){
    kernel_row = polykernel
  }else if(kernel_row =="exp"){
    kernel_row = expkernel
  }else if(kernel_row == "const"){
    kernel_row = constkernel
  }
  
  kernel_col <- match.arg(kernel_col)
  if (kernel_col == "linear") {
    kernel_col = linearkernel
  }else if(kernel_col == "poly"){
    kernel_col = polykernel
  }else if(kernel_col =="exp"){
    kernel_col = expkernel
  }else if(kernel_col == "const"){
    kernel_col = constkernel
  }
  
  d1 = nrow(X[[1]]); d2 = ncol(X[[1]]); n = length(X)
  K_row = Karray(X,kernel_row,type="row")
  K_col = Karray(X,kernel_col,type="col")
  compareobj = 10^10
  
  
  
  for(nsim in 1:rep){
    error = 10; iter = 0; obj = 10^10
    # initialize P_row,P_col
    P_row = randortho(d1)[,1:r,drop = F]; P_col = randortho(d2)[,1:r,drop = F]  
    
    
    
    while((iter < 20)&(error >10^-3)){
      # update C
      W_row = P_row%*%t(P_row); W_col = P_col%*%t(P_col)
      #Dmat = matrix(nrow = n,ncol = n)
      #for(i in 1:n){
      #for(j in 1:n){
      #     Dmat[i,j] = sum(W_row*K_row[i,j,,])+sum(W_col*K_col[i,j,,])
      #  }
      #}
      ## loop is slow. Replace loop by tensor-matrix multiplication
      
      Dmat=matrix(unfold(K_row,c(1,2),c(3,4))@data%*%c(W_row)+unfold(K_col,c(1,2),c(3,4))@data%*%c(W_col),nrow=n,ncol=n)
      
      dvec = rep(1,n)
      Dmat = Makepositive((y%*%t(y))*Dmat)
      Amat = cbind(y,diag(1,n),-diag(1,n))
      bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
      res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      alpha=res$solution
      
      # update P
      #CC_row = 0; CPh_row = array(0,dim = c(n,r,d1))
      #CC_col = 0; CPh_col = array(0,dim = c(n,r,d2))
      #for(i in 1:n){
      # for(j in 1:n){
      #   CC_row = CC_row + y[i]*alpha[i]*y[j]*alpha[j]*K_row[i,j,,]
      #   CC_col = CC_col + y[i]*alpha[i]*y[j]*alpha[j]*K_col[i,j,,]
      # }
      #}
      
      ## no need to find explicite inverse
      #CCi_row = solve(t(P_row)%*%CC_row%*%P_row)
      #CCi_col = solve(t(P_col)%*%CC_col%*%P_col)
      
      ## CPh_row is an intermediate step for CC_row. No need to compute twice.
      #for(i in 1:n){
      # cph_row = 0; cph_col = 0;
      # for(j in 1:n){
      #   cph_row = cph_row + y[j]*alpha[j]*K_row[j,i,,]
      #   cph_col = cph_col + y[j]*alpha[j]*K_col[j,i,,]
      # }
      # CPh_row[i,,] = t(P_row)%*%cph_row
      # CPh_col[i,,] = t(P_col)%*%cph_col
      #}
      
      
      
      CPh_row=ttl(K_row,list(t(as.matrix(y*alpha)),t(P_row)),ms=c(1,3))
      CPh_col=ttl(K_col,list(t(as.matrix(y*alpha)),t(P_col)),ms=c(1,3))
      
      CC_row=ttl(CPh_row,list(t(as.matrix(y*alpha)),t(P_row)),ms=c(2,4))
      CC_col=ttl(CPh_col,list(t(as.matrix(y*alpha)),t(P_col)),ms=c(2,4))
      
      CC_row=as.matrix(CC_row@data[1,1,,])
      CC_col=as.matrix(CC_col@data[1,1,,])
      
      factor_row=unfold(ttm(CPh_row,sqrtm(Makepositive(CC_row))$Binv,3),2,c(1,3,4))@data
      factor_col=unfold(ttm(CPh_col,sqrtm(Makepositive(CC_col))$Binv,3),2,c(1,3,4))@data
      Dmat=factor_row%*%t(factor_row)+factor_col%*%t(factor_col)
      
      #Dmat = matrix(nrow = n,ncol = n)
      #for(i in 1:n){
      # for(j in 1:n){
      #   Dmat[i,j] = sum(CCi_row%*%CPh_row[i,,]*CPh_row[j,,])+
      #                            sum(CCi_col%*%CPh_col[i,,]*CPh_col[j,,])
      # }
      #}
      ### Use Gram matrix to efficiently compute Dmat
      
      
      
      dvec = rep(1,n)
      Dmat = Makepositive((y%*%t(y))*Dmat)
      Amat = cbind(y,diag(1,n),-diag(1,n))
      bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
      res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      alpha=res$solution
      obj=c(obj,-res$value)
      iter = iter+1
      error = abs(-obj[iter+1]+obj[iter])/obj[iter]
      
      # P formula
      #P_row = 0; P_col = 0
      #for(i in 1:n){
      # P_row = P_row + alpha[i]*y[i]*rbind(CPh_row[i,,])
      # P_col = P_col + alpha[i]*y[i]*rbind(CPh_col[i,,])
      #}
      
      #P_row = svd(t(P_row)%*%CCi_row)$u; P_col = svd(t(P_col)%*%CCi_col)$u
      ## No need to multiply CCi_row. Right multiplying a full-rank matrix will not change the left singular vector
      
      P_row=ttm(CPh_row,t(as.matrix(alpha*y)),2)[1,1,,]@data
      P_col=ttm(CPh_col,t(as.matrix(alpha*y)),2)[1,1,,]@data
      
      P_row=matrix(P_row,nrow=r)
      P_col=matrix(P_col,nrow=r)
      
      P_row = svd(P_row)$v
      P_col = svd(P_col)$v
      
    }
    if(compareobj>obj[iter+1]){
      P_row_optimum=P_row; P_col_optimum=P_col;
      obj_optimum=obj;
      compareobj=obj[iter+1]
    }
  }
  
  
  ## move outside the loop. Only need to compute once for the optimal replicate.
  P_row= P_row_optimum; P_col= P_col_optimum;
  W_row = P_row%*%t(P_row); W_col = P_col%*%t(P_col)
  #Dmat = matrix(nrow = n,ncol = n)
  #for(i in 1:n){
  # for(j in 1:n){
  #   Dmat[i,j] = sum(W_row*K_row[i,j,,])+sum(W_col*K_col[i,j,,])
  # }
  #}
  Dmat=K=matrix(unfold(K_row,c(1,2),c(3,4))@data%*%c(W_row)+unfold(K_col,c(1,2),c(3,4))@data%*%c(W_col),nrow=n,ncol=n)
  
  dvec = rep(1,n)
  Dmat = Makepositive((y%*%t(y))*Dmat)
  Amat = cbind(y,diag(1,n),-diag(1,n))
  bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
  res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
  alpha=res$solution
  
  
  
  slope = function(Xnew){
    
    newK = rep(0,n)
    for( i in 1:n){
      newK[i] = sum(W_row*kernel_row(t(Xnew),t(X[[i]])))+
        sum(W_col*kernel_col(Xnew,X[[i]]))
    }
    
    return(sum(alpha*y*newK))
  }
  
  # intercept part estimation (update b)
  #positive = min(unlist(lapply(X,slope))[which(y==1)])
  #negative = max(unlist(lapply(X,slope))[which(y==-1)])
  
  yfit=K%*%(alpha*y) ## faster than lapply
  positive=min(yfit[y==1])
  negative=max(yfit[y==-1])
  intercept = -(positive+negative)/2
  
  compareobj = obj[iter+1]
  predictor = function(Xnew) sign(slope(Xnew)+intercept)
  
  result$alpha = alpha
  result$slope = slope; result$predict = predictor
  result$intercept = intercept
  result$P_row = P_row; result$P_col = P_col
  result$obj = obj[-1]; result$iter = iter; result$error = error; result$fitted=yfit+intercept ## add fitted value as a criterium to select cost
  return(result)
}








load(file = "bbnet68_spatial_orientation.RData")
hinge = function(y) ifelse(1-y>0,1-y,0)
X = list()
for(n in 1:212){
  X[[n]] = A[,,n]
}

mX = matrix(0,nrow(X[[1]]),ncol(X[[1]]))
for(n in 1:212){
  mX = mX+X[[n]]
}
mX = mX/212

# Centerize
Xc = list()
for(n in 1:212){
  Xc[[n]]= X[[n]]-mX
}


# Normalize
Xn = list()
for(n in 1:212){
  Xn[[n]] = X[[n]]/sqrt(sum((X[[n]])^2))
}



# Normalize + centerize
Xnc = list()
for(n in 1:212){
  Xnc[[n]] = Xc[[n]]/sqrt(sum((Xc[[n]])^2))
}

y = ifelse(VSPLOT=="low",-1,+1)


set.seed(1)
l1 = split(sample(1:106,106),as.factor(1:5))
l2 = split(sample(107:212,106),as.factor(1:5))
cindex = list()
for (k in 1:5) {
  cindex[[k]] = c(l1[[k]],l2[[k]])
}

indset = matrix(0,nrow =320 ,ncol = 3)
s = 0
for(k in 1:4){
  for(i in 1:20){
    for(j in 1:4){
      s = s+1
      indset[s,]= c(i,0.5*j,k)
    }
  }
}




cvresult = matrix(nrow = 1, ncol = 5)
colnames(cvresult) = c("zloss","hloss","rank","cost","method")
fitresult = matrix(nrow = 5,ncol = 212)


r = indset[BATCH,][1]
cost = indset[BATCH,][2]
method = indset[BATCH,][3]

if (method ==1) {
  X = X
}else if (method == 2) {
  X = Xc
}else if (method ==3) {
  X = Xn
}else{
  X = Xnc
}




cvresult[1,3:5] = c(r,cost,method)
hloss = 0;  zloss = 0
for(k in 1:5){
  test_index = cindex[[k]]
  train_index = setdiff(1:212,test_index)
  train_X = X[train_index]; train_y =  y[train_index]
  con=SMMK_con(train_X,train_y,r,kernel_row="linear",kernel_col="linear",cost, rep = 1, p = .5)
  
  #0-1 loss
  y_predict = unlist(lapply(X[test_index],con$predict))
  zloss = zloss + length(which(y_predict-y[test_index]!=0))/length(y_predict)
  
  
  #hinge loss
  hloss = hloss + 
    sum(hinge(unlist(lapply(X[test_index],function(x) con$slope(x)+con$intercept))*y[test_index]))
  
  #fitresult 
  fitresult[k,] = unlist(lapply(X,function(x) con$slope(x)+con$intercept))
  
  
  
  print(paste("cv",r,"-",cost,"-",k," is done",sep = ""))
  
  # get error for each CV
}
cvresult[1,1] = zloss/5
cvresult[1,2] = hloss/5


result = list()
result[[1]] = cvresult
result[[2]] = fitresult

save(result,file = paste("CV_",r,"_",method,"_",cost,".RData",sep = ""))
