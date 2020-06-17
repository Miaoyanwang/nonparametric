library(pracma)
library(quadprog)

eps = 10^-5



objv = function(B,b0,X,y,cost = 10,prob = F){
  if (prob == F) {
    value = sum(B*B)/2+cost*sum(pmax(1-y*unlist(lapply(X,function(x) sum(B*x)+b0)),0))
  }else{
    ind = which(y==1)
    value = sum(B*B)/2 +
      (1-prob)*cost*sum(pmax(1-y[ind]*unlist(lapply(X[ind],function(x) sum(B*x)+b0)),0)) +
      prob*cost*sum(pmax(1-y[-ind]*unlist(lapply(X[-ind],function(x) sum(B*x)+b0)),0))

  }
  return(value)
}

# Generating dataset
gendat = function(m,n,r,N,b0,type="non"){
  result = list()
  # simulation
  # Weight
  if(type=="non"){
  rU = matrix(runif(m*r,-1,1),nrow = m)
  rV = matrix(runif(n*r,-1,1),nrow = n)
  B = rU%*%t(rV)

  # predictor matrix
  X = list()
  for (i in 1:N) {
    X[[i]] <- matrix(runif(m*n,-1,1),nrow = m,ncol=n)
  }
  }else if(type=="sym"){
      rU = matrix(runif(m*r,-1,1),nrow = m)
      B = rU%*%t(rU)
      
      # predictor matrix
      X = list()
      for (i in 1:N) {
          X[[i]] <- matrix(runif(m*n,-1,1),nrow = m,ncol=n)
          X[[i]]=(X[[i]]+t(X[[i]]))/2
      }
  }

  # classification
  y = list()
  for (i in 1:N) {
    y[[i]] = sign(sum(B*X[[i]])+b0)
  }
  y = unlist(y)

  # predictor vector
  x = matrix(nrow =N,ncol = m*n)
  for(i in 1:N){
    x[i,] = as.vector(X[[i]])
  }
  dat = data.frame(y = factor(y), x)

  result$B = B
  result$X = X; result$y = y; result$dat = dat
  return(result)
}


kernelm = function(X,H,y,type = c("u","v")){
  n = length(X)
  x = matrix(unlist(X),nrow = length(X),byrow = T)
  if (type == "u") {
    hx = matrix(unlist(lapply(X,function(x) x%*%H)),nrow = length(X),byrow = T)
  } else {
    hx = matrix(unlist(lapply(X,function(x) H%*%x)),nrow = length(X),byrow = T)
  }
  Q = matrix(nrow = n,ncol = n)
  for (i in 1:n) {
    for(j in i:n){
      Q[i,j] = sum(x[i,]*hx[j,])*y[i]*y[j]
      Q[j,i] = Q[i,j]
    }
  }
  h = eigen(Q)
  Q = (h$vectors)%*%diag(pmax(h$values,eps))%*%t(h$vectors)
  return(Q)
}



# smm = function(X,y,r,cost = 10){
#   result = list()
#   error = 10
#   iter = 0
#   # SMM
#   m= nrow(X[[1]]); n = ncol(X[[1]]); N = length(X)
#
#   #initialization
#   U = randortho(m)[,1:r]
#   # U = matrix(runif(m*r,-1,1),nrow = m)
#   V = randortho(n)[,1:r]
#   # V = matrix(runif(n*r,-1,1),nrow = n)
#   obj = objv(U%*%t(V),0,X,y,cost);obj
#
#   while((iter <20)&(error>10^-4)){
#     # update U fixing V
#     Vs = V%*%solve(t(V)%*%V)
#     H = Vs%*%t(V)
#     dvec = rep(1,length(X))
#     Dmat = kernelm(X,H,y,"u")
#     Amat = cbind(y,diag(1,N),-diag(1,N))
#     bvec = c(rep(0,1+N),rep(-cost,N))
#     alpha = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
#     Bpart=matrix(t(y*alpha$solution)%*%matrix(unlist(X),nrow = length(X),byrow = T),nrow = m)
#     U = Bpart%*%Vs
#
#
#     # update V fixing U
#     Us = U%*%solve(t(U)%*%U)
#     H = Us%*%t(U)
#     Dmat = kernelm(X,H,y,"v")
#     alpha = solve.QP(Dmat,dvec,Amat,bvec,meq = 1)
#     Bpart=matrix(t(y*alpha$solution)%*%matrix(unlist(X),nrow = length(X),byrow = T),nrow = m)
#     V = t(Bpart)%*%Us
#
#
#     ## intercept estimation
#     Bhat = U%*%t(V);Bhat
#     positiv = min(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==1)])
#     negativ = max(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==-1)])
#     if ((1-positiv)<(-1-negativ)) {
#       b0hat = -(positiv+negativ)/2
#     }else{
#       gridb0 = seq(from = -1-negativ,to = 1-positiv,length = 100)
#       b0hat = gridb0[which.min(sapply(gridb0,function(b) objv(Bhat,b,X,y)))]
#     }
#     obj = c(obj,objv(Bhat,b0hat,X,y,cost));obj
#     iter = iter+1
#     error = abs(-obj[iter+1]+obj[iter])/obj[iter];error
#
#   }
#   predictor = function(x) sign(sum(Bhat*x)+b0hat)
#   result$B = Bhat; result$b0 = b0hat; result$obj = obj; result$iter = iter
#   result$error = error; result$predict = predictor
#   return(result)
# }




# svm = function(X,y,cost = 10){
#   result = list()
#   error = 10
#   iter = 0
#   # SVM
#   m= nrow(X[[1]]); n = ncol(X[[1]]); N = length(X)
#
#   H = diag(1,n)
#   dvec = rep(1,length(X))
#   Dmat = kernelm(X,H,y,"u")
#   Amat = cbind(y,diag(1,N),-diag(1,N))
#   bvec = c(rep(0,1+N),rep(-cost,N))
#   alpha = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
#   Bhat=matrix(t(y*alpha$solution)%*%matrix(unlist(X),nrow = length(X),byrow = T),nrow = m)
#   b0hat = -(min(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==1)])+
#               max(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==-1)]))/2
#   obj = objv(Bhat,b0hat,X,y,cost)
#
#   predictor = function(x) sign(sum(Bhat*x)+b0hat)
#   result$B = Bhat; result$b0 = b0hat; result$obj = obj;
#   result$predict = predictor
#   return(result)
# }


## SMM with multiple initialization and probability
smm = function(X,y,r,cost = 10,rep = 10, p = .5,drow){
  result = list()
  if (p==.5){
    cost = 2*cost
  }
  # SMM
  m= nrow(X[[1]]); n = ncol(X[[1]]); N = length(X)
  dcol=m-drow
  
  compareobj = 10^100
  for (nsim in 1:rep) {
    error = 10
    iter = 0
    #initialization
    U_temp = randortho(drow)[,1:r]
    # U = matrix(runif(m*r,-1,1),nrow = m)
    V_temp = randortho(dcol)[,1:r]
    U=matrix(0,nrow=m,ncol=2*r)
    U[1:drow,(r+1):(2*r)]=U_temp
    U[(drow+1):m,1:r]=V_temp
    
    V=matrix(0,nrow=m,ncol=2*r)
    V[1:drow,1:r]=U_temp
    V[(drow+1):m,(r+1):(2*r)]=V_temp
    
    # V = matrix(runif(n*r,-1,1),nrow = n)
    obj = objv(U%*%t(V),0,X,y,cost,prob = p);obj

    while((iter <20)&(error>10^-3)){
      # update U fixing V
      V=svd(V)$u
      Vs = V%*%solve(t(V)%*%V)
      H = Vs%*%t(V)
      dvec = rep(1,length(X))
      Dmat = kernelm(X,H,y,"u")
      Amat = cbind(y,diag(1,N),-diag(1,N))
      bvec = c(rep(0,1+N),ifelse(y==1,-cost*(1-p),-cost*p))
      alpha = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      Bpart=matrix(t(y*alpha$solution)%*%matrix(unlist(X),nrow = length(X),byrow = T),nrow = m)
      U = Bpart%*%Vs
      #val=c(val,alpha$value)
      ### U_new=matrix(0,nrow=5,ncol=2)
      #for(i in 1:N){
      #   U_new=U_new+alpha$solution[i]*y[i]*X[[i]]%*%V%*%solve(t(V)%*%V)
      #}

      # update V fixing U
      Us = U%*%solve(t(U)%*%U)
      H = Us%*%t(U)
      Dmat = kernelm(X,H,y,"v")
      alpha = solve.QP(Dmat,dvec,Amat,bvec,meq = 1)
      Bpart=matrix(t(y*alpha$solution)%*%matrix(unlist(X),nrow = length(X),byrow = T),nrow = m)
      V = t(Bpart)%*%Us
      #val=c(val,alpha$value)

      ## intercept estimation
      Bhat = U%*%t(V);Bhat
      positiv = min(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==1)])
      negativ = max(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==-1)])
      if ((1-positiv)<(-1-negativ)) {
        b0hat = -(positiv+negativ)/2
      }else{
        gridb0 = seq(from = -1-negativ,to = 1-positiv,length = 100)
        b0hat = gridb0[which.min(sapply(gridb0,function(b) objv(Bhat,b,X,y,cost,prob = p)))]
      }
      obj = c(obj,objv(Bhat,b0hat,X,y,cost,prob = p));obj
      iter = iter+1
      error = abs(-obj[iter+1]+obj[iter])/obj[iter];error

    }
    if (compareobj>obj[iter+1]) {
      compareobj = obj[iter+1]
      predictor = function(x) sign(sum(Bhat*x)+b0hat)
      result$B = Bhat; result$b0 = b0hat; result$obj = obj; result$iter = iter
      result$error = error; result$predict = predictor
    }

  }
  return(result)
}



kernelmat = function(x,y,kernels = function(x1,x2) sum(x1*x2)){
  N = length(y)
  Q = matrix(nrow = N,ncol = N)
  for (i in 1:N) {
    for(j in i:N){
      Q[i,j] =kernels(x[i,],x[j,])*y[i]*y[j]
      Q[j,i] = Q[i,j]
    }
  }
  h = eigen(Q)
  Q = (h$vectors)%*%diag(pmax(h$values,eps))%*%t(h$vectors)
  return(Q)
}


# SVM with kernel functions and weighted cost function
svm = function(X,y,cost = 10, kernels = function(x1,x2) sum(x1*x2), p = .5,minimization = FALSE){
  if (p==.5) {
    cost = 2*cost
  }
  result = list()
  error = 10
  iter = 0
  # SVM
  m= nrow(X[[1]]); n = ncol(X[[1]]); N = length(X)

  x = matrix(unlist(X),nrow = N,byrow = T)
  dvec = rep(1,length(X))
  Dmat = kernelmat(x,y,kernels)
  Amat = cbind(y,diag(1,N),-diag(1,N))
  bvec = c(rep(0,1+N),ifelse(y==1,-cost*(1-p),-cost*p))
  alpha = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
  coef = y*alpha$solution

  Bhat=matrix(t(coef)%*%x,nrow = m)
  # b0hat = -(min(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==1)])+
  #             max(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==-1)]))/2
  b0hat = -(min((Dmat%*%alpha$solution)[which(y==1)]*y[which(y==1)])+
              max((Dmat%*%alpha$solution)[which(y==-1)]*y[which(y==-1)]))/2
  if(minimization==TRUE){
    positiv = min(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==1)])
    negativ = max(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==-1)])
    if ((1-positiv)<(-1-negativ)) {
      b0hat = -(positiv+negativ)/2
    }else{
      gridb0 = seq(from = -1-negativ,to = 1-positiv,length = 100)
      b0hat = gridb0[which.min(sapply(gridb0,function(b) objv(Bhat,b,X,y,prob = p)))]
    }
  }
  obj = objv(Bhat,b0hat,X,y,cost,prob = p)

  predictor = function(y){
    return(sign(sum(coef*apply(x,1, function(x) kernels(x,y)))+b0hat))
  }
  # predictor = function(x) sign(sum(Bhat*x)+b0hat)
  result$B = Bhat; result$b0 = b0hat; result$obj = obj;
  result$predict = predictor
  return(result)
}


posteriorSMM = function(X,y,r,test,cost = 10,rep = 10, precision = 0.1){
  a = 1:(1/precision-1)
  for(i in 1:(1/precision-1)){
    fit = smm(X,y,r,cost,rep, p = i*precision)$predict
    a[i] = fit(test)
  }
  if (all(a==1)) {
    return(1)
  }else if(all(a==-1)){
    return(0)
  }else{
    return((max(which(a==1))+min(which(a==-1)))/(2/precision))
  }
}

posterior = function(X,y,cost = 10,test,precision=0.1,kernels = function(x1,x2) sum(x1*x2)){
  a = 1:(1/precision-1)
  for(i in 1:(1/precision-1)){
    fit = svm(X,y,cost, kernels, p = i*precision)$predict
    a[i] = fit(test)
  }
  if (all(a==1)) {
    return(1)
  }else if(all(a==-1)){
    return(0)
  }else{
    return((max(which(a==1))+min(which(a==-1)))/(2/precision))
  }
}

posterior2 = function(dat,test,precision=0.1){
  a = 1:(1/precision-1)
  for(i in 1:(1/precision-1)){
    classcosts <- table(as.factor(dat$y))  # the weight vector must be named with the classes names
    classcosts[1] <- i# a class -1 mismatch has a terrible cost
    classcosts[2] <- (1/precision) - i   # a class +1 mismatch not so much...

    fit = e1071::svm(factor(y) ~ ., data = dat, scale = FALSE, kernel = "radial",
                     class.weights = classcosts)
    a[i] = ifelse(predict(fit, test)==1,1,-1)
  }
  if (all(a==1)) {
    return(1)
  }else if(all(a==-1)){
    return(0)
  }else{
    return((max(which(a==1))+min(which(a==-1)))/(2*(1/precision)))
  }
}
