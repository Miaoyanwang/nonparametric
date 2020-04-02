library(pracma)
library(quadprog)

eps = 10^-5

sqrtH = function(Us,U){
  h = eigen(Us%*%t(U))
  return(h$vectors%*%diag(sqrt(pmax(h$values,eps))))
}

objv = function(B,b0,X,y,cost = 10){
  return(sum(B*B)/2+cost*sum(pmax(1-y*unlist(lapply(X,function(x) sum(B*x)+b0)),0)))
}

# Generating dataset
gendat = function(m,n,r,N,b0){
  result = list()
  # simulation
  # Weight
  rU = matrix(runif(m*r,-1,1),nrow = m)
  rV = matrix(runif(n*r,-1,1),nrow = n)
  B = rU%*%t(rV)

  # predictor matrix
  X = list()
  for (i in 1:N) {
    X[[i]] <- matrix(runif(m*n,-1,1),nrow = m,ncol=n)
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



smm = function(X,y,r,cost = 10){
  result = list()
  error = 10
  iter = 0
  # SMM
  m= nrow(X[[1]]); n = ncol(X[[1]]); N = length(X)

  #initialization
  U = randortho(m)[,1:r]
  # U = matrix(runif(m*r,-1,1),nrow = m)
  V = randortho(n)[,1:r]
  # V = matrix(runif(n*r,-1,1),nrow = n)
  obj = objv(U%*%t(V),0,X,y,cost);obj

  while((iter <20)&(error>10^-4)){
    # update U fixing V
    Vs = V%*%solve(t(V)%*%V)
    H = Vs%*%t(V)
    dvec = rep(1,length(X))
    Dmat = kernelm(X,H,y,"u")
    Amat = cbind(y,diag(1,N),-diag(1,N))
    bvec = c(rep(0,1+N),rep(-cost,N))
    alpha = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
    Bpart=matrix(t(y*alpha$solution)%*%matrix(unlist(X),nrow = length(X),byrow = T),nrow = m)
    U = Bpart%*%Vs


    # update V fixing U
    Us = U%*%solve(t(U)%*%U)
    H = Us%*%t(U)
    Dmat = kernelm(X,H,y,"v")
    alpha = solve.QP(Dmat,dvec,Amat,bvec,meq = 1)
    Bpart=matrix(t(y*alpha$solution)%*%matrix(unlist(X),nrow = length(X),byrow = T),nrow = m)
    V = t(Bpart)%*%Us


    ## intercept estimation
    Bhat = U%*%t(V);Bhat
    positiv = min(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==1)])
    negativ = max(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==-1)])
    if ((1-positiv)<(-1-negativ)) {
      b0hat = -(positiv+negativ)/2
    }else{
      gridb0 = seq(from = -1-negativ,to = 1-positiv,length = 100)
      b0hat = gridb0[which.min(sapply(gridb0,function(b) objv(Bhat,b,X,y)))]
    }
    obj = c(obj,objv(Bhat,b0hat,X,y,cost));obj
    iter = iter+1
    error = abs(-obj[iter+1]+obj[iter])/obj[iter];error

  }
  predictor = function(x) sign(sum(Bhat*x)+b0hat)
  result$B = Bhat; result$b0 = b0hat; result$obj = obj; result$iter = iter
  result$error = error; result$predict = predictor
  return(result)
}







svm = function(X,y,cost = 10){
  result = list()
  error = 10
  iter = 0
  # SVM
  m= nrow(X[[1]]); n = ncol(X[[1]]); N = length(X)

  H = diag(1,n)
  dvec = rep(1,length(X))
  Dmat = kernelm(X,H,y,"u")
  Amat = cbind(y,diag(1,N),-diag(1,N))
  bvec = c(rep(0,1+N),rep(-cost,N))
  alpha = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
  Bhat=matrix(t(y*alpha$solution)%*%matrix(unlist(X),nrow = length(X),byrow = T),nrow = m)
  b0hat = -(min(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==1)])+
              max(unlist(lapply(X,function(x) sum(Bhat*x)))[which(y==-1)]))/2
  obj = objv(Bhat,b0hat,X,y,cost)

  predictor = function(x) sign(sum(Bhat*x)+b0hat)
  result$B = Bhat; result$b0 = b0hat; result$obj = obj;
  result$predict = predictor
  return(result)
}



