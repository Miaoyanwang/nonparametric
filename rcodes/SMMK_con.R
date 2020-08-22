#SMMK
library(pracma)
Makepositive = function(mat){
  h = eigen(mat,symmetric = T)
  nmat = (h$vectors)%*%diag(pmax(h$values,10^-4))%*%t(h$vectors)
  return(nmat)
}

Karray2 = function(X, kernel = function(a,b) sum(a*b),type = c("row","col")){
  d1= nrow(X[[1]]); d2 = ncol(X[[1]]); n = length(X)
  if(type == "row"){
    KX = array(dim = c(n,n,d1,d1))
    for(i in 1:n){
      for(j in 1:n){
        KX[i,j, , ] = K(X[[i]],X[[j]],kernel,type)
      }
    }
  }else{
    KX = array(dim = c(n,n,d2,d2))
    for(i in 1:n){
      for(j in 1:n){
        KX[i,j, , ] = K(X[[i]],X[[j]],kernel,type)
      }
    }
  }
  return(KX)
}


K = function(A,B,kernel = function(a,b) sum(a*b),type = c("row","col")){
  if(nrow(A) != nrow(B)) stop("two matrix sizes are different")
  if(ncol(A) != ncol(B)) stop("two matrix sizes are different")
  d1 = nrow(A); d2 = ncol(A)
  if (type == "row") {
    Kab = matrix(nrow= d1,ncol = d1)
    for(i in 1:d1){
      for(j in 1:d1){
        Kab[i,j] = kernel(A[i,],B[j,])
      }
    }
  }else{
    Kab = matrix(nrow= d2,ncol = d2)
    for(i in 1:d2){
      for(j in 1:d2){
        Kab[i,j] = kernel(A[,i],B[,j])
      }
    }
  }
  return(Kab)
}


SMMK_con = function(X,y,r,kernel = function(a,b) sum(a*b), cost = 10, rep = 1, p = .5){
  result = list()
  
  d1 = nrow(X[[1]]); d2 = ncol(X[[1]]); n = length(X)
  K_row = Karray2(X,kernel,type = "row")
  K_col = Karray2(X,kernel,type = "col")
  compareobj = 10^10
  for(nsim in 1:rep){
    error = 10; iter = 0; obj = 10^10
    # initialize P_row,P_col
    P_row = randortho(d1)[,1:r,drop = F]; P_col = randortho(d2)[,1:r,drop = F]  
    
    
    
    while((iter < 20)&(error >10^-5)){
      # update C
      W_row = P_row%*%t(P_row); W_col = P_col%*%t(P_col)
      Dmat = matrix(nrow = n,ncol = n)
      for(i in 1:n){
        for(j in 1:n){ 
          Dmat[i,j] = sum(W_row*K_row[i,j,,])+sum(W_col*K_col[i,j,,])
        }
      }  
      dvec = rep(1,n)
      Dmat = Makepositive((y%*%t(y))*Dmat)
      Amat = cbind(y,diag(1,n),-diag(1,n))
      bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
      res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      alpha=res$solution
      
      
      # update P
      CC_row = 0; CPh_row = array(0,dim = c(n,r,d1))
      CC_col = 0; CPh_col = array(0,dim = c(n,r,d2))
      for(i in 1:n){
        for(j in 1:n){
          CC_row = CC_row + y[i]*alpha[i]*y[j]*alpha[j]*K_row[i,j,,]
          CC_col = CC_col + y[i]*alpha[i]*y[j]*alpha[j]*K_col[i,j,,]
        }
      }
      CCi_row = solve(t(P_row)%*%CC_row%*%P_row)
      CCi_col = solve(t(P_col)%*%CC_col%*%P_col)
      
      for(i in 1:n){
        cph_row = 0; cph_col = 0;
        for(j in 1:n){
          cph_row = cph_row + y[j]*alpha[j]*K_row[j,i,,]
          cph_col = cph_col + y[j]*alpha[j]*K_col[j,i,,]
        }
        CPh_row[i,,] = t(P_row)%*%cph_row
        CPh_col[i,,] = t(P_col)%*%cph_col
      }
      
      Dmat = matrix(nrow = n,ncol = n)
      for(i in 1:n){
        for(j in 1:n){ 
          Dmat[i,j] = sum(CCi_row%*%CPh_row[i,,]*CPh_row[j,,])+
                                   sum(CCi_col%*%CPh_col[i,,]*CPh_col[j,,])
        }
      }  
      
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
      P_row = 0; P_col = 0
      for(i in 1:n){
        P_row = P_row + alpha[i]*y[i]*rbind(CPh_row[i,,])
        P_col = P_col + alpha[i]*y[i]*rbind(CPh_col[i,,])
      }
      P_row = svd(t(P_row)%*%CCi_row)$u; P_col = svd(t(P_col)%*%CCi_col)$u
      
    }
    if(compareobj>obj[iter+1]){
      W_row = P_row%*%t(P_row); W_col = P_col%*%t(P_col)
      Dmat = matrix(nrow = n,ncol = n)
      for(i in 1:n){
        for(j in 1:n){ 
          Dmat[i,j] = sum(W_row*K_row[i,j,,])+sum(W_col*K_col[i,j,,])
        }
      }  
      dvec = rep(1,n)
      Dmat = Makepositive((y%*%t(y))*Dmat)
      Amat = cbind(y,diag(1,n),-diag(1,n))
      bvec = c(rep(0,1+n),ifelse(y==1,-cost*(1-p),-cost*p))
      res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      alpha=res$solution
      
      fx = function(Xnew){
        
        value = t(as.matrix(alpha*y))%*%
          unlist(lapply(X,
                        function(x) sum(W_row*K(x,Xnew,kernel,type = "row"))+
                          sum(W_col*K(x,Xnew,kernel,type = "col"))))
        return(value)
      }
      # intercept part estimation (update b)
      positiv = min(unlist(lapply(X,fx))[which(y==1)])
      negativ = max(unlist(lapply(X,fx))[which(y==-1)])
      intercept = -(positiv+negativ)/2
      
      compareobj = obj[iter+1]
      dfunc = function(Xnew) fx(Xnew)+intercept ; classifier = function(Xnew) sign(fx(Xnew)+intercept)
      result$alpha = alpha
      result$dfunc = dfunc; result$classifier = classifier
      result$P_row = P_row; result$P_col = P_col
      result$obj = obj[-1]; result$iter = iter; result$error = error
    }
  }
  return(result)
}

