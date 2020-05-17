#SMMK
library(pracma)
library(quadprog)
# Make sure Q matrix is positive definite
Makepositive = function(mat){
  h = eigen(mat,symmetric = T)
  nmat = (h$vectors)%*%diag(pmax(h$values,10^-4))%*%t(h$vectors)
  return(nmat)
}

# Make symmnetric matrix [0,X^T;X,0]
makesym = function(mat){
  m = nrow(mat); n = ncol(mat)
  nmat = rbind(cbind(matrix(0,n,n),t(mat)),cbind(mat,matrix(0,m,m)))
  return(nmat)
}

# save kernel values
Karray = function(X,kernel = function(X1,X2) t(X1)%*%X2){
  m= nrow(X[[1]]); n = ncol(X[[1]]); N = length(X)
  K = array(dim = c(N,N,n,n))
  for(i in 1:N){
    for(j in 1:N){
      K[i,j, , ] = kernel(X[[i]],X[[j]])
    }
  }
  return(K)
}


# objective value function
objm = function(X,y,alpha,V,b,K,cost = 10){
  m= nrow(X[[1]]); n = ncol(X[[1]]); N = length(X)
  Hv = V%*%solve(t(V)%*%V)%*%t(V)
  Kv = matrix(nrow =N,ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      Kv[i,j] = sum(Hv*K[i,j,,])
    }
  }
  coef = as.matrix(alpha*y)
  obj = t(coef)%*%Kv%*%coef/2 + cost*sum(pmax(1-y*(Kv%*%coef+b),0))
  return(obj)
}



# Main function (not available for weighed loss yet)
smmk = function(X,y,r,kernel = function(X1,X2) t(X1)%*%X2, cost = 10,rep = 1,p = .5){
  result = list()
  m= nrow(X[[1]]); n = ncol(X[[1]]); N = length(X)
  K = Karray(X,kernel)

  compareobj = 10^10
  for(nsim in 1:rep){
    error = 10; iter = 0; V = randortho(n)[,1:r,drop = F]
    ##obj = compareobj


    ### update U fixing V  ###
    ## normalize V
    V=svd(V)$u
    ## vtvi = solve(t(V)%*%V) ## commented out by Miaoyan
    Hv = V%*%t(V)
    Kv = matrix(nrow =N,ncol = N)
    for(i in 1:N){
      for(j in 1:N){
        Kv[i,j] = sum(Hv*K[i,j,,])
      }
    }
    dvec = rep(1,length(X))
    Dmat = Makepositive((y%*%t(y))*Kv)
    Amat = cbind(y,diag(1,N),-diag(1,N))
    bvec = c(rep(0,1+N),ifelse(y==1,-cost*(1-p),-cost*p))
    res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
    alpha=res$solution
    obj=-res$value

    while((iter < 20)&(error >10^-3)){
      ### update V fixing U ###

      # sum_{i,j} alpha_ialpha_jy_iy_jK(X_i,X_j)
      aayyK = 0
      for(i in 1:N){
        for(j in 1:N){
          aayyK = aayyK + alpha[i]*alpha[j]*y[i]*y[j]*K[i,j,,]
        }
      }

      # sum_j alpha_iy_iK(X_i,X_j)
      ayK = array(0,dim = c(N,n,n))
      for(i in 1:N){
        for(j in 1:N){
          ayK[i,,] = ayK[i,,] + alpha[j]*y[j]*K[j,i,,]
        }
      }

      # (U^TU)^{-1}
      utui = solve(t(V)%*%aayyK%*%V)

      # U^th(X_i)
      uth = array(dim = c(N,r,n))
      for(i in 1:N){
        uth[i,,] = t(V)%*%ayK[i,,]
      }

      # Ku[i,j] = tr(H_uh(X_i),H_uh(X_j))
      Ku = matrix(nrow = N,ncol = N)
      for(i in 1:N){
        for(j in 1:N){
          Ku[i,j] = sum(uth[i,,]*(utui%*%uth[j,,]))
        }
      }
      dvec = rep(1,length(X))
      Dmat = Makepositive((y%*%t(y))*Ku)
      Amat = cbind(y,diag(1,N),-diag(1,N))
      bvec = c(rep(0,1+N),ifelse(y==1,-cost*(1-p),-cost*p))
      beta = solve.QP(Dmat,dvec,Amat,bvec,meq =1)$solution

      # update V
      V = 0
      for(i in 1:N){
        V = V+rbind(uth[i,,])*beta[i]*y[i]
      }
      V = t(V)%*%utui


      ### update U fixing V  ###
      V=svd(V)$u
      #vtvi = solve(t(V)%*%V)
      Hv = V%*%t(V)
      Kv = matrix(nrow =N,ncol = N)
      for(i in 1:N){
        for(j in 1:N){
          Kv[i,j] = sum(Hv*K[i,j,,])
        }
      }
      dvec = rep(1,length(X))
      Dmat = Makepositive((y%*%t(y))*Kv)
      Amat = cbind(y,diag(1,N),-diag(1,N))
      bvec = c(rep(0,1+N),ifelse(y==1,-cost*(1-p),-cost*p))
      res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
      alpha=res$solution

      obj = c(obj,-res$value)
      iter = iter+1
      error = abs(-obj[iter+1]+obj[iter])/obj[iter]
    }

    # slope part estimation
    if (compareobj>obj[iter+1]) {
      slope = function(nX){
        coef <-  as.matrix(alpha*y)
        sp <-  t(coef)%*%unlist(lapply(X,function(x) sum(Hv*kernel(nX,x))))
        return(sp)
      }

      # intercept part estimation (update b)
      positiv = min(unlist(lapply(X,slope))[which(y==1)])
      negativ = max(unlist(lapply(X,slope))[which(y==-1)])
      if ((1-positiv)<(-1-negativ)) {
        b0hat = -(positiv+negativ)/2
      }else{
        gridb0 = seq(from = -1-negativ,to = 1-positiv,length = 100)
        b0hat = gridb0[which.min(sapply(gridb0,function(b) objm(X,y,alpha,V,b,K,cost)))]
      }
      ##obj = c(obj,objm(X,y,alpha,V,b0hat,K,cost));obj
      #iter = iter+1
      #error = abs(-obj[iter+1]+obj[iter])/obj[iter];error
      compareobj = obj[iter+1]
      predictor = function(nX) sign(slope(nX)+b0hat)
      result$slope = slope; result$b0 = b0hat; result$obj = obj[-1]; result$iter = iter
      result$error = error; result$predict = predictor; result$V = V
    }


  }
  return(result)
}


## Some kernels (Expkernel does not work)
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


