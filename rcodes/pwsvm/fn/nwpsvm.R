psi.function <- function(value, x, v, w, l, kernel.function, kernel.param){
value <- matrix(value, 1, length(value))
temp <- kernel.function(value, x, kernel.param)
psi.value <- apply(w * c(temp - mean(temp)), 2, sum)/l
psi.value %*% v
}


kernel.function <- function (x, y = x, param.kernel = 1/p) 
{
    n <- nrow(x)
    m <- nrow(y)
    p <- ncol(x)
    normx <- drop((x^2) %*% rep(1, p))
    normy <- drop((y^2) %*% rep(1, p))
    a <- x %*% t(y)
    a <- (-2 * a + normx) + outer(rep(1, n), normy, "*")
    exp(-a * param.kernel)
}



nwpsvm <- function(x, y, lambda, H = 10, prop.d = 1/2) {
require(kernlab)
n <- length(y)
tau <- mean(as.numeric(dist(x)))
kernel.param <- 1/tau^2

Kn <- kernel.function(x, x, kernel.param)
Qn <- diag(n) - matrix(1/n, n, n)

temp <- Qn %*% Kn %*% Qn
kk <- floor(prop.d * n)

e.result <- eigen(temp)
w <- psi <- e.result$vectors[,1:kk, drop = F]
l <- e.result$values[1:kk]

P.psi <- psi %*% solve(t(psi) %*% psi) %*% t(psi)

P.script <- P.psi * outer(y, y) 

pi.grid <- seq(0, 1, length = H+2)[-c(1,H+2)]
Mn <- matrix(0, kk, kk)
h <- 1
for (h in 1:H) {
    u <- lambda * n * matrix(wvec(pi.grid[h], y))
    qp <- ipop(c = -matrix(rep(1, n)), H = (1/2) * P.script, 
               A = y, b = 0, r = 0, l = matrix(0, n), u = u)
    alpha <- qp@primal
    cvec <- 1/2 * solve(t(psi) %*% psi) %*% t(psi) %*% (alpha * y)
    Mn <- Mn + cvec %*% t(cvec)
    }

    result <- eigen(Mn)
    v <- result$vectors
    obj <- list(psi.function = psi.function, x = x, v = v, w = w, l = l, 
                kernel.function = kernel.function, kernel.param = kernel.param)
    obj
}


phix <- function(value, obj, order = 2) {
    psi.function <- obj$psi.function
    x <- obj$x
    v <- obj$v
    w <- obj$w
    l <- obj$l
    kernel.function <- obj$kernel.function
    kernel.param <- obj$kernel.param
    
    p <- ncol(x)
    if (length(value) == p) {
    temp <- psi.function(value, x, v[,1:order, drop = F], w, l, kernel.function, kernel.param)
    } else if (ncol(value) == p) {
    temp <- t(apply(value, 1, psi.function, x, v[,1:order, drop = F], w, l, kernel.function, kernel.param))
    } else if (nrow(value) == p) {
    temp <- t(apply(value, 2, psi.function, x, v[,1:order, drop = F], w, l, kernel.function, kernel.param))
    } else stop("check `str(value)`")
    temp
}
