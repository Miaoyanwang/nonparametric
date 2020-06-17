lwpsvm <- function(x, y, H, lambda) 
{
require(kernlab)
n <- length(y)
p <- ncol(x)

step <- 1/H
pi.grid <- seq(step, 1-step, by = step)

bar.x <- apply(x, 2, mean)
cov.x <- cov(x)

# standardization
temp <- eigen(cov.x)
D <- diag(sqrt(temp$values))
V <- temp$vectors
sd.x <-  V %*% D
inv.sd.x <- diag(1/sqrt(temp$values)) %*% t(V)
z <- t(inv.sd.x %*% (t(x) - bar.x))

    w <- matrix(0, p, H)
    for (h in 1:(H-1)) {
    alpha <- rep(0, n)
    temp <- ksvm(x = z, y = as.factor(y), type = "C-svc", kernel = "vanilladot", kpar = list(), C = lambda, class.weights = c("1" = 1-pi.grid[h], "-1" = pi.grid[h]))
    alpha[temp@SVindex] <- unlist(temp@alpha)
    w[,h] <- apply(unlist(temp@coef) * z[temp@SVindex,], 2, sum)
    }

psi <- solve(t(sd.x)) %*% w
Mn <- matrix(0, p, p)
for (h in 1:H) Mn <- Mn + psi[,h, drop = F] %*% t(psi[,h, drop = F])

obj <- eigen(Mn)
obj
}
