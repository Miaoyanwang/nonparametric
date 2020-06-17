Gn <- function(rho, obj, n){
p <- nrow(obj$vectors)
v <- obj$values

value <- NULL
for (k in 1:p) {
value[k] <- sum(v[1:k]) - rho * (k * log(n) / n^(1/2)) * v[1]} 
value
}
