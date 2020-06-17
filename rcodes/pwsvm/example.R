rm(list = ls())
library(kernlab)

setwd("~/Desktop/pwsvm_R")

sourceDir <- function(path, trace = TRUE, ...)
{
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$"))
  {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir("fn")

wisc <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data", 
                   sep = ",")

x <- matrix(unlist(wisc[,-c(1,2)]), ncol = 30)
y <- 2*as.numeric(unlist(wisc[,2])) - 3
n <- length(y)
p <- ncol(x)
lambda <- 1
H <- 20



################
# linear PWSVM #
################

obj <- temp.lsvm <- lwpsvm(x, y, H, lambda) 
lsvm <- temp.lsvm$vectors
value.lsvm <- temp.lsvm$values
x.lsvm <- x %*% lsvm

plot(x.lsvm[,1], x.lsvm[,2], type = "n", xlab = "1st predictor", ylab  = "2nd predictor")
points(x.lsvm[y == 1,1], x.lsvm[y == 1,2], col = 4, pch = "+")
points(x.lsvm[y != 1,1], x.lsvm[y != 1,2], col = 2)

###################
# nonlinear PWSVM #
###################

temp.nsvm <- nwpsvm(x, y, lambda, prop.d = 1/5)
x.nsvm <- phix(x, temp.nsvm)

boxplot(x.nsvm[y == 1,1], x.nsvm[y != 1,1], 
        xlab = "Y", axes = F, 
        ylab = expression(hat(phi)[1](x)))
axis(1, seq(0.5, 2.5, by = 0.5), c(NA, "+1", NA, "-1", NA))
axis(2, las = 1)
