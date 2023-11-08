library(doMC)
library(doParallel)
library(bigmemory)
library(ff)
Sys.setenv('R_MAX_VSIZE'=32000000000)
registerDoParallel(cores=8)
registerDoMC(cores = 8)
Rcpp::sourceCpp("src/qpadmslack.cpp")
N = 10000
p = 100
tau = 0.7
rho = 0.5
a = 3.7
pho = 5
beta_true = rep(0, p)
beta_true[6] = beta_true[12] = beta_true[15] = beta_true[20] = 1
gcov = function(p, rho){
  cov = matrix(1, p, p);
  for(i in 1:p){
    for(j in 1:p){
      if(i < j) cov[i,j] = rho^{j-i}
      else cov[i,j] = cov[j,i]
    }
  }
  cov
}
cov = gcov(p, rho)
set.seed(66)
X = matrix(rnorm(N*p), nrow=N)
X = X%*%chol(cov)
X[,1] = pnorm(X[,1])
e = rnorm(N)
Y = X[,6]+X[,12]+X[,15]+X[,20]+0.7*X[,1]*e
beta_true[1] = 0.7*qnorm(tau)
##we consider five different values for the partition number K, i.e., 1, 10, 20, 50 and 100
K = c(1, 10, 20, 50, 100)
##the corresponding suitable value for lambda under different K
lambda = c(20, 30, 60, 100, 200)
##AE and Time respectively record the estimation accuracy and computational time (pseudo-parallel) of the algorithm under different K
AE = Time = rep(0, length(K))
for(k in 1:length(K)){
  qpadmslack = paraQPADMslackcpp(Y, X, K[k], tau, "scad", a, lambda[k], pho)
  beta = qpadmslack$Estimation
  ##beta_true is the true value of the coefficient
  AE[k] = sum(abs(beta-beta_true))
  Time[k] = qpadmslack$Time
  print(k)
}
