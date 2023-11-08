#This file is largely the same as the currentExample file
#It defines a bunch of values and constructs the X and Y matrixes.
#The major difference is that it will output X and Y as csv files
#The files are named X and Y, respectively.
#The purpose was to output this for testing the Python scripts
N = 100000
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

write.table(X, "X", append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = FALSE)

write.table(Y, "Y", append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = FALSE)
