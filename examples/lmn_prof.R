# generate data
n <- 50
q <- 2
Y <- matrix(rnorm(n*q),n,q) # response matrix
X <- matrix(1,n,1) # covariate matrix
V <- exp(-(1:n)/n) # diagonal variance specification
suff <- lmn_suff(Y, X = X, V = V, Vtype = "diag") # sufficient statistics

# profile loglikelihood
lmn_prof(suff)

# check that it's the same as loglikelihood at MLE
lmn_loglik(Beta = suff$Bhat, Sigma = suff$S/suff$n, suff = suff)
