# generate data
n <- 50
q <- 3
Y <- matrix(rnorm(n*q),n,q) # response matrix
X <- 1 # intercept covariate
V <- 0.5 # scalar variance specification
suff <- lmn_suff(Y, X = X, V = V) # sufficient statistics

# calculate loglikelihood
Beta <- matrix(rnorm(q),1,q)
Sigma <- diag(rexp(q))
lmn_loglik(Beta = Beta, Sigma = Sigma, suff = suff)
