# generate data
n <- 50
q <- 2
p <- 3
Y <- matrix(rnorm(n*q),n,q) # response matrix
X <- matrix(rnorm(n*p),n,p) # covariate matrix
V <- .5 * exp(-(1:n)/n) # Toeplitz variance specification

suff <- lmn_suff(Y = Y, X = X, V = V, Vtype = "acf") # sufficient statistics
# default noninformative prior pi(Beta, Sigma) ~ |Sigma|^(-(q+1)/2)
prior <- lmn_prior(p = suff$p, q = suff$q)
post <- lmn_post(suff, prior = prior) # posterior MNIW parameters
lmn_marg(suff, prior = prior, post = post)
