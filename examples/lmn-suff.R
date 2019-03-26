# Data
n <- 50
q <- 3
Y <- matrix(rnorm(n*q),n,q)

# No intercept, diagonal V input
X <- 0
V <- exp(-(1:n)/n)
lmn.suff(Y, X = X, V = V, Vtype = "diag")

# X = (scaled) Intercept, scalar V input (no need to specify Vtype)
X <- 2
V <- .5
lmn.suff(Y, X = X, V = V)

# X = dense matrix, Toeplitz variance matrix
p <- 2
X <- matrix(rnorm(n*p), n, p)
acf <- SuperGauss::Toeplitz(acf = 0.5*exp(-seq(1:n)/n))
lmn.suff(Y, X = X, V = acf, Vtype = "acf")

