#-------------------------------------------------------------------------------
#
# test functions for LMN package
#
#-------------------------------------------------------------------------------

source("lmn-package.R")
require(mvtnorm)

#-------------------------------------------------------------------------------
#
# test lmn.suff
#
#-------------------------------------------------------------------------------

# using full matrix V

n <- 10
q <- 3
p <- 2
npred <- 5

Y <- rMnorm(n, q)
XX <- rMnorm(n+npred, p)
VV <- crossprod(rMnorm(n+npred))
X <- XX[1:n,]
x <- XX[n+1:npred,]
V <- VV[1:n,1:n]
w <- VV[1:n,n+1:npred]
v <- VV[n+1:npred,n+1:npred]

suff <- lmn.suff(Y = Y, X = XX, V = VV, npred = npred, debug = FALSE)

Beta.hat <- solve(crossprod(X, solve(V, X))) %*% crossprod(X, solve(V, Y))
range(Beta.hat - suff$Beta.hat)

T <- crossprod(X, solve(V, X))
range(T - suff$T)

S <- crossprod((Y - X%*%Beta.hat), solve(V, Y - X%*%Beta.hat))
range(S - suff$S)

ldV <- ldet(V)
range(ldV - suff$ldV)

Ap <- crossprod(w, solve(V, Y))
range(Ap - suff$Ap)

Xp <- x - crossprod(w, solve(V, X))
range(Xp - suff$Xp)

Vp <- v - crossprod(w, solve(V, w))
range(Vp - suff$Vp)


# using diagonal matrix diagV

n <- 10
q <- 3
p <- 2
npred <- 5

Y <- rMnorm(n, q)
XX <- rMnorm(n+npred, p)
diagV <- rchisq(1, df = 2)
VV <- diag(rep(diagV, len = n+npred))
X <- XX[1:n,]
x <- XX[n+1:npred,]
V <- VV[1:n,1:n]
w <- VV[1:n,n+1:npred]
v <- VV[n+1:npred,n+1:npred]

suff <- lmn.suff(Y = Y, X = XX, diagV = diagV, npred = npred, debug = FALSE)

Beta.hat <- solve(crossprod(X, solve(V, X))) %*% crossprod(X, solve(V, Y))
range(Beta.hat - suff$Beta.hat)

T <- crossprod(X, solve(V, X))
range(T - suff$T)

S <- crossprod((Y - X%*%Beta.hat), solve(V, Y - X%*%Beta.hat))
range(S - suff$S)

ldV <- ldet(V)
range(ldV - suff$ldV)

Ap <- crossprod(w, solve(V, Y))
range(Ap - suff$Ap)

Xp <- x - crossprod(w, solve(V, X))
range(Xp - suff$Xp)

Vp <- v - crossprod(w, solve(V, w))
range(Vp - suff$Vp)

#-------------------------------------------------------------------------------
#
# test lmn.loglik
#
#-------------------------------------------------------------------------------

llMNorm <- function(Y, Mu, V, Sigma) {
  n <- nrow(Y)
  p <- ncol(Y)
  Z <- Y - Mu
  Z <- solve(Sigma, t(Z)) %*% solve(V, Z)
  lp <- n*p*log(2*pi) + sum(diag(Z))
  lp <- lp + n * ldet(Sigma) + p * ldet(V)
  -.5 * lp
}


n <- 10
q <- 3
p <- 2

Y <- rMnorm(n, q)
X <- rMnorm(n, p)
V <- crossprod(rMnorm(n))

suff <- lmn.suff(Y = Y, X = X, V = V)

Beta <- rMnorm(p,q)
Sigma <- crossprod(rMnorm(q))
ll1 <- lmn.loglik(Beta = Beta, Sigma = Sigma, suff = suff)
ll2 <- dmvnorm(c(Y), c(X %*% Beta), kronecker(Sigma, V), log = TRUE)
ll3 <- llMNorm(Y, X %*% Beta, V, Sigma)
ll1-ll2


#-------------------------------------------------------------------------------
#
# test lmn.lprof
#
#-------------------------------------------------------------------------------

n <- 10
q <- 3
p <- 2

Y <- rMnorm(n, q)
X <- rMnorm(n, p)
V <- crossprod(rMnorm(n))
suff <- lmn.suff(Y = Y, X = X, V = V)
Beta.hat <- suff$Beta.hat
Sigma.hat <- suff$S/suff$n
ll1 <- lmn.loglik(Beta.hat, Sigma.hat, suff = suff)
ll2 <- lmn.lprof(suff)
ll1-ll2
