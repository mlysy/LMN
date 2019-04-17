#--- code for vignette ----------------------------------------------------

require(LMN)
require(mniw) # won't be required later

#--- toy example ----------------------------------------------------------

# problem dimensions
qq <- 2 # number of responses
n <- 200 # number of observations

# parameters
Beta <- matrix(c(.3, .5, .7, .2), 2, qq)
Sigma <- matrix(c(.005, -.001, -.001, .002), qq, qq)
alpha <- .4
lambda <- .1

# simulate data
xseq <- seq(0, 10, len = n) # x vector
X <- cbind(1, xseq^alpha) # covariate matrix
V <- exp(-(outer(xseq, xseq, "-")/lambda)^2) # between-row variance
Mu <- X %*% Beta # mean matrix
Y <- mniw::rMNorm(n = 1, Mu = Mu, RowV = V, ColV = Sigma) # response matrix

# plot data
par(mfrow = c(1,1))
plot(x = 0, type = "n", xlim = range(xseq), ylim = range(Mu, Y),
     xlab = "x", ylab = "y")
lines(x = xseq, y = Mu[,1], col = "red")
lines(x = xseq, y = Mu[,2], col = "blue")
points(x = xseq, y = Y[,1], pch = 16, col = "red")
points(x = xseq, y = Y[,2], pch = 16, col = "blue")
legend("bottomright",
       legend = c("Response 1", "Response 2", "Expected", "Observed"),
       lty = c(NA, NA, 1, NA), pch = c(22, 22, NA, 16),
       pt.bg = c("red", "blue", "black", "black"), seg.len = 1)

# MLE of (Beta, Sigma) for theta = true value
suff <- lmn.suff(Y = Y, X = X, V = V, Vtype = "full")
Beta_mle <- suff$Bhat
Sigma_mle <- suff$S/suff$n

cov2sigrho <- function(Sigma) {
  sig <- sqrt(diag(Sigma))
  n <- length(sig) # dimensions of Sigma
  names(sig) <- paste0("sigma",1:n)
  # indices of upper triangular elements
  iupper <- matrix(1:(n^2),n,n)[upper.tri(Sigma, diag = FALSE)]
  rho <- cov2cor(Sigma)[iupper]
  rnames <- apply(expand.grid(1:n, 1:n), 1, paste0, collapse = "")
  names(rho) <- paste0("rho",rnames[iupper])
  c(sig,rho)
}

# estimate theta = (alpha, lambda)
toy_acf <- function(lambda) exp(-((xseq-xseq[1])/lambda)^2)

# check that it's indeed Toeplitz
all.equal(suff,
          lmn.suff(Y = Y, X = X, V = toy_acf(lambda), Vtype = "acf"))

# check that it's faster to use Vtype = "acf"
system.time({
  replicate(100, lmn.suff(Y = Y, X = X, V = V, Vtype = "full"))
})
system.time({
  replicate(100, lmn.suff(Y = Y, X = X, V = toy_acf(lambda), Vtype = "acf"))
})

# estimate theta = (alpha, lambda)
Tz <- SuperGauss::Toeplitz(n = n)

# sufficient statistics for toy model
toy_suff <- function(theta) {
  X <- cbind(1, xseq^theta[1])
  Tz$setAcf(acf = toy_acf(theta[2]))
  lmn.suff(Y = Y, X = X, V = Tz, Vtype = "acf")
}

# _negative_ profile likelihood for theta
toy_prof <- function(theta) {
  if(theta[2] < 0) return(-Inf) # range restriction lambda > 0
  suff <- toy_suff(theta)
  -lmn.prof(suff = suff)
}

# MLE of theta
opt <- optim(par = c(alpha, lambda), # starting value
                 fn = toy_prof) # objective function
theta_mle <- opt$par

# MLE of (Beta, Sigma)
suff <- toy_suff(theta_mle)
Beta_mle <- suff$Bhat
Sigma_mle <- suff$S/suff$n
Theta_mle <- c(theta_mle, t(Beta_mle), cov2sigrho(Sigma_mle))
names(Theta_mle) <- c("alpha", "lambda",
                      "Beta_01", "Beta_02", "Beta_11", "Beta_12",
                      "sigma_1", "sigma_2", "rho_12")
Theta_true <- c(alpha, lambda, t(Beta), cov2sigrho(Sigma))
names(Theta_true) <- names(Theta_mle)
rbind(true = Theta_true, mle = Theta_mle)

# full _negative_ loglikelihood for the toy model
toy_nll <- function(Theta) {
  theta <- Theta[1:2]
  Beta <- rbind(Theta[3:4], Theta[5:6])
  Sigma <- diag(Theta[7:8]^2)
  Sigma[1,2] <- Sigma[2,1] <- Theta[7]*Theta[8]*Theta[9]
  # calculate loglikelihood
  suff <- toy_suff(theta)
  -lmn.loglik(Beta = Beta, Sigma = Sigma, suff = suff)
}

# variance estimator
Theta_ve <- solve(numDeriv::hessian(func = toy_nll, x = Theta_mle))
Theta_se <- sqrt(diag(Theta_ve)) # standard errors

# display
disp <- rbind(true = c(alpha, lambda, t(Beta), cov2sigrho(Sigma)),
              mle = Theta_mle,
              se = Theta_se)


#--- gcir example --------------------------------------------------------------

gcir_sim <- function(N, dt, Theta, x0) {
  # parameters
  gamma <- Theta[1]
  mu <- Theta[2]
  sigma <- Theta[3]
  lambda <- Theta[4]
  Rt <- rep(NA, N+1)
  Rt[1] <- x0
  for(ii in 1:N) {
    Rt[ii+1] <- rnorm(1, mean = Rt[ii] - gamma * (Rt[ii] - mu) * dt,
                      sd = sigma * Rt[ii]^lambda * sqrt(dt))
  }
  Rt
}

gcir_suff <- function(lambda) {
  lmn.suff(Y = Y, X = X,
           V = exp(lambda * lR2) * dt, Vtype = "diag")
}

gcir_prof <- function(lambda) {
  if(lambda <= 0) return(Inf)
  -lmn.prof(suff = gcir_suff(lambda))
}

plot(Rt)


# plausible values for gcir model (Ait-Sahalia 1999)
#gamma <- .0876
#mu <- .0844
#sigma <- .7791
#beta <- 1.48

#deltat <- 1/12
#N <- (98-63+1)*12

Theta <- c(gamma = .07, mu = .01, sigma = .6, lambda = .9)
dt <- 1/12
N <- 12 * 20


set.seed(7)
# simulate data
Rt <- gcir_sim(N = N, dt = dt, Theta = Theta, x0 = Theta["mu"])
anyNA(Rt)
# precompute these
Y <- matrix(diff(Rt))
X <- cbind(-Rt[1:N], 1) * dt
# since Rt^(2*lambda) is calculated as exp(2*lambda * log(Rt)),
# precompute 2*log(Rt) to speed up calculations
lR2 <- 2 * log(Rt[1:N])
opt <- optimize(f = gcir_prof,
                interval = c(.001, 10))
lambda_mle <- opt$minimum
suff <- gcir_suff(lambda_mle)
suff$Bhat

# posterior inference
prior <- lmn.prior(p = 2, q = 1)

gcir_marg <- function(lambda) {
  suff <- gcir_suff(lambda)
  post <- lmn.post(suff, prior)
  lmn.marg(suff = suff, prior = prior, post = post)
}

lambda_mode <- optimize(f = gcir_marg,
                        interval = c(.01, 10),
                        maximum = TRUE)$maximum
lambda_quad <- -numDeriv::hessian(func = gcir_marg, x = lambda_mode)[1]
lambda_rng <- lambda_mode + c(-5,5) * 1/sqrt(lambda_quad)

npts <- 1000
lambda_seq <- seq(lambda_rng[1], lambda_rng[2], len = npts)
lambda_lpdf <- sapply(lambda_seq, gcir_marg)
lambda_pdf <- exp(lambda_lpdf - max(lambda_lpdf))
lambda_pdf <- lambda_pdf / sum(lambda_pdf) / (lambda_seq[2]-lambda_seq[1])

plot(lambda_seq, lambda_pdf, type = "l")

npost <- 5e4 # number of posterior draws

# marginal sampling from p(lambda | R)
lambda_post <- sample(lambda_seq, size = npost, prob = lambda_pdf,
                      replace = TRUE)

# conditional sampling from p(B, Sigma | lambda, R)
BSig_post <- lapply(lambda_post, function(lambda) {
  lmn.post(gcir_suff(lambda), prior)
})
BSig_post <- list2mniw(BSig_post) # convert to vectorized mniw format
BSig_post <- mniw::rmniw(npost,
                          Lambda = BSig_post$Lambda,
                          Omega = BSig_post$Omega,
                          Psi = BSig_post$Psi,
                          nu = BSig_post$nu)
# convert to Theta = (gamma, mu, sigma, lambda)
Theta_post <- cbind(gamma = BSig_post$X[1,1,],
                    mu = BSig_post$X[2,1,]/BSig_post$X[1,1,],
                    sigma = sqrt(BSig_post$V[1,1,]),
                    lambda = lambda_post)


# keep only draws for which gamma, mu > 0
ikeep <- pmin(Theta_post[,1], Theta_post[,2]) > 0
Theta_post <- Theta_post[ikeep,]
# convert mu to log scale for plotting purposes
Theta_post[,"mu"] <- log10(Theta_post[,"mu"])

# posterior distributions and true parameter values
Theta_names <- c("gamma", "log[10](mu)", "sigma", "lambda")
Theta_true <- Theta
Theta_true["mu"] <- log10(Theta_true["mu"])
par(mfrow = c(2,2), mar = c(2,2,3,.5)+.5)
for(ii in 1:ncol(Theta_post)) {
  hist(Theta_post[,ii], breaks = 40, freq = FALSE,
       xlab = "", ylab = "",
       main = parse(text = paste0("p(",
                                  Theta_names[ii], "*\" | \"*bold(R))")))
  abline(v = Theta_true[ii], col = "red", lwd = 2)
  if(ii == 1) {
    legend("topright", inset = .05,
           legend = c("Posterior Distribution", "True Parameter Value"),
           lwd = c(NA, 2), pch = c(22, NA), seg.len = 1.5,
           col = c("black", "red"), bg = c("white", NA), cex = .85)
  }
}

#--- fbm example ----------------------------------------------------------


# autocorrelation of fBM increments
fbm_acf <- function(alpha, dT, N) {
  if (N == 1) {
    acf <- dT^alpha
  }
  else {
    acf <- (dT * (0:N))^alpha
    acf <- 0.5 * (acf[1:N+1] + c(acf[2], acf[1:(N-1)]) - 2 * acf[1:N])
  }
  acf
}

# simulate fbm with drift

qq <- 2 # number of dimensions
N <- 1800 # number of observations
dT <- 1/60 # interobservation time

alpha <- .7 # subdiffusion exponent
mu <- rnorm(qq) # linear drift term
## mu <- rnorm(qq, sd = .1) # linear drift term
## eta <- rnorm(qq, sd = .1) # quadratic drift term
Sigma <- crossprod(matrix(rnorm(qq^2), qq, qq)) # scale matrix
## mu <- rep(0, qq)
## eta <- rep(0, qq)
## Sigma <- diag(qq)

# fbm increments
acf <- fbm_acf(alpha = alpha, dT = dT, N = N-1)
dZ <- mniw::rMNorm(n = 1, RowV = toeplitz(acf), ColV = Sigma)
# drift term
tseq <- (0:(N-1)) * dT
## dr <- tseq %o% mu + tseq^2 %o% eta
dr <- tseq %o% mu
# observations
Xt <- apply(rbind(0, dZ), 2, cumsum) + dr

# plot it
par(mfrow = c(1,3), mar = c(4, 4, .5, .5)+.1)
x1lab <- expression(X[1*t])
x2lab <- expression(X[2*t])
tlab <- expression(t)
plot(Xt, type = "l", xlab = x1lab, ylab = x2lab)
plot(tseq, Xt[,1], type = "l", xlab = tlab, ylab = x1lab)
plot(tseq, Xt[,2], type = "l", xlab = tlab, ylab = x2lab)

# parameter inference
require(SuperGauss)

dX <- apply(Xt, 2, diff) # observation increments
Tz <- SuperGauss::Toeplitz(n = N-1)
ddr <- apply(cbind(tseq, tseq^2), 2, diff) # drift increments

# sufficient statistics of nuisance parameters given alpha
fbm_suff <- function(alpha) {
  acf <- fbm_acf(alpha = alpha, dT = dT, N = N-1)
  Tz$setAcf(acf = acf)
  lmn.suff(Y = dX, X = ddr, V = Tz, Vtype = "acf")
}

# profile likelihood function
fbm_prof <- function(alpha) {
  lmn.prof(suff = fbm_suff(alpha))
}

# plot it
alpha_seq <- seq(.4, 1, len = 200)
ll_seq <- sapply(alpha_seq, fbm_prof)
par(mfrow = c(1,1))
plot(alpha_seq, ll_seq, type = "l",
     xlab = expression(alpha),
     ylab = expression("\u2113"[prof](alpha*" | "*bold(X)[t])))
abline(v = alpha, lty = 2) # add true value

# calculate the mle
# alpha by profile likelihood
alpha_mle <- optimize(f = fbm_prof,
                      interval = c(.01, .99), maximum = TRUE)$maximum
# remaining parameter mles
suff_mle <- fbm_suff(alpha = alpha_mle)
mu_mle <- suff_mle$Bhat[1,]
eta_mle <- suff_mle$Bhat[2,]
Sigma_mle <- suff_mle$S/suff_mle$n

# calculate standard errors using hessian and parametrization:
# theta = (alpha, mu, eta, sigma_1, sigma_2, rho)

# fbm + drift _negative_ loglikelihood for all parameters
fbm_nll <- function(theta) {
  # extract (alpha, B, Sigma) from theta
  alpha <- theta[1]
  Beta <- rbind(theta[2:3], theta[4:5])
  Sigma <- diag(theta[6:7]^2)
  Sigma[1,2] <- Sigma[2,1] <- theta[6]*theta[7]*theta[8]
  # calculate loglikelihood
  suff <- fbm_suff(alpha)
  -lmn.loglik(Beta = Beta, Sigma = Sigma, suff = suff)
}

theta_mle <- c(alpha_mle, mu_mle, eta_mle,
               sqrt(Sigma_mle[1,1]), sqrt(Sigma_mle[2,2]),
               Sigma_mle[1,2]/sqrt(Sigma_mle[1,1]*Sigma_mle[2,2]))

