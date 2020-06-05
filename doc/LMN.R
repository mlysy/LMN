## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(out.width = "\\textwidth")
suppressMessages({
  require(LMN)
  require(kableExtra)
})
cran_link <- function(pkg) paste0("[**", pkg, "**](https://CRAN.R-project.org/package=", pkg, ")")
build_cache <- FALSE

## ----toy_sim_calc, include = FALSE--------------------------------------------
require(LMN)

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
Y <- mniw::rMNorm(n = 1, Lambda = Mu,
                  SigmaR = V, SigmaC = Sigma) # response matrix
## # response matrix
## Y <- matrix(rnorm(n*qq), n, qq)
## Y <- crossprod(chol(V), Y) %*% chol(Sigma) + Mu

# plot data
par(mfrow = c(1,1), mar = c(4,4,.5,.5)+.1)
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

## ----toy_sim, ref.label = "toy_sim_calc", fig.height = 4, fig.width = 6.5, fig.cap = paste0("Simulated data from the toy model \\\\eqref{eq:toy} with $n = ", n, "$ and $q = ", qq, "$.")----
require(LMN)

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
Y <- mniw::rMNorm(n = 1, Lambda = Mu,
                  SigmaR = V, SigmaC = Sigma) # response matrix
## # response matrix
## Y <- matrix(rnorm(n*qq), n, qq)
## Y <- crossprod(chol(V), Y) %*% chol(Sigma) + Mu

# plot data
par(mfrow = c(1,1), mar = c(4,4,.5,.5)+.1)
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

## -----------------------------------------------------------------------------
suff <- lmn_suff(Y = Y, X = X, V = V, Vtype = "full")
sapply(suff, function(x) paste0(dim(as.array(x)), collapse = "x"))

## ---- results = "hold"--------------------------------------------------------
# check than dense matrix and Toeplitz matrix calculations are the same

# autocorrelation function, or first row of V_theta
toy_acf <- function(lambda) exp(-((xseq-xseq[1])/lambda)^2)

# check that calculation of suff is the same
all.equal(suff,
          lmn_suff(Y = Y, X = X, V = toy_acf(lambda), Vtype = "acf"))

# check that it's much faster to use Vtype = "acf"
system.time({
  # using dense variance matrix
  replicate(100, lmn_suff(Y = Y, X = X, V = V, Vtype = "full"))
})
system.time({
  # using Toeplitz variance matrix
  replicate(100, lmn_suff(Y = Y, X = X, V = toy_acf(lambda), Vtype = "acf"))
})

## -----------------------------------------------------------------------------
# pre-allocate memory for Toeplitz matrix calcuations
Tz <- SuperGauss::Toeplitz$new(N = n)

# sufficient statistics for the toy model
# sufficient statistics for toy model
toy_suff <- function(theta) {
  X <- cbind(1, xseq^theta[1])
  Tz$set_acf(acf = toy_acf(theta[2]))
  lmn_suff(Y = Y, X = X, V = Tz, Vtype = "acf")
}

# _negative_ profile likelihood for theta
toy_prof <- function(theta) {
  if(theta[2] < 0) return(-Inf) # range restriction lambda > 0
  suff <- toy_suff(theta)
  -lmn_prof(suff = suff)
}

# MLE of theta
opt <- optim(par = c(alpha, lambda), # starting value
             fn = toy_prof) # objective function
theta_mle <- opt$par

# MLE of (Beta, Sigma)
suff <- toy_suff(theta_mle)
Beta_mle <- suff$Bhat
Sigma_mle <- suff$S/suff$n

# display:
# convert variance matrix to vector of standard deviations and correlations.
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

Theta_mle <- c(theta_mle, t(Beta_mle), cov2sigrho(Sigma_mle))
names(Theta_mle) <- c("alpha", "lambda",
                      "beta_01", "beta_02", "beta_11", "beta_12",
                      "sigma_1", "sigma_2", "rho_12")
signif(Theta_mle, 2)

## -----------------------------------------------------------------------------
# full _negative_ loglikelihood for the toy model
toy_nll <- function(Theta) {
  theta <- Theta[1:2]
  Beta <- rbind(Theta[3:4], Theta[5:6])
  Sigma <- diag(Theta[7:8]^2)
  Sigma[1,2] <- Sigma[2,1] <- Theta[7]*Theta[8]*Theta[9]
  # calculate loglikelihood
  suff <- toy_suff(theta)
  -lmn_loglik(Beta = Beta, Sigma = Sigma, suff = suff)
}

# uncertainty estimate:
# variance estimator
Theta_ve <- solve(numDeriv::hessian(func = toy_nll, x = Theta_mle))
Theta_se <- sqrt(diag(Theta_ve)) # standard errors

# display
tab <- rbind(true = c(alpha, lambda, t(Beta), cov2sigrho(Sigma)),
             mle = Theta_mle, se = Theta_se)
colnames(tab) <- paste0("$\\", gsub("([0-9]+)", "{\\1}", names(Theta_mle)), "$")
rownames(tab) <- c("True Value", "MLE", "Std. Error")
kableExtra::kable(as.data.frame(signif(tab,2)))

## ----gcir_setup, include = FALSE----------------------------------------------
set.seed(7) # for reproducible results

# simulate data from the gcir model
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

# true parameter values
Theta <- c(gamma = .07, mu = .01, sigma = .6, lambda = .9)
dt <- 1/12 # interobservation time (in years)
N <- 12 * 20 # number of observations (20 years)

## ----gcir_sim, include = FALSE------------------------------------------------

Rt <- gcir_sim(N = N, dt = dt, Theta = Theta, x0 = Theta["mu"])

par(mar = c(4,4,.5,.5))
plot(x = 0:N*dt, y = 100*Rt, pch = 16, cex = .8,
     xlab = "Time (years)", ylab = "Interest Rate (%)")

## ---- gcir_data, ref.label = c("gcir_setup", "gcir_sim"), echo = -1, fig.height = 4, fig.width = 6.5, fig.cap = paste0("Simulated interest rates from discrete-time approximation \\\\eqref{eq:euler} to the gCIR model \\\\eqref{eq:chan}, with $N = ", N, "$, $\\dt = 1/12$, and $\\TTh = (", paste0(Theta, collapse = ", "), ")$.")----
set.seed(7) # for reproducible results

# simulate data from the gcir model
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

# true parameter values
Theta <- c(gamma = .07, mu = .01, sigma = .6, lambda = .9)
dt <- 1/12 # interobservation time (in years)
N <- 12 * 20 # number of observations (20 years)

Rt <- gcir_sim(N = N, dt = dt, Theta = Theta, x0 = Theta["mu"])

par(mar = c(4,4,.5,.5))
plot(x = 0:N*dt, y = 100*Rt, pch = 16, cex = .8,
     xlab = "Time (years)", ylab = "Interest Rate (%)")

## ----gcir_mle-----------------------------------------------------------------
# precomputed values
Y <- matrix(diff(Rt))
X <- cbind(-Rt[1:N], 1) * dt
# since Rt^(2*lambda) is calculated as exp(2*lambda * log(Rt)),
# precompute 2*log(Rt) to speed up calculations
lR2 <- 2 * log(Rt[1:N])

# sufficient statistics for gCIR model
gcir_suff <- function(lambda) {
  lmn_suff(Y = Y, X = X,
           V = exp(lambda * lR2) * dt, Vtype = "diag")
}

# _negative_ profile likelihood for gCIR model
gcir_prof <- function(lambda) {
  if(lambda <= 0) return(Inf)
  -lmn_prof(suff = gcir_suff(lambda))
}

# MLE of Theta via profile likelihood

# profile likelihood for lambda
opt <- optimize(f = gcir_prof, interval = c(.001, 10))
lambda_mle <- opt$minimum
# conditional MLE for remaining parameters
suff <- gcir_suff(lambda_mle)
Theta_mle <- c(gamma = suff$Bhat[1,1],
               mu = suff$Bhat[2,1]/suff$Bhat[1,1],
               sigma = sqrt(suff$S[1,1]/suff$n),
               lambda = lambda_mle)
Theta_mle

## ----gcir_prior---------------------------------------------------------------
prior <- lmn_prior(p = 2, q = 1) # default prior
prior

## ----gcir_post_lambda, fig.height = 4, fig.width = 6.5------------------------

# log of marginal posterior p(lambda | R)
gcir_marg <- function(lambda) {
  suff <- gcir_suff(lambda)
  post <- lmn_post(suff, prior)
  lmn_marg(suff = suff, prior = prior, post = post)
}

# grid sampler for lambda ~ p(lambda | R)

# estimate the effective support of lambda by taking
# mode +/- 5 * sqrt(quadrature)
lambda_mode <- optimize(f = gcir_marg,
                        interval = c(.01, 10),
                        maximum = TRUE)$maximum
lambda_quad <- -numDeriv::hessian(func = gcir_marg, x = lambda_mode)[1]
lambda_rng <- lambda_mode + c(-5,5) * 1/sqrt(lambda_quad)

# plot posterior on this range
lambda_seq <- seq(lambda_rng[1], lambda_rng[2], len = 1000)
lambda_lpdf <- sapply(lambda_seq, gcir_marg) # log-pdf
# normalized pdf
lambda_pdf <- exp(lambda_lpdf - max(lambda_lpdf))
lambda_pdf <- lambda_pdf / sum(lambda_pdf) / (lambda_seq[2]-lambda_seq[1])
par(mar = c(2,4,2,.5))
plot(lambda_seq, lambda_pdf, type = "l",
     xlab = expression(lambda), ylab = "",
     main = expression(p(lambda*" | "*bold(R))))

## ----gcir_post, cache = build_cache-------------------------------------------
npost <- 5e4 # number of posterior draws

# marginal sampling from p(lambda | R)
lambda_post <- sample(lambda_seq, size = npost, prob = lambda_pdf,
                      replace = TRUE)

# conditional sampling from p(B, Sigma | lambda, R)
BSig_post <- lapply(lambda_post, function(lambda) {
  lmn_post(gcir_suff(lambda), prior)
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
apply(Theta_post, 2, min)

## ----gcir_hist, fig.width = 6.5, fig.height = 5, fig.cap = "Posterior distribution and true value for each gCIR parameter $\\TTh = (\\gamma, \\mu, \\sigma, \\lambda)$."----
# keep only draws for which gamma, mu > 0
ikeep <- pmin(Theta_post[,1], Theta_post[,2]) > 0
mean(ikeep) # a good number of draws get discarded
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

