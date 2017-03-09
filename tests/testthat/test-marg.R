library(LMN)
source("lmn-testfunctions.R")
context("Marginal")

test_that("Lik * Prior = Marg * Cond.", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(-1, 0, 1, 3, 4), q = c(1, 3, 5),
                          type = c("scalar", "diag", "acf", "V"),
                          noSigma = c(TRUE, FALSE),
                          noOmega = c(TRUE, FALSE),
                          noPsi = c(TRUE, FALSE),
                          preSuff = c(TRUE, FALSE),
                          noPrior = c(TRUE, FALSE))
  # remove noBeta && (noSigma || noOmega)
  case.par <- case.par[with(case.par, {p != 0 | (!noSigma & !noOmega)}),]
  # when noPrior always have noOmega = noPsi = TRUE
  case.par <- case.par[with(case.par, {!noPrior | (noOmega & noPsi)}),]
  ncases <- nrow(case.par)
  rownames(case.par) <- 1:ncases
  n <- 20
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    # get parameters
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    noBeta <- p == 0
    type <- as.character(cp$type)
    noSigma <- cp$noSigma
    noOmega <- cp$noOmega
    noPsi <- cp$noPsi
    preSuff <- cp$preSuff
    noPrior <- cp$noPrior
    if(noBeta) noOmega <- TRUE
    # set up data
    Y <- rMnorm(n,q)
    if(p == 0) {
      XX <- 0
    } else if(p == -1) {
      p <- 1
      XX <- rnorm(1)
      XR <- matrix(XX, n, 1)
    } else {
      XX <- rMnorm(n,p)
      XR <- XX
    }
    acf <- exp(-2*(1:n)/n)
    diagV <- rexp(n)
    if(type == "scalar") {
      VV <- acf[1]
      VR <- diag(rep(acf[1], n))
    } else if(type == "diag") {
      VV <- diagV
      VR <- diag(diagV)
    } else {
      VV <- toeplitz(acf)
      VR <- VV
    }
    # set up hyperparameters
    if(!noBeta) {
      Lambda <- rMnorm(p,q)
    } else {
      Lambda <- NULL
    }
    if(!noOmega) {
      Omega <- crossprod(rMnorm(p))
    } else {
      Omega <- ifelse(noBeta, NA, 0)
    }
    if(!noPsi) {
      Psi <- crossprod(rMnorm(q))
    } else {
      Psi <- 0
    }
    if(!noSigma) {
      nu <- ifelse(noPrior, 0, rexp(1,rate=1/q) + q)
    } else {
      nu <- NA
    }
    # setup parameters
    if(!noBeta) {
      Beta <- rMnorm(p,q)
    } else {
      Beta <- NULL
    }
    if(!noSigma) {
      Sigma <- crossprod(rMnorm(q))
    } else {
      Sigma <- diag(q)
    }
    # prior
    prior <- list(Lambda = Lambda, Omega = Omega, Psi = Psi, nu = nu)
    #if(noPrior) {
    #  lpi <- 0
    #} else {
      lpi <- lMNIW(X = Beta, V = Sigma,
                   prior$Lambda, prior$Omega, prior$Psi, prior$nu)
    #}
    # loglik
    if(noBeta) {
      Mu <- matrix(0,n,q)
    } else {
      Mu <- XR %*% Beta
    }
    ll <- lMnorm(Y, Mu, VR, Sigma)
    # logpost & logmarg
    if(type == "acf") {
      Data <- list(Y = Y, X = XX, acf = acf)
    } else {
      Data <- list(Y = Y, X = XX, V = VV)
    }
    suff <- do.call(lmn.suff, args = Data)
    if(noPrior) {
      if(preSuff) {
        post <- lmn.post(suff, noSigma = noSigma)
        # get prior from post
        lpm <- lmn.marg(suff, post = post)
      } else {
        post <- do.call(lmn.post,
                        args = c(Data, list(noSigma = noSigma)))
        # explicitly pass prior
        lpm <- do.call(lmn.marg,
                       args = c(Data, list(noSigma = noSigma,
                         prior = NULL)))
      }
    } else {
      # infer noSigma from prior
      if(preSuff) {
        post <- lmn.post(suff, prior = prior)
        # get prior from post
        lpm <- lmn.marg(suff, post = post)
      } else {
        post <- do.call(lmn.post, args = c(Data, list(prior = prior)))
        # explicitly pass prior
        lpm <- do.call(lmn.marg,
                       args = c(Data, list(post = post, prior = prior)))
      }
    }
    lpc <- lMNIW(Beta, Sigma,
                 post$Lambda, post$Omega, post$Psi, post$nu)
    # smallest of relative and absolute error
    mxd <- abs((lpi+ll) - (lpc+lpm))
    mxd <- min(mxd, mxd/abs(lpi+ll))
    if(calc.diff) {
      MaxDiff[ii] <- mxd
    } else {
      expect_equal(mxd, 0, info = paste0("Case ", ii))
    }
  }
})
