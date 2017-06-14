#--- basic testing --------------------------------------------------------------

require(LMN)

# DurbinLevinson

calcMode <- 0
d <- sample(5, 1)
n <- sample(100, 1)
k <- ifelse(calcMode == 2, d, sample(10,1))
acf <- exp(-2*(1:n)/n)
T <- toeplitz(acf)
X <- matrix(rnorm(d*n),n,d)
Y <- matrix(rnorm(k*n),n,k)
if(calcMode == 1) {
  IP <- crossprod(X, solve(T, X))
} else {
  IP <- crossprod(X, solve(T, Y))
  if(calcMode == 2) IP <- diag(IP)
}
ldV <- 2*sum(log(diag(chol(T))))
if(calcMode == 1) {
  DLE <- LMN:::DurbinLevinsonEigen(X = t(X), Y = matrix(0), acf = acf,
                                   calcMode = calcMode)
} else {
  DLE <- LMN:::DurbinLevinsonEigen(X = t(X), Y = t(Y), acf = acf,
                                   calcMode = calcMode)
}
if(calcMode == 1) {
  DLB <- LMN:::DurbinLevinsonBase(X = X, Y = X, acf = acf,
                            calcMode = calcMode)
} else {
  DLB <- LMN:::DurbinLevinsonBase(X = X, Y = Y, acf = acf,
                            calcMode = calcMode)
}
cbind(Eigen = c(IP = max(abs(IP - DLE$IP)),
        ldV = abs(ldV - DLE$ldV)),
      Base = c(IP = max(abs(IP - DLB$IP)),
        ldV = abs(ldV - DLB$ldV)))

#--- speed comparisons ----------------------------------------------------------

case.par <- expand.grid(d = c(1, 2, 5, 10, 50), k = c(1, 2, 5, 10, 50),
                        n = c(100, 1e3, 2e3), calcMode = 0:2)
# for calcMode = 2 need d = k
case.par <- case.par[with(case.par, {(calcMode != 2) | (d == k)}),]
ncases <- nrow(case.par)

nreps <- 10
Time.Eigen <- rep(NA, ncases)
Time.Base <- rep(NA, ncases)
for(ii in 1:ncase) {
  message("ii = ", ii)
  cp <- case.par[ii,]
  n <- cp$n
  d <- cp$d
  k <- cp$k
  calcMode <- cp$calcMode
  acf <- exp(-2*(1:n)/n - .01*rnorm(n))
  X <- matrix(rnorm(n*d),n,d)
  Y <- matrix(rnorm(n*k),n,k)
  # Eigen
  tm <- system.time({
    for(jj in 1:nreps) {
      if(calcMode == 1) {
        DLE <- LMN:::DurbinLevinsonEigen(X = t(X), Y = matrix(0),
                                         acf = acf,
                                         calcMode = calcMode)
      } else {
        DLE <- LMN:::DurbinLevinsonEigen(X = t(X), Y = t(Y), acf = acf,
                                         calcMode = calcMode)
      }
    }
  })
  Time.Eigen[ii] <- tm[3]
  # Base
  tm <- system.time({
    for(j in 1:nreps) {
      if(calcMode == 1) {
        DLB <- LMN:::DurbinLevinsonBase(X = X, Y = X, acf = acf,
                                        calcMode = calcMode)
      } else {
        DLB <- LMN:::DurbinLevinsonBase(X = X, Y = Y, acf = acf,
                                        calcMode = calcMode)
      }
    }
  })
  Time.Base[ii] <- tm[3]
}

#--- one case at a time ---------------------------------------------------------

dl.calc <- function(n, k, d, calcMode, nreps = 1) {
  if((calcMode == 2) && (k != d)) stop("d == k for calcMode == 2.")
  acf <- exp(-2*(1:n)/n - .01*rnorm(n))
  X <- matrix(rnorm(n*d),n,d)
  Y <- matrix(rnorm(n*k),n,k)
  # Eigen
  tme <- system.time({
    for(jj in 1:nreps) {
      if(calcMode == 1) {
        DLE <- LMN:::DurbinLevinsonEigen(X = t(X), Y = matrix(0),
                                         acf = acf,
                                         calcMode = calcMode)
      } else {
        DLE <- LMN:::DurbinLevinsonEigen(X = t(X), Y = t(Y), acf = acf,
                                         calcMode = calcMode)
      }
    }
  })
  # Base
  tmb <- system.time({
    for(j in 1:nreps) {
      if(calcMode == 1) {
        DLB <- LMN:::DurbinLevinsonBase(X = X, Y = X, acf = acf,
                                        calcMode = calcMode)
      } else {
        DLB <- LMN:::DurbinLevinsonBase(X = X, Y = Y, acf = acf,
                                        calcMode = calcMode)
      }
    }
  })
  c(Eigen = tme[3], Base = tmb[3])
}

n <- 2e3
k <- 1
d <- 1
calcMode <- 1
nreps <- 100
dl.calc(n, k, d, calcMode, nreps)

#--- toeplitz2 ------------------------------------------------------------------

R <- .1 * 1:3
C <- -(2:5) * .27
R[1] <- C[1]
LMN:::toeplitz2(col = C, row = R)

#--- check that suff recovers the MLE -------------------------------------------

require(LMN)
require(numDeriv)

calc.grad <- TRUE
case.par <- expand.grid(p = c(-1, 0, 1, 3, 5), q = c(1, 3, 5),
                        type = c("scalar", "diag", "acf", "V"),
                        noSigma = c(TRUE, FALSE))
# remove noBeta && noSigma
case.par <- case.par[with(case.par, {p != 0 | !noSigma}),]
ncases <- nrow(case.par)
n <- 20
if(calc.grad) {
  MaxGrad <- rep(NA, ncases)
}

for(ii in 1:ncases) {
  # get parameters
  cp <- case.par[ii,]
  p <- cp$p
  q <- cp$q
  noBeta <- p == 0
  type <- as.character(cp$type)
  noSigma <- cp$noSigma
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
  # calculate with lmn.suff
  if(type == "acf") {
    suff <- lmn.suff(Y = Y, X = XX, acf = acf)
  } else {
    suff <- lmn.suff(Y = Y, X = XX, V = VV)
  }
  Beta.hat <- suff$Beta.hat
  Sigma.hat <- suff$S/suff$n
  # check that gradient equals 0
  if(!noBeta && !noSigma) {
    logdens <- function(theta) {
      lMnorm(X = Y, RowV = VR,
             Mu = XR %*% matrix(theta[1:(p*q)],p,q),
             ColV = ltri2Sig(theta[p*q + (1:(q*(q+1)/2))]))
    }
    ld.grad <- grad(logdens, x = c(Beta.hat, Sig2ltri(Sigma.hat)))
  } else if(noBeta && !noSigma) {
    logdens <- function(theta) {
      lMnorm(X = Y, RowV = VR,
             Mu = matrix(0,n,q), ColV = ltri2Sig(theta))
    }
    ld.grad <- grad(logdens, x = Sig2ltri(Sigma.hat))
  } else if(!noBeta && noSigma) {
    logdens <- function(theta) {
      lMnorm(X = Y, RowV = VR,
             Mu = XR %*% matrix(theta[1:(p*q)],p,q),
             ColV = diag(q))
    }
    ld.grad <- grad(logdens, x = Beta.hat)
  } else {
    ld.grad <- NA
  }
  if(calc.grad) {
    MaxGrad[ii] <- max(abs(ld.grad))
  } else {
    expect_equal(max(abs(ld.grad)),0)
  }
}

#--- test profile likelihood ----------------------------------------------------

calc.diff <- TRUE
case.par <- expand.grid(p = c(-1, 0, 1, 3, 5), q = c(1, 3, 5),
                        type = c("scalar", "diag", "acf", "V"),
                        noSigma = c(TRUE, FALSE))
# remove noBeta && noSigma
case.par <- case.par[with(case.par, {p != 0 | !noSigma}),]
ncases <- nrow(case.par)
n <- 20
if(calc.grad) {
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
  # calculate with lmn.suff
  if(type == "acf") {
    suff <- lmn.suff(Y = Y, X = XX, acf = acf)
  } else {
    suff <- lmn.suff(Y = Y, X = XX, V = VV)
  }
  Beta.hat <- suff$Beta.hat
  Sigma.hat <- suff$S/suff$n
  # check that profile likelihood equals long hand calculation
  if(!noBeta) {
    Mu <- XR %*% Beta.hat
  } else {
    Mu <- matrix(0,n,q)
  }
  if(!noSigma) {
    Sigma <- Sigma.hat
  } else {
    Sigma <- diag(q)
  }
  llR <- lMnorm(X = Y, Mu = Mu, RowV = VR, ColV = Sigma)
  llp <- lmn.prof(suff = suff, noSigma = noSigma)
  if(calc.diff) {
    MaxDiff[ii] <- abs(llR - llp)
  } else {
    expect_equal(llR, llp)
  }
}

#--- marginal and conditional ---------------------------------------------------

# convert a matrix to matlab format by printing assignment commands
matrixR2M <- function(X, nm, debug = FALSE) {
  if(debug) browser()
  vp <- function(x) {
    cat("[")
    cat(x, sep = ", ")
    cat("]")
  }
  if(!missing(nm)) {
    cat(paste0(nm, " = "))
  }
  if(is.vector(X) || nrow(X) == 1) {
    if(length(X) == 1) {
      cat(X)
    } else {
      vp(X)
    }
    cat(";\n")
  } else {
    if(ncol(X) == 1) {
      vp(X)
      cat("';\n")
    } else {
      cat("[")
      for(ii in 1:nrow(X)) {
        if(ii > 1) {
          cat("  ")
        }
        vp(X[ii,])
        if(ii < nrow(X)) {
          cat("; ...\n")
        }
      }
      cat("];\n")
    }
  }
}

matrixR2M(Lambda, debug = FALSE)

# problem specification
cp <- list(p = 2, q = 3, type = "full",
           noSigma = FALSE, noOmega = TRUE, noPsi = TRUE,
           noPrior = TRUE, preSuff = TRUE)
n <- 10
p <- cp$p
q <- cp$q
noBeta <- p == 0
type <- as.character(cp$type)
noSigma <- cp$noSigma
noOmega <- cp$noOmega
noPsi <- cp$noPsi
noPrior <- cp$noPrior
#preSuff <- cp$preSuff
# if(noBeta && noSigma) noSigma <- FALSE
if(noBeta) noOmega <- TRUE

# set up data
Y <- rMnorm(n,q)
if(p == 0) {
  X <- 0
} else if(p == -1) {
  p <- 1
  X <- rnorm(1)
  XR <- matrix(X, n, 1)
} else {
  X <- rMnorm(n,p)
  XR <- X
}
acf <- exp(-2*(1:n)/n)
diagV <- rexp(n)
if(type == "scalar") {
  V <- acf[1]
  VR <- diag(rep(acf[1], n))
} else if(type == "diag") {
  V <- diagV
  VR <- diag(diagV)
} else {
  V <- toeplitz(acf)
  VR <- V
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
  nu <- rexp(1,rate=1/q) + q
} else {
  nu <- NA
}

# set up parameters
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

# ok print all
Vars <- c("Y", "X", "V",
          "Beta", "Sigma",
          "Lambda", "Omega", "Psi", "nu")
for(ii in Vars) {
  eval(parse(text = paste0(ii, " = round(", ii, ",2)")))
  eval(parse(text = paste0("matrixR2M(X = ", ii,
               ", nm = \"", ii, "\")")))
}
XR <- round(XR,2)
VR <- round(VR,2)


# prior
prior <- list(Lambda = Lambda, Omega = Omega, Psi = Psi, nu = nu)
lpi <- lMNIW(X = Beta, V = Sigma,
             Lambda = Lambda, Omega = Omega, Psi = Psi, nu = nu)
# loglik
if(noBeta) {
  Mu <- matrix(0,n,q)
} else {
  Mu <- XR %*% Beta
}
ll <- lMnorm(Y, Mu, VR, Sigma)
# logpost
if(type == "acf") {
  suff <- lmn.suff(Y = Y, X = X, acf = acf)
} else {
  suff <- lmn.suff(Y = Y, X = X, V = V)
}
post <- lmn.post(suff, prior = prior, debug = FALSE)
lpc <- lMNIW(Beta, Sigma,
             post$Lambda, post$Omega, post$Psi, post$nu)
# logmarg
lpm <- lmn.marg(suff, post = post, debug = FALSE)

(ll+lpi) - (lpc+lpm)
