# problem dimensions
p <- 2
q <- 4

# default noninformative prior pi(Beta, Sigma) ~ |Sigma|^(-(q+1)/2)
lmn_prior(p, q)

# pi(Sigma) ~ |Sigma|^(-(q+1)/2)
# Beta | Sigma ~ Matrix-Normal(0, I, Sigma)
lmn_prior(p, q, Lambda = 0, Omega = 1)

# Sigma = diag(q)
# Beta ~ Matrix-Normal(0, I, Sigma = diag(q))
lmn_prior(p, q, Lambda = 0, Omega = 1, nu = NA)
