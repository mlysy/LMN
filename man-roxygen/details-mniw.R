#' @details The Matrix-Normal Inverse-Wishart (MNIW) distribution \eqn{(\boldsymbol{B}, \boldsymbol{\Sigma}) \sim \textrm{MNIW}(\boldsymbol{\Lambda}, \boldsymbol{\Omega}, \boldsymbol{\Psi}, \nu)}{(B, \Sigma) ~ MNIW(\Lambda, \Omega, \Psi, \nu)} on random matrices \eqn{\boldsymbol{X}_{p \times q}}{X_(p x q)} and symmetric positive-definite \eqn{\boldsymbol{\Sigma}_{q \times q}}{\Sigma_(q x q)} is defined as
#' \deqn{
#' \begin{array}{rcl}
#' \boldsymbol{\Sigma} & \sim & \textrm{Inverse-Wishart}(\boldsymbol{\Psi}, \nu) \\
#' \boldsymbol{B} \mid \boldsymbol{\Sigma} & \sim & \textrm{Matrix-Normal}(\boldsymbol{\Lambda}, \boldsymbol{\Omega}^{-1}, \boldsymbol{\Sigma}),
#' \end{array}
#' }{
#' \Sigma ~ Inverse-Wishart(\Psi, \nu)
#' }
#' \deqn{
#' \vspace{-1em}
#' }{
#' B | \Sigma ~ Matrix-Normal(\Lambda, \Omega^{-1}, \Sigma),
#' }
#' where the Matrix-Normal distribution is defined in \code{\link{lmn.suff}}.
#'
