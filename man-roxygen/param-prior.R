#' @param prior A list specifying the prior MNIW distribution, with elements:
#' \describe{
#'   \item{\code{Lambda}}{A \code{p x q} mean matrix for \code{Beta}.}
#'   \item{\code{Omega}}{A \code{q x q} column-wise precision matrix for \code{Beta}.}
#'   \item{\code{Psi}}{A \code{p x p} scale matrix for \code{Sigma}.}
#'   \item{\code{nu}}{A scalar degrees-of-freedom parameter for \code{Sigma}.}
#' }
#' Any omitted or null elements default to zero.  \code{prior = NULL} sets all MNIW parameters to zero.
