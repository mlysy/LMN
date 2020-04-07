#' @param prior A list specifying the prior MNIW distribution, with elements:
#' \describe{
#'   \item{`Lambda`}{A `p x q` mean matrix for `Beta`.}
#'   \item{`Omega`}{A `q x q` column-wise precision matrix for `Beta`.}
#'   \item{`Psi`}{A `p x p` scale matrix for `Sigma`.}
#'   \item{`nu`}{A scalar degrees-of-freedom parameter for `Sigma`.}
#' }
#' Any omitted or null elements default to zero.  `prior = NULL` sets all MNIW parameters to zero.
