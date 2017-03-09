#' Random matrix of iid normals
#' @param n Number or rows
#' @param p Number of columns.  When omitted defaults to \code{p = n}.
#' @return A \code{n x p} matrix of of iid elements each \code{N(0,1)}.
rMnorm <- function(n, p) {
  if(missing(p)) p <- n
  matrix(rnorm(n*p), n, p)
}
