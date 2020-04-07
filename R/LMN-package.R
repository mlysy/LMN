#' Inference for Linear Models with Nuisance Parameters.
#'
#' Efficient profile likelihood and marginal posteriors when nuisance parameters are those of linear regression models.
#'
#' @details Consider a model \eqn{p(\boldsymbol{Y} \mid \boldsymbol{B}, \boldsymbol{\Sigma}, \boldsymbol{\theta})}{p(Y | B, \Sigma, \theta)} of the form
#' \deqn{
#' \boldsymbol{Y} \sim \textrm{Matrix-Normal}(\boldsymbol{X}(\boldsymbol{\theta})\boldsymbol{B}, \boldsymbol{V}(\boldsymbol{\theta}), \boldsymbol{\Sigma}),
#' }{
#' Y ~ Matrix-Normal(X(\theta) B, V(\theta), \Sigma),
#' }
#' where \eqn{\boldsymbol{Y}_{n \times q}}{Y_(n x q)} is the response matrix, \eqn{\boldsymbol{X}(\theta)_{n \times p}}{X(\theta)_(n x p)} is a covariate matrix which depends on \eqn{\boldsymbol{\theta}}{\theta}, \eqn{\boldsymbol{B}_{p \times q}}{B_(p x q)} is the coefficient matrix, \eqn{\boldsymbol{V}(\boldsymbol{\theta})_{n \times n}}{V(\theta)_(n x n)} and \eqn{\boldsymbol{\Sigma}_{q \times q}}{\Sigma_(q x q)} are the between-row and between-column variance matrices, and (suppressing the dependence on \eqn{\boldsymbol{\theta}}{\theta}) the Matrix-Normal distribution is defined by the multivariate normal distribution
#' \eqn{
#' \textrm{vec}(\boldsymbol{Y}) \sim \mathcal{N}(\textrm{vec}(\boldsymbol{X}\boldsymbol{B}), \boldsymbol{\Sigma} \otimes \boldsymbol{V}),
#' }{
#' vec(Y) ~ N( vec(X B), \Sigma \%x\% V ),
#' }
#' where \eqn{\textrm{vec}(\boldsymbol{Y})}{vec(Y)} is a vector of length \eqn{nq} stacking the columns of of \eqn{\boldsymbol{Y}}{Y}, and \eqn{\boldsymbol{\Sigma} \otimes \boldsymbol{V}}{\Sigma \%x\% V} is the Kronecker product.
#'
#' The model above is referred to as a Linear Model with Nuisance parameters (LMN) \eqn{(\boldsymbol{B}, \boldsymbol{\Sigma})}{(B,\Sigma)}, with parameters of interest \eqn{\boldsymbol{\theta}}{\theta}.  That is, the \pkg{LMN} package provides tools to efficiently conduct inference on \eqn{\boldsymbol{\theta}}{\theta} first, and subsequently on \eqn{(\boldsymbol{B}, \boldsymbol{\Sigma})}{(B,\Sigma)}, by Frequentist profile likelihood or Bayesian marginal inference with a Matrix-Normal Inverse-Wishart (MNIW) conjugate prior on \eqn{(\boldsymbol{B}, \boldsymbol{\Sigma})}{(B,\Sigma)}.
#' @importFrom SuperGauss Toeplitz
#' @importFrom Rcpp evalCpp
#' @useDynLib LMN, .registration = TRUE
"_PACKAGE"
