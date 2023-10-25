#' aftosmac: Optimal Subsampling Methods for Accelerated Failure Time Models
#'
#' A package that use optimal subsampling method to fast estimate the Accelerated
#' Failure Time (AFT) models. This package implements developed inference 
#' procedures for parametric AFT models and semi-parametric AFT models
#' with both the rank-based approach and the least squares approach by the optimal 
#' subsampling approach. For the parametric AFT models, the package only 
#' supports the Weibull AFT model and estimates both the scale parameter and 
#' regression coefficients. For the rank-based approach, the package currently 
#' only support the Gehan's weight with the smoothed estimating function 
#' derived by the induced smoothed procedure. We developed two types of optimal
#' subsampling probabilities (SSPs) based on the A-optimal criteria and the
#' L-optimal criteria from the optimal design of experiment.

#'
#' @aliases aftosmac-packages
#' @references Yang, Z., Wang, H., Yan, J. (2023) Subsampling Approach for Least Squares
#' Fitting of Semi-parametric Accelerated Failure Time Models to Massive
#' Survival Data.
#' \emph{Research Square preprint}, https://www.researchsquare.com/article/rs-2866415/v1.
#' @references Yang Z., Wang, H., Yan, J. (2022) Optimal Subsampling for
#' Parametric Accelerated Failure Time Models with Massive
#' Survival Data.  2022; 27): 5421–5431.
#' \emph{Statistics in Medicine}, \bold{41}(27): 5421–5431.

#'
#' @importFrom survival Surv
#' @importFrom methods getClass
#' @importFrom stats approx lsfit model.matrix rbinom rexp var
#' @import Rcpp nleqslv
#'
#' @useDynLib aftosmac, .registration = TRUE
#' @docType package
"_PACKAGE"
NULL
