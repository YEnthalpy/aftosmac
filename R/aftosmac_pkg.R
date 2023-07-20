#' aftosmac: Optimal Subsampling Methods for Accelerated Failure Time Models
#'
#' A package that use optimal subsampling method to fast estimate the
#' regression coefficients and their statistical inferences for Accelerated
#' Failure Time (AFT) models.
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
