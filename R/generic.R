#' @export
#' @noRd
coef.aftosmac <- function(object, ...){
  z <- object
  if ("1" %in% names(z$convergence) | "2" %in% names(z$convergence) |
      ! 0 %in% z$convergence) {
    stop("Failed to get a converging result.")
  }
  if (!is.aftosmac(z)) stop("Most be aftosmac class")
  ans <- z["call"]
  out <- z$coefficients
  names(out) <- z$vari.name
  return(out)
}

#' @export
#' @noRd
residuals.aftosmac <- function(object, ...){
  z <- object
  if ("1" %in% names(z$convergence) | "2" %in% names(z$convergence) |
      ! 0 %in% z$convergence) {
    stop("Failed to get a converging result.")
  }
  if (!"(Intercept)" %in% colnames(z$x)){
    z$x <- as.matrix(cbind(1, z$x))
  }
  if ("Scale" %in% names(z$coefficients)){
    z$coefficients <- z$coefficients[-which(names(z$coefficients) == "Scale")]
  }
  if (!is.aftosmac(z)) stop("Most be aftosmac class")
  ans <- z["call"]
  out <- log(z$y[, 1]) - as.matrix(z$x) %*% z$coefficients
  out
}

#' @export
#' @noRd
vcov.aftosmac <- function(object, ...){
  z <- object
  if ("1" %in% names(z$convergence) | "2" %in% names(z$convergence) |
      ! 0 %in% z$convergence) {
    stop("Failed to get a converging result.")
  }
  if (!is.aftosmac(z)) stop("Most be aftosmac class")
  ans <- z["call"]
  out <- z$covmat
  out
}


#' @export
#' @noRd
predict.aftosmac <- function(object, newdata = NULL, se.fit = FALSE, ...){
  z <- object
  if ("1" %in% names(z$convergence) | "2" %in% names(z$convergence) |
      ! 0 %in% z$convergence) {
    stop("Failed to get a converging result.")
  }
  z$x <- as.matrix(z$x)
  out <- NULL
  if (!"(Intercept)" %in% colnames(z$x)){
    z$x <- as.matrix(cbind(1, z$x))
  }
  if ("Scale" %in% names(z$coefficients)){
    z$coefficients <- z$coefficients[-which(names(z$coefficients) == "Scale")]
  }
  if (!is.aftosmac(z)) stop("Most be aftosmac class")
  ans <- z["call"]
  if (is.null(newdata)) {
    out$fit <- exp(drop(z$x %*% z$coefficients)) # originally on log scale
  }
  if (!is.null(newdata)) {
    newdata0 <- as.matrix(newdata, ncol = length(z$coefficients))
    out$fit <- exp(drop(newdata0 %*% z$coefficients))
  }
  out$se.fit <- NA
  out
}
