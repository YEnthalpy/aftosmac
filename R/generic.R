#' @export
#' @noRd
coef.aftosmac <- function(object, ...){
  z <- object
  if (!is.aftosmac(z)) stop("Most be aftosmac class")
  ans <- z["call"]
  out <- z$coefficients
  names(out) <- z$vari.name
  out
}

#' @export
#' @noRd
residuals.aftosmac <- function(object, ...){
  z <- object
  if (!"(Intercept)" %in% colnames(z$x)){
    z$x <- as.matrix(cbind(1, z$x))
  }
  if ("Scale" %in% names(z$coefficients)){
    z$coefficients <- z$coefficients[-which(names(z$coefficients) == "Scale")]
  }
  if (!is.aftosmac(z)) stop("Most be aftosmac class")
  ans <- z["call"]
  if (!"(Intercept)" %in% colnames(z$x)){
    z$x <- as.matrix(cbind(1, z$x))
  }
  out <- log(z$y[, 1]) - as.matrix(z$x) %*% z$coefficients
  out
}

#' @export
#' @noRd
vcov.aftosmac <- function(object, ...){
  z <- object
  if (!is.aftosmac(z)) stop("Most be aftosmac class")
  ans <- z["call"]
  out <- z$covmat
  out
}


#' @export
#' @noRd
predict.aftosmac <- function(object, newdata = NULL, se.fit = FALSE, ...){
  z <- object
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
