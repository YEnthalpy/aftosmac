is.aftosmac <- function(x) {
  inherits(x, "aftosmac")
}

#' @export
summary.aftosmac <- function(object,...){
  z <- object
  if (!is.aftosmac(z)) stop("Most be aftosmac class")
  ans <- z["call"]
  if (is.na(z$covmat[1])) {
    ind.scale <- which("Scale" %in% names(z$coefficients))
    est.smry <- c(z$coefficients[-ind.scale], z$coefficients[ind.scale])
    temp <- list(Estimate = est.smry, StdErr = rep(NA, length(est.smry)),
                 z.value = rep(NA, length(est.smry)),
                 p.value = rep(NA, length(est.smry)))
    res <- list(call = z$call, coefficients = temp, model = z$model)
    class(res) <- "summary.aftosmac"
    return(res)
  }
  if ("Scale" %in% names(z$coefficients)) {
    ind.scale <- which("Scale" %in% names(z$coefficients))
    est.smry <- c(z$coefficients[-ind.scale], z$coefficients[ind.scale])
    ind.scale.se <- which("Scale" %in% colnames(z$covmat))
    se.smry <- as.numeric(c(sqrt(diag(z$covmat))[-ind.scale.se],
                            sqrt(diag(z$covmat))[ind.scale.se]))
  }else {
    est.smry <- z$coefficients
    se.smry <- as.numeric(c(NA, sqrt(diag(z$covmat))))
  }

  z.val.smry <- est.smry / se.smry
  temp <- cbind(Estimate = est.smry, StdErr = se.smry,
                z.value = round(z.val.smry, 3),
                p.value = round(2 * pnorm(-abs(z.val.smry)), 3))

  res <- list(call = z$call, coefficients = temp, model = z$model,
              covg = z$convergence)
  class(res) <- "summary.aftosmac"
  return(res)
}

#' @export
print.summary.aftosmac <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat(as.character(paste0("model: ", x$model)))
  cat("\n")
  if ("1" %in% names(x$covg) | "2" %in% names(x$covg) | !0 %in% x$covg) {
    cat("Failed to get a converging result.")
  }else{
    printCoefmat(as.data.frame(x$coefficients), P.values = TRUE, has.Pvalue = TRUE)
  }
}
