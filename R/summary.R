is.aftosmac <- function(x) {
  inherits(x, "aftosmac")
}

#' @export
summary.aftosmac <- function(object,...){
  z <- object
  if (!is.aftosmac(z)) stop("Most be aftosmac class")
  ans <- z["call"]
  if ("Scale" %in% names(z$coefficients)) {
    ind.scale <- which("Scale" %in% names(z$coefficients))
    est.smry <- c(z$coefficients[-ind.scale], z$coefficients[ind.scale])
    ind.scale.se <- which("Scale" %in% colnames(z$covmat[[1]]))
    se.smry <- lapply(z$covmat, function(t){
      as.numeric(c(sqrt(diag(t))[-ind.scale.se],
                   sqrt(diag(t))[ind.scale.se]))
    })
  }else {
    est.smry <- z$coefficients
    se.smry <- lapply(z$covmat, function(t){
      as.numeric(c(NA, sqrt(diag(t))))
    })
  }
  temp <- lapply(se.smry, function(t){
    z.val.smry <- est.smry / t
    out <- cbind(Estimate = est.smry, StdErr = t,
                  z.value = round(z.val.smry, 3),
                  p.value = round(2 * pnorm(-abs(z.val.smry)), 3))
    return(out)
  })
  res <- list(call = z$call, coefficients = temp, method = z$method,
              covg = z$convergence)
  class(res) <- "summary.aftosmac"
  return(res)
}

#' @export
print.summary.aftosmac <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat(as.character(paste0("method: ", x$method)))
  cat("\n")
  if ("1" %in% names(x$covg) | "2" %in% names(x$covg) | !0 %in% x$covg) {
    cat("Failed to get a converging result.")
  }else{
    for (i in seq_along(x$coefficients)) {
      cat(as.character(paste0("Statistical Inferences: ", names(x$coefficients)[i])))
      cat("\n")
      printCoefmat(as.data.frame(x$coefficients[[i]]), P.values = TRUE, has.Pvalue = TRUE)
    }
  }
}
