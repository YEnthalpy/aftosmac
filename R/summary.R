is.aftosmac <- function(x) {
  inherits(x, "aftosmac")
}

#' @export
summary.aftosmac <- function(object,...){
  z <- object
  if (!is.aftosmac(z)) stop("Most be aftosmac class")
  ans <- z["call"]
  var.meth <- z$var.meth
  
  if (z$model == "par.weibull") {
    est.smry <- c(z$coefficients[-1], z$coefficients[1])
    se.smry <- as.numeric(c(NA, sqrt(diag(z$covmat))[-1], sqrt(diag(z$covmat))[1]))
  }
  
  z.val.smry <- est.smry / se.smry
  temp <- cbind(Estimate = est.smry, StdErr = se.smry, 
                z.value = round(z.val.smry, 3), 
                p.value = round(2 * pnorm(-abs(z.val.smry)), 3))
  
  res <- list(call = z$call, coefficients = temp, model = z$model)
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
  printCoefmat(as.data.frame(x$coefficients), P.values = TRUE, has.Pvalue = TRUE)
}
