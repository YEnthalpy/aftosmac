lsSlp <- function(DF, engine) {

}

rankSlp.gehan.s <- function(DF, engine) {
  xmat <- as.matrix(DF$covaraites)
  y <- log(DF$time)
  delta <- DF$delta
  ssps <- DF$ssps
  gehan_s_jaco(xmat, y, delta, ssps, engine@b, engine@n)
}

parSlp.weibull <- function(DF, engine) {
  xmat <- as.matrix(DF$covaraites)
  y <- log(DF$time)
  delta <- DF$delta
  ssps <- DF$ssps
  sigma <- engine@b[1]
  er <- drop((y - xmat %*% engine@b[-1]) / sigma)
  exper <- exp(er)
  d2.beta <- -t(xmat * exper  / ssps) %*% xmat / (sigma^2)
  d2.sig <- sum((delta + 2 * er * (delta - exper) - (er^2)
                 * exper) / ssps) / (sigma^2)
  d2.sigbeta <- colSums((xmat * (delta - exper * (1 + er)))) / (sigma ^ 2)
  cbind(c(d2.sig, d2.sigbeta), rbind(t(d2.sigbeta), d2.beta))
}

estSlp <- function(DF, engine, fitMtd = c("rank", "ls"),
                   rankWt = c("gehan")) {
  xmat <- as.matrix(DF$covaraites)
  y <- log(DF$time)
  delta <- DF$delta
  ssps <- DF$ssps
  p <- ncol(xmat) - 1
  r <- length(y)
  z <- matrix(rexp(engine@B * r), nrow = r, byrow = FALSE)
  fitMtd <- match.arg(fitMtd)
  rankWt <- match.arg(rankWt)
  if (fitMtd == "ls") {
    method <- "semi.ls"
  }else if (fit.Mtd == "rank") {
    method <- paste("semi", fitMtd, rankWt, "s", sep = ".")
  }
  engine <- do.call("new", list(Class = method))

  g <- uls(x, y, delta, beta, pi, n, seq_along(y))
  v <- var(crossprod(z, g)) / r
  zbs <- matrix(ifelse(rbinom(b * p, 1, 0.5) == 1, 1, -1),
                ncol = p, byrow = TRUE
  )
  zb <- zbs %*% with(eigen(v), vectors %*% (values^(-0.5) * t(vectors)))
  beta1 <- zb / sqrt(r) + matrix(rep(beta[-1], b),
                                 ncol = p,
                                 byrow = TRUE
  )
  beta1 <- cbind(rep(beta[1], b), beta1)
  response <- matrix(NA, nrow = b, ncol = p)
  for (i in seq_len(b)) {
    response[i, ] <- colSums(uls(
      x, y, delta, beta1[i, ],
      pi, n, seq_along(y)
    )) / sqrt(r)
  }
  return(t(lsfit(zb, response, intercept = FALSE)$coefficients))
}
