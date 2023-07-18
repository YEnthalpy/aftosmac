lsEst <- function(DF, engine) {
  xmat <- as.matrix(DF$covaraites)
  y <- log(DF$time)
  delta <- DF$delta
  ssps <- DF$ssps
  beta <- engine@b
  n <- engine@n
  cxmat <- center(xmat[, -1], ssps, n)
  py <- xmat %*% beta
  tmp <- eres((y - py), delta, ssps, engine@ind_sub)[[1]]
  hy <- delta * y + (1 - delta) * (tmp + py)
  cy <- hy - mean(hy / ssps) / n
  beta <- beta[-1]
  cxmat * drop(cy - cxmat %*% beta) / ssps / n
}

rankEst.gehan.s <- function(DF, engine) {
  xmat <- as.matrix(DF$covaraites)
  y <- log(DF$time)
  delta <- DF$delta
  ssps <- DF$ssps
  gehan_s_mtg(xmat, y, delta, ssps, engine@b, engine@ind_sub-1, engine@n)
}

rankEst.gehan.ns <- function(DF, engine) {
  xmat <- as.matrix(DF$covaraites)
  y <- log(DF$time)
  delta <- DF$delta
  ssps <- DF$ssps
  gehan_ns_mtg(xmat, y, delta, ssps, engine@b, engine@ind_sub-1, engine@n)
}

parEst.weibull <- function(DF, engine) {
  xmat <- as.matrix(DF$covaraites)
  y <- log(DF$time)
  delta <- DF$delta
  ssps <- DF$ssps
  sigma <- engine@b[1]
  er <- drop((y - xmat %*% engine@b[-1]) / sigma)
  exper <- exp(er)
  d.beta <- xmat * (exper - delta) / ssps / sigma
  d.sig <- (er * exper - delta * er - delta) / ssps / sigma
  cbind(d.sig, d.beta)
}



