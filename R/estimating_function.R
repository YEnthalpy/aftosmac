eres <- function(e, delta, pi, ind_km) {
  e_sub <- e[ind_km]
  ord_sub <- order(e_sub)
  km_sub <- km(e_sub[ord_sub], delta[ind_km][ord_sub], pi[ind_km][ord_sub])
  es_sub <- km_sub[[1]]
  s_sub <- km_sub[[2]]
  edif <- c(diff(es_sub), 0)
  int <- rev(cumsum(rev(edif * s_sub)))
  es_int <- int + s_sub * es_sub
  int <- int / s_sub + es_sub
  ehat <- approx(es_sub, int, e, method = "constant", ties = "ordered")$y
  ehat[is.na(ehat)] <- e[is.na(ehat)]
  return(list(ehat, es_int))
}

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



