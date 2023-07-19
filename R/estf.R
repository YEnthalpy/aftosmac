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
  xmat <- as.matrix(DF[, -c(1:2, ncol(DF))])
  y <- log(DF$time)
  cxmat <- center(xmat[, -1], DF$ssps, engine@n)
  py <- xmat %*% engine@b
  tmp <- eres((y - py), DF$status, DF$ssps, engine@ind_sub)[[1]]
  hy <- DF$status * y + (1 - DF$status) * (tmp + py)
  cy <- hy - mean(hy / DF$ssps) / engine@n
  cxmat * drop(cy - cxmat %*% engine@b[-1]) / DF$ssps / engine@n
}

rankEst.gehan.s <- function(DF, engine) {
  gehan_s_mtg(as.matrix(DF[, -c(1:2, ncol(DF))]), log(DF$time), DF$status,
              DF$ssps, engine@b, engine@ind_sub-1, engine@n)
}

rankEst.gehan.ns <- function(DF, engine) {
  gehan_ns_mtg(as.matrix(DF[, -c(1:2, ncol(DF))]), log(DF$time),
               DF$status, DF$ssps, engine@b, engine@ind_sub-1, engine@n)
}

parEst.weibull <- function(DF, engine) {
  xmat <- as.matrix(DF[, -c(1:2, ncol(DF))])
  er <- drop((log(DF$time) - xmat %*% engine@b[-1]) / engine@b[1])
  exper <- exp(er)
  d.beta <- xmat * (exper - DF$status) / DF$ssps / engine@b[1]
  d.sig <- (er * exper - DF$status * er - DF$status) / DF$ssps / engine@b[1]
  cbind(d.sig, d.beta)
}



