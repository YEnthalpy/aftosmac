aftosmac.ssps <- function(DF, engine, fitMtd = c("rank", "ls"),
                        rankWt = c("gehan"),
                        sspType = c("uniform", "optA", "optL")) {
  xmat <- as.matrix(DF$covaraites)
  y <- log(DF$time)
  delta <- DF$delta
  ssps <- DF$ssps
  n <- length(y)
  r0 <- engine@r0
  if (sspType == "uniform") {
    return(list(ssp = rep(1/n, n), ind_pt = sample(n, r0, TRUE),
                converge = 0))
  }else {
    ssp_pt <- rep(1/n, r0)
    ind_pt <- sample(n, r0, TRUE)
    xmat_pt <- xmat[ind_pt, ]
    y_pt <- y[ind_pt]
    dt_pt <- delta[ind_pt]
    DF_pt <- list(x = xmat_pt, )
    mle_pt <- aftosmac.est(DF, engine)
    if (mle_pt$converge != 0) {
      return(list(bt_pt = NA, ssp = NA, index.pilot = NA,
                  converge = mle_pt$converge))
    }
    bt_pt <- mle_pt$coe
    g <- aftosmac.est(DF, engine)
    if (ssp_type == "optL") {
      g_nm <- sqrt(rowSums(g^2))
      ssp <- g_nm / sum(g_nm) * (1 - engine@alpha) + engine@alpha / n
    } else if (ssp_type == "optA") {
      m_inv <- solve(aftosmac.slope(DF, engine))
      m_mse <- sqrt(colSums((tcrossprod(m_inv, g))^2))
      ssp <- m_mse / sum(m_mse) * (1 - engine@alpha) + engine@alpha / n
    }
    return(list(ssp = ssp, index.pilot = ind_pt, converge = 0))
  }
}
