aftosmac.ssps <- function(DF, engine, sspType = c("uniform", "optA", "optL")) {
  xmat <- as.matrix(DF[, c(1, 2)])
  y <- log(DF$time)
  delta <- DF$status
  DF$ssps <- rep(1/engine@n, engine@n)
  if (sspType == "uniform") {
    return(list(ssp = rep(1/engine@n, engine@n),
                ind.pt = sample(engine@n, engine@r0, TRUE),
                converge = 0, iteration = 0))
  }else {
    ind_pt <- sample(engine@n, engine@r0, TRUE)
    mle_pt <- aftosmac.fit(DF = DF[ind_pt, ], engine = engine)
    if (mle_pt$converge != 0) {
      return(list(coe.pt = NA, ssp = NA, ind.pt = NA,
                  converge = mle_pt$converge))
    }
    # contribution of each observation to the estimating function
    engine@ind_sub <- ind_pt
    engine@b <- mle_pt$coe
    g <- aftosmac.est(DF = DF, engine = engine)
    engine@ind_sub <- seq_len(engine@r0)
    m <- aftosmac.slope(DF = DF[ind_pt, ], engine = engine)
    # optimal ssps
    if (sspType == "optL") {
      g_nm <- sqrt(rowSums(g^2))
      ssp <- g_nm / sum(g_nm) * (1 - engine@alpha) + engine@alpha / engine@n
    } else if (sspType == "optA") {
      m_inv <- solve(m)
      m_mse <- sqrt(colSums((tcrossprod(m_inv, g))^2))
      ssp <- m_mse / sum(m_mse) * (1 - engine@alpha) + engine@alpha / engine@n
    }
    return(list(ssp = ssp, M.pt = m, ind.pt = ind_pt,
                converge = 0, coe.pt = mle_pt$coe, iteration = mle_pt$iter))
  }
}
