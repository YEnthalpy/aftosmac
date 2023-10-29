parllk.weibull <- function(DF, engine) {
  xmat <- as.matrix(DF[, -c(1, 2, ncol(DF))])
  y <- log(DF$time)
  ssps <- DF$ssps
  er <- drop((y - xmat %*% engine@b[-1]) / engine@b[1])
  exper <- exp(er)
  sum((exper - DF$status * (er - log(engine@b[1]))) / ssps) / engine@n / nrow(DF)
}


parFit.weibull <- function(DF, engine) {
  xmat <- as.matrix(DF[, -c(1, 2, ncol(DF))])
  y <- log(DF$time)
  ssps <- DF$ssps
  beta <- engine@b0[-1]
  sigma <- engine@b0[1]
  for (i in seq_len(engine@maxit)) {
    engine@b <- c(sigma, beta)
    llk.old <- parllk.weibull(DF, engine)
    # update beta
    er <- drop((y - xmat %*% beta) / sigma)
    exper <- exp(er)
    # first derivative with respect to beta
    d.beta <- -colSums(xmat * (exper - DF$status) / ssps / sigma) / engine@n / nrow(DF)
    # second derivative with respec to beta
    d2.beta <- t(xmat * exper / ssps) %*% xmat / (sigma^2) / engine@n / nrow(DF)
    updBeta <- tryCatch(
      solve(d2.beta, -d.beta),
      error = function(e) NA,
      warning = function(w) NA
    )
    if (is.na(updBeta[1])) {
      return(list(coe = rep(NA, ncol(xmat)+1),
                  converge = 1, iter = NA))
    }
    engine@b <- c(sigma, beta + updBeta)
    # backtracking linear search
    llk.new <- parllk.weibull(DF, engine)
    t <- 1
    norm.dbeta <- sum(d.beta ^ 2) / 2
    while (llk.new > llk.old + t * norm.dbeta) {
      t <- t / 2
      engine@b <- c(sigma, beta + t * updBeta)
      llk.new <- parllk.weibull(DF, engine)
    }
    beta <- beta + t * updBeta

    engine@b <- c(sigma, beta)
    llk.old <- parllk.weibull(DF, engine)
    # update sigma
    er <- drop((y - xmat %*% beta) / sigma)
    exper <- exp(er)
    # first derivative with respect to sigma
    d.sig <- sum((er * exper - DF$status * er - DF$status) / ssps / sigma) / engine@n / nrow(DF)
    # second derivative with respect to sigma
    d2.sig <- sum((DF$status + 2 * er * (DF$status - exper) - (er^2)
                   * exper) / ssps) / (sigma^2) / engine@n / nrow(DF)
    updSig <- -d.sig / d2.sig
    engine@b <- c(sigma + updSig, beta)
    # backtracking linear search
    llk.new <- parllk.weibull(DF, engine)
    t <- 1
    norm.dsig <- d.sig ^ 2 / 2
    while (llk.new > llk.old + t * norm.dsig) {
      t <- t / 2
      engine@b <- c(sigma + t * updSig, beta)
      llk.new <- parllk.weibull(DF, engine)
    }
    sigma <- sigma + t * updSig

    if (sqrt(sum(updSig^2 + updBeta^2)) <= engine@tol) {
      return(list(coe = c(sigma, beta), converge = 0, iter = i))
    }
    if (i == engine@maxit) {
      return(list(coe = c(sigma, beta), converge = 2, iter = i))
    }
  }
}

lsFit <- function(DF, engine) {
  xmat <- as.matrix(DF[, -c(1, 2, ncol(DF))])
  y <- log(DF$time)
  beta <- engine@b0
  r <- length(y)
  cxmat <- center(xmat[, -1], DF$ssps, engine@n)
  for (i in seq_len(engine@maxit)) {
    py <- xmat %*% beta
    tmp <- eres((y - py), DF$status, DF$ssps, seq_len(r))[[1]]
    hy <- DF$status * y + (1 - DF$status) * (tmp + py)
    cy <- drop(hy - mean(hy / DF$ssps) / engine@n)
    newBeta <- tryCatch(
      solve(crossprod(cxmat, cxmat / DF$ssps), colSums(cxmat * cy / DF$ssps)),
      error = function(e) NA,
      warning = function(w) NA
    )
    if (is.na(newBeta[1])) {
      return(list(coe = rep(NA, ncol(xmat)), converge = 1, iter = NA))
    }
    newBeta <- c(NA, drop(newBeta))
    newBeta[1] <- max(eres((y - xmat[, -1] %*% newBeta[-1]),
                           DF$status, DF$ssps, seq_len(r))[[2]])
    e <- sqrt(sum((newBeta - beta)^2))
    if (e < engine@tol) {
      return(list(coe = newBeta, converge = 0, iter = i))
    } else {
      beta <- newBeta
    }
    if (i == engine@maxit) {
      return(list(coe = newBeta, converge = 2, iter = i))
    }
  }
}

rankFit.gehan.s <- function(DF, engine) {
  xmat <- as.matrix(DF[, -c(1:2, ncol(DF))])
  y <- log(DF$time)
  init <- lsfit(xmat, y, intercept = FALSE)$coefficient
  out <- nleqslv(
    x = init, fn = function(b) {
      colSums(gehan_smth(xmat, log(DF$time), DF$status, DF$ssps, b, engine@n)) / nrow(xmat)
    }, jac = function(b) {
      gehan_s_jaco(xmat, log(DF$time), DF$status, DF$ssps, b, engine@n)
    }, method = "Broyden", jacobian = FALSE,
    control = list(ftol = engine@tol, xtol = 1e-20,
                   maxit = engine@maxit)
  )
  conv <- out$termcd
  coe <- as.vector(out$x)
  if (conv == 1) {
    conv <- 0
  } else if (conv %in% c(2, 4)) {
    conv <- 2
    coe <- rep(NA, ncol(xmat))
  } else {
    conv <- 1
    coe <- rep(NA, ncol(xmat))
  }
  return(list(coe = coe, converge = conv, iter = out$nfcnt))
}

rankFit.gehan.ns <- function(DF, engine) {
  xmat <- as.matrix(DF[, -c(1:2, ncol(DF))])
  out <- nleqslv(
    x = engine@b0, fn = function(b) {
      colSums(gehan_ns(xmat, log(DF$time),
                       DF$status, b, DF$ssps, engine@n))
    }, method = "Broyden", jacobian = FALSE,
    control = list(ftol = engine@tol, xtol = 1e-20,
                   maxit = engine@maxit)
  )
  conv <- out$termcd
  coe <- out$x
  if (conv == 1) {
    conv <- 0
  } else if (conv %in% c(2, 4)) {
    conv <- 2
    coe <- rep(NA, ncol(xmat))
  } else {
    conv <- 1
    coe <- rep(NA, ncol(xmat))
  }
  return(list(coe = coe, converge = conv, iter = out$nfcnt))
}


## Generic function -- subsample estimates
setGeneric("aftosmac.fit", function(DF, engine) standardGeneric("aftosmac.fit"))

setMethod("aftosmac.fit", signature(engine = "semi.ls"), lsFit)
setMethod("aftosmac.fit", signature(engine = "semi.rank.gehan.s"), rankFit.gehan.s)
setMethod("aftosmac.fit", signature(engine = "semi.rank.gehan.ns"), rankFit.gehan.ns)
setMethod("aftosmac.fit", signature(engine = "par.weibull"), parFit.weibull)
