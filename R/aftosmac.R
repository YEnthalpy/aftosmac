aftosmac <- function(formula, data, r0, r, sspType, B = 1, R = 20,
                     rankWt = c("gehan"), parDist = c("weibull"),
                     eqType = c("s", "ns"),
                     fitMtd = c("rank", "ls", "par"),
                     estMtd = c("cs", "ce"),
                     se = c("NULL", "parTrue", "parFull"),
                     control = list()) {
  rankWt <- match.arg(rankWt)
  eqType <- match.arg(eqType)
  fitMtd <- match.arg(fitMtd)
  estMtd <- match.arg(estMtd)
  parDist <- match.arg(parDist)
  se <- match.arg(se)
  scall <- match.call() # call function, including all inputs
  mnames <- c("", "formula", "data") # variable related to data
  cnames <- names(scall) # names of input arguments
  cnames <- cnames[match(mnames, cnames, 0)] # common arguments of inputs and data
  mcall <- scall[cnames] # input arguments related to model
  mcall[[1]] <- as.name("model.frame")
  m <- eval(mcall, parent.frame()) # make the input data as a dataframe
  mterms <- attr(m, "terms") # attributes of input data
  obj <- unclass(m[, 1]) # survival time and censoring indicator
  if (class(m[[1]]) != "Surv" || ncol(obj) > 2)
    stop("aftsrr only supports Surv object with right censoring.", call. = FALSE)
  formula[[2]] <- NULL
  DF <- as.data.frame(cbind(obj, model.matrix(mterms, m))) # add covariates
  if (fitMtd == "ls") {
    method <- "semi.ls"
  }else if (fitMtd == "rank") {
    # delete intercept for rank-based approach
    DF <- DF[,-which(colnames(DF) == "(Intercept)")]
    method <- paste("semi", fitMtd, rankWt, eqType, sep = ".")
  }else if (fitMtd == "par") {
    method <- paste(fitMtd, parDist, sep = ".")
  }
  # create engine
  engine.control <- control[names(control) %in% names(attr(getClass(method), "slots"))]
  engine <- do.call("new", c(list(Class = method), engine.control))
  if (engine@b0 == 0) {
    engine@b0 <- as.numeric(rep(0, ncol(DF)-2))
  }else if (engine@b0 == 1) {
    if (fitMtd == "rank") {
      engine@b0 <- as.numeric(lsfit(DF[, -(1:2)], DF[, 1])$coefficient)[-1]
    }else {
      engine@b0 <- as.numeric(lsfit(DF[, -(1:2)], DF[, 1], intercept = FALSE)$coefficient)
    }
  }
  if (length(engine@b0) != ncol(DF)-2) {
    stop("Initial value length does not match with the numbers of covariates",
         call. = FALSE)
  }
  if (fitMtd == "par") {
    engine@b0 <- c(1, engine@b0)
  }
  # get optimal SSPs
  engine@r0 <- r0
  engine@n <- nrow(DF)
  optSSPs <- aftosmac.ssps(DF, engine, fitMtd = fitMtd,
                           rankWt = rankWt, sspType = sspType)
  if (optSSPs$converge != 0) {
    stop(paste0("Fail to get a converging pilot estimator. The converging code is ",
                optSSPs$converge))
  }
  DF$ssps <- optSSPs$ssp
  indPt <- optSSPs$ind.pt
  # sample with replacement
  indSec <- sample(engine@n, r, prob = optSSPs$ssp, replace = TRUE)
  if (estMtd == "cs") {
    indSamp <- c(indPt, indSec)
    sspSamp <- c(rep(1/engine@n, r0), optSSPs$ssp[indSec])
    DF.Samp <- DF[indSamp, ]
    DF.Samp$ssps <- sspSamp
    est.out <- aftosmac.fit(DF.Samp, engine)
    coe.out <- est.out$coe
    itr.out <- est.out$iter
    if (est.out$converge != 0) {
      stop(paste0("Fail to get a converging second-step estimator.
                  The converging code is ", est.out$converge))
    }
  }else if (estMtd == "ce") {
    DF.snd <- DF[indSec, ]
    est.snd <- aftosmac.fit(DF.snd, engine)
    if (est.snd$converge != 0) {
      stop(paste0("Fail to get a converging subsample estimator.
                  The converging code is ", est.snd$converge))
    }
    engine@b <- est.snd$coe
    M.snd <- r * aftosmac.slope(DF.snd, engine)
    engine@b <- optSSPs$coe.pt
    coe.out <- drop(solve(r0 * optSSPs$M.pt + M.snd) %*%
                     (r0 * optSSPs$M.pt %*% optSSPs$coe.pt + M.snd %*% est.snd$coe))
    itr.out <- est.snd$iter
  }
  if (se != "NULL") {
    indSamp <- c(indPt, indSec)
    sspSamp <- c(rep(1/engine@n, r0), optSSPs$ssp[indSec])
    engine@b <- coe.out
    engine@ind_sub <- seq_len(r + r0)
    g <- aftosmac.est(DF[indSamp, ], engine)
    vc <- crossprod(g) * (r + r0)
    # inverse of the hessian matrix estimated by the second step subsample
    m_inv <- solve(aftosmac.slope(DF[indSamp, ], engine))
    # sandwich estimator for full estimator
    vx <- m_inv %*% vc %*% m_inv / (r + r0)
    std <- c(sqrt(diag(vx)))
    # sandwich estimator for true coe
    if (se == "parTrue") {
      vc_add <- crossprod(sspSamp * g, g) * (r + r0) * engine@n
      vc_amend <- vc + ((r + r0) / engine@n) * vc_add
      vx_amend <- m_inv %*% vc_amend %*% m_inv / (r + r0)
      std <- c(sqrt(diag(vx_amend)))
    }
    names(coe.out) <- names(std) <- colnames(DF[, -c(1, 2, ncol(DF))])
  }else {
    names(coe.out) <- colnames(DF[, -c(1, 2, ncol(DF))])
    std <- NA
  }
  return(list(coe = coe.out, std = std, converge = 0, iter = itr.out))
}


