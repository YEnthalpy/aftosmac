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
  DF <- as.data.frame(cbind(obj, model.matrix(mterms, m, contrasts))) # add covariates
  # intercept is added in method functions if needed
  DF <- DF[, -which(colnames(DF) == "(Intercept)")]
  if (fitMtd == "ls") {
    method <- "semi.ls"
  }else if (fit.Mtd == "rank") {
    method <- paste("semi", fitMtd, rankWt, eqType, sep = ".")
  }else if (fit.Mtd == "par") {
    method <- paste(fit.Mtd, parDist, sep = ".")
  }
  # create engine
  engine.control <- control[names(control) %in% names(attr(getClass(method), "slots"))]
  engine <- do.call("new", c(list(Class = method), engine.control))
  if (engine@b0 == 0) {
    engine@b0 <- as.numeric(rep(0, ncol(DF)-2))
  }else if (engine@b0 == 1) {
    engine@b0 <- as.numeric(lsfit(DF[, -(1:2)], DF[, 1])$coefficient)[-1]
  }
  if (length(engine@b0) != ncol(DF) - 2)
    stop("Initial value length does not match with the numbers of covariates",
         call. = FALSE)
  # get optimal SSPs
  engine@r0 <- r0
  engine@n <- nrow(DF)
  optSSPs <- aftosmac.ssps(DF, engine, fitMtd = fitMtd,
                           rankWt = rankWt, sspType = sspType)
  if (optSSPs$converge != 0) {
    stop(paste0("Fail to get a converging pilot estimator. The converging code is ", ssps$converge))
  }
  indPt <- optSSPs$ind.pt
  # sample with replacement
  indSec <- sample(engine@n, r, prob = optSSPs$ssp, replace = TRUE)
  if (estMtd == "cs") {
    indSamp <- c(indPt, indSec)
    sspSamp <- c(rep(1/engine@n, r0), optSSPs$ssp[indSec])
    DF.Samp <- DF[indSamp, ]
    DF.Samp$ssps <- sspSamp
    estOut <- aftosmac.fit(DF.Samp, engine)
    if (estOut$converge != 0) {
      stop(paste0("Fail to get a converging second-step estimator.
                  The converging code is ", estOut$converge))
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
    M.pt <- r0 * aftosmac.slope(DF[indPt], engine) # wrong
    estOut <- drop(solve(M.pt + M.snd) %*% (M.pt %*% optSSPs$coe.pt + M.snd %*% est.snd$coe))
  }
  if (!is.NULL(se)) {
    indSamp <- c(indPt, indSec)
    sspSamp <- c(rep(1/engine@n, r0), optSSPs$ssp[indSec])
    engine@b <- estOut
    engine@ind_sub <- seq_len(r + r0)
    g <- aftosmac.est(DF[indSamp, ], engine)
    vc <- crossprod(g) * (r + r0)
    vc_add <- crossprod(sspSamp * g, g) * (r + r0) * engine@n
    vc_amend <- vc + ((r + r0) / engine@n) * vc_add
    # inverse of the hessian matrix estimated by the second step subsample
    m_inv <- solve(aftosmac.slope(DF[indSamp, ], engine))
    # sandwich estimator for full estimator
    vx <- m_inv %*% vc %*% m_inv / (r + r0)
    std <- c(sqrt(diag(vx)))
    # sandwich estimator for true coe
    vx_amend <- m_inv %*% vc_amend %*% m_inv / (r + r0)
    std_amend <- c(sqrt(diag(vx_amend)))
  }
  return(list(
    coe = coe_out, std = std, std_amend = std_amend,
    iter = iter, converge = 0, time = time
  ))
}


