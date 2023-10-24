# get variance matrix
vcovm <- function(DF.Samp, engine, se) {
  if (se == "NULL") {
    return(NA)
  }
  r.Sec <- nrow(DF.Samp)
  g <- aftosmac.est(DF.Samp, engine)
  # estimate vc
  # z <- matrix(rexp(engine@B * r.Sec), nrow = r.Sec, byrow = FALSE)
  # vc <- var(crossprod(z, g)) / r.Sec
  vc <- crossprod(g) / r.Sec
  # inverse of the hessian matrix estimated by the second step subsample
  m_inv <- solve(aftosmac.slope(DF.Samp, engine))
  # sandwich estimator for full estimator
  vx <- m_inv %*% vc %*% t(m_inv) / r.Sec
  out <- vx
  # sandwich estimator for true coe
  if (se == "parTrue") {
    vc_add <- crossprod(DF.Samp$ssps * g, g) / r.Sec * engine@n
    vc_amend <- vc + (r.Sec / engine@n) * vc_add
    vx_amend <- m_inv %*% vc_amend %*% m_inv / r.Sec
    out <- vx_amend
  }
  return(out)
}

# get intercept for semi-parametric models
intcp <- function(DF.Samp, engine) {
  xmat.Samp <- as.matrix(DF.Samp[, -c(1, 2, ncol(DF.Samp))])
  delta.Samp <- DF.Samp[, 2]
  y.Samp <- log(DF.Samp[, 1])
  coe.icpt <- max(eres(
    (y.Samp - xmat.Samp %*% engine@b), delta.Samp, DF.Samp[, ncol(DF.Samp)],
    seq_len(nrow(DF.Samp))
  )[[2]])
  out <- c(coe.icpt, engine@b)
  names(out) <- c("Intercept", colnames(DF.Samp)[-c(1, 2, ncol(DF.Samp))])
  return(out)
}


#' Optimal Subsampling Method for Accelerated Failure Time Models.
#'
#' Estimate the full data statistical inferences based on a subsample. The
#' subsample is derived by sampling with probabilities. The sampling probabilities
#' is defined by the A-optimality and L-optimality from optimal design of
#' experiments.
#'
#' @param formula  a formula expression, of the form \code{response ~ predictors}.
#'     The \code{response} is a \code{Surv} object with right censoring.
#' @param data an optional data.frame in which to interpret the variables occurring
#'     in the \code{formula}.
#' @param control controls maxiter and tolerance.
#' @param subsample.size the pilot and second-subsample size
#' @param sspType the type of optimal subsampling probabilities.
#' @param repeated number of subsmaples needed to derive the final estimator
#' @param R the number of noises used in the resampling approach.
#' @param se method to calculate the standard errors.
#' @param method methods to fit the regression model.
#' @param combine how to combine the pilot and second-step subsample.
#' @param constrasts ...
#'
#' @return The function return a list containing at least the following components:
#' \describe{
#'   \item{coe}{a vector of point estimates}
#'   \item{std}{a vector of estimated standard errors}
#'   \item{converge}{indicator of convergence}
#'   \item{iter}{iteration needed to get the converging result}

#' }
#'
#' @export
#' @keywords aftosmac
#'

aftosmac <- function(formula, data, subsample.size,
                     sspType = c("optA", "optL", "uniform"),
                     method = c("weibull", "ls", "gehan"),
                     se = c("NULL", "parTrue", "parFull"),
                     combine = c("estimator", "sample"),
                     repeated = 0, R = 0,
                     constrasts = NULL,
                     control = list()) {
  method <- match.arg(method)
  se <- match.arg(se)
  sspType <- match.arg(sspType)
  scall <- match.call() # call function, including all inputs
  mnames <- c("", "formula", "data") # variable related to data
  cnames <- names(scall) # names of input arguments
  cnames <- cnames[match(mnames, cnames, 0)] # common arguments of inputs and data
  mcall <- scall[cnames] # input arguments related to model
  mcall[[1]] <- as.name("model.frame")
  m <- eval(mcall, parent.frame()) # make the input data as a dataframe
  mterms <- attr(m, "terms") # attributes of input data
  obj <- unclass(m[, 1]) # survival time and censoring indicator
  if (!inherits(m[[1]], "Surv") || ncol(obj) > 2)
    stop("aftosmac only supports Surv object with right censoring.", call. = FALSE)
  formula[[2]] <- NULL
  DF <- as.data.frame(cbind(obj, model.matrix(mterms, m))) # add covariates
  r0 <- subsample.size[1]
  r <- subsample.size[2]
  if (method == "ls") {
    method <- "semi.ls"
  }else if (method == "gehan") {
    # delete intercept for rank-based approach
    DF <- DF[,-which(colnames(DF) == "(Intercept)")]
    method <- "semi.rank.gehan.s"
    combine <- "estimator"
  }else if (method == "weibull") {
    method <- "par.weibull"
  }
  # create engine
  engine.control <- control[names(control) %in% names(attr(getClass(method), "slots"))]
  engine <- do.call("new", c(list(Class = method), engine.control))
  if (engine@b0 == 0) {
    if (method == "par.weibull") {
      engine@b0 <- c(1, as.numeric(rep(0, ncol(DF)-2)))
    }else if (method == "semi.rank.gehan.s") {
      engine@b0 <- as.numeric(lsfit(DF[, -(1:2)], DF[, 1])$coefficient)[-1]
    }else if (method == "semi.ls") {
      engine@b0 <- as.numeric(rep(0, ncol(DF)-2))
    }
  }
  if (length(engine@b0) != ncol(DF)-2 && method != "par.weibull") {
    stop("Initial value length does not match with the numbers of covariates",
         call. = FALSE)
  }
  if (length(engine@b0) != ncol(DF)-1 && method == "par.weibull") {
    stop("Initial value length does not match with the numbers of covariates",
         call. = FALSE)
  }

  # get optimal SSPs
  engine@r0 <- r0
  engine@n <- nrow(DF)
  optSSPs <- aftosmac.ssps(DF, engine, sspType = sspType)
  if (optSSPs$converge != 0) {
    stop(paste0("Fail to get a converging pilot estimator. The converging code is ",
                optSSPs$converge))
  }
  DF$ssps <- optSSPs$ssp
  indPt <- optSSPs$ind.pt

  # sample with replacement
  indSec <- sample(engine@n, r, prob = optSSPs$ssp, replace = TRUE)

  # combined subsample
  indSamp <- c(indSec, indPt)
  sspSamp <- c(optSSPs$ssp[indSec], rep(1/engine@n, r0))
  DF.Samp <- DF[indSamp, ]
  DF.Samp$ssps <- sspSamp

  if (combine == "sample") {
    est.out <- aftosmac.fit(DF.Samp, engine)
    coe.out <- est.out$coe
    itr.out <- est.out$iter
    if (est.out$converge != 0) {
      stop(paste0("Fail to get a converging second-step estimator.
                  The converging code is ", est.out$converge))
    }
  }else if (combine == "estimator") {
    DF.snd <- DF[indSec, ]
    est.snd <- aftosmac.fit(DF.snd, engine)
    if (est.snd$converge != 0) {
      stop(paste0("Fail to get a converging subsample estimator.
                  The converging code is ", est.snd$converge))
    }
    engine@b <- est.snd$coe
    engine@ind_sub <- seq_len(r)
    M.snd <- r * aftosmac.slope(DF.snd, engine)
    if (method == "semi.ls") {
      engine@b <- est.snd$coe[-1]
      optSSPs$coe.pt <- optSSPs$coe.pt[-1]
    }
    coe.out <- drop(solve(r0 * optSSPs$M.pt + M.snd) %*%
                     (r0 * optSSPs$M.pt %*% optSSPs$coe.pt + M.snd %*% engine@b))
    if (method == "semi.ls") {
      coe.out <- intcp(DF.Samp[, -which(colnames(DF.Samp) == "(Intercept)")], engine)
    }
    itr.out <- est.snd$iter
  }

  # update engine
  engine@b <- coe.out
  engine@ind_sub <- seq_len(r + r0)

  # get variance matrix
  covmat <- vcovm(DF.Samp, engine, se)

  if (method == "par.weibull") {
    names(coe.out) <- c("Scale", colnames(DF)[-c(1, 2, ncol(DF))])
    covmat <- covmat[-2, -2]
    colnames(covmat) <- rownames(covmat) <- names(coe.out)[-2]
  }else if (method == "semi.ls") {
    names(coe.out) <- colnames(DF)[-c(1, 2, ncol(DF))]
    if (!is.na(covmat[1])){
      colnames(covmat) <- rownames(covmat) <- names(coe.out)[-1]
    }
  }else if (method == "semi.rank.gehan.s") {
    # add intercept term for the rank based estimator
    coe.out <- intcp(DF.Samp, engine)
    names(coe.out) <- c("(Intercept)", colnames(DF)[-c(1, 2, ncol(DF))])
    if (!is.na(covmat[1])){
      colnames(covmat) <- rownames(covmat) <- names(coe.out)[-1]
    }
  }

  out <- list(call = scall, vari.name = names(coe.out),
              coefficients = coe.out, covmat = covmat, convergence = 0,
              var.meth = se, ssp.type = sspType, model = method,
              combine = combine, repeated = repeated, R = R,
              ind.ptsample = indPt, ind.sndsample = indSec)
  out$x <- DF[-c(1, 2, ncol(DF))]
  out$y <- DF[, c(1, 2)]
  class(out) <- "aftosmac"
  return(out)
}


