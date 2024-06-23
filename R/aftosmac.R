# get variance matrix
vcovm <- function(DF.Samp, engine) {
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
  # sandwich estimator for true coe
  vc_add <- crossprod(DF.Samp$ssps * g, g) / r.Sec * engine@n
  vc_amend <- vc + (r.Sec / engine@n) * vc_add
  vx_amend <- m_inv %*% vc_amend %*% m_inv / r.Sec
  out <- list("Full Data Estimates" = vx, "True Regression Coefficients" = vx_amend)
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
  names(out) <- c("(Intercept)", colnames(DF.Samp)[-c(1, 2, ncol(DF.Samp))])
  return(out)
}

# Function to get one subsample estimator
onefit <- function(DF, engine, optSSPs, combine, method, n.repeat) {
  indPt <- optSSPs$ind.pt
  # sample with replacement
  indSec <- sample(engine@n, engine@r, prob = optSSPs$ssp, replace = TRUE)

  # combined subsample
  indSamp <- c(indSec, indPt)
  sspSamp <- c(optSSPs$ssp[indSec], rep(1/engine@n, engine@r0))
  DF.Samp <- DF[indSamp, ]
  DF.Samp$ssps <- sspSamp

  if (combine == "sample") {
    est.out <- aftosmac.fit(DF.Samp, engine)
    covg.out <- est.out$converge
    if (covg.out != 0) {
      return(list(coe = NA, itr = NA,
                  covg = covg.out, DF.Samp = NA))
    }
    coe.out <- est.out$coe
    itr.out <- est.out$iter
  }else if (combine == "estimator") {
    DF.snd <- DF[indSec, ]
    est.snd <- aftosmac.fit(DF.snd, engine)
    covg.out <- est.snd$converge
    if (covg.out != 0) {
      return(list(coe = NA, itr = NA,
                  covg = covg.out, DF.Samp = NA))
    }
    engine@b <- est.snd$coe
    engine@ind_sub <- seq_len(engine@r)
    M.snd <- engine@r * aftosmac.slope(DF.snd, engine)
    if (method == "semi.ls") {
      engine@b <- est.snd$coe[-1]
      optSSPs$coe.pt <- optSSPs$coe.pt[-1]
    }
    coe.out <- drop(solve(engine@r0 * optSSPs$M.pt + M.snd) %*%
                      (engine@r0 * optSSPs$M.pt %*% optSSPs$coe.pt + M.snd %*% engine@b))
    if (method == "semi.ls") {
      engine@b <- coe.out
      coe.out <- intcp(DF.Samp[, -which(colnames(DF.Samp) == "(Intercept)")], engine)
    }
    itr.out <- est.snd$iter
  }
  if (method == "semi.rank.gehan.s") {
    engine@b <- coe.out
    coe.out <- intcp(DF.Samp, engine)
  }
  if (n.repeat != 1) {
    DF.Samp <- NA
  }
  return(list(coe = coe.out, iter = itr.out, covg = covg.out, DF.Samp = DF.Samp))
}


#' Optimal Subsampling Method for Accelerated Failure Time Models.
#'
#' Fit a accelerated failure time (AFT) model based on a small subsample.
#' The subsample is chosen by sampling with resampling. The subsampling probabilities (SSPs)
#' are derived by the A-optimal and L-optimal criteria from the optimal design of
#' experiment. Three types of models are supported, the Weibull parametric AFT model,
#' the rank-based semi-parametric AFT model and the least-squares based semi-parametric
#' AFT model with the Gehan's weight.
#'
#' The estimating equation is solved by
#' a block gradient decent method for the Weibull AFT model and by an iterative approach
#' for the least-squares based semi-prarmetric AFT model.
#' For the rank-based semi-parametric AFT model, the non-smooth estimating equation
#' is smoothed by an induced smooth procedure which is solved the
#' by the Quasi Newton's method implemented as \code{nleqslv} in the package \pkg{nleqslv}.
#'
#' When \code{n.repeat = 1}, the variance matrix is estimated by a sandwich form
#' \deqn{\Sigma = A^{-1}V(A^{-1})^T,} where \eqn{V} is the asymptotic variance of
#' the estimating function and \eqn{A} is the slope matrix. The sandwich estimator
#' estimates the variance matrix well for the Weibull parametric AFT model and the rank-based
#' semi-prarmetric AFT model. For the least-squares based semi-prarmetric AFT
#' model whose slope matrix
#' is estimated by a resampling method. And the corresponding
#' sandwich estimator will overestimate the variance when the censoring
#' rate is high censoring rate. We fix this specific problem for the least-squares based approach
#' by selecting more than one subsample (\code{n.repeat > 1}) to estimate the variance matrix
#' directly by multiple subsample estimators.
#' The variance matrix with respect to the full data estimator
#' is calculated used these estimators instead of using the sandwich form.
#' Moreover, the variance matrix with respect to the full data estimator and true
#' regression coefficients are both calculated.
#'
#' The optimal SSPs are estimated by a pilot estimator using a small pilot subsample.
#' The subsample estimator is derived by the subsample chosed based on the estimated
#' optimal SSPs. To enhance the estimation accuracy, the pilot estimator and the subsample
#' estimator need to be combined. We provided two methods for combination. When
#' \code{combine = "sample"}, the pilot sample and the subsample sample are combined to
#' derived the final estimator. When \code{combine = "estimator"}, the pilot estimator
#' and the subsample estimator is combined to derived the final estimator. The second
#' method is less time-consuming and is always used in the divide-and-conquer approach.
#'
#'
#'
#' @param formula  a formula expression, of the form \code{response ~ predictors}.
#'     The \code{response} is a \code{Surv} object with right censoring.
#' @param data an data.frame in which to interpret the variables occurring
#'     in the \code{formula}.
#' @param n.pilot the pilot subsample size
#' @param n.sub the subsample size which is always much large than
#'     the pilot subsample size
#' @param n.repeat number of subsmaples used to derive the final estimator. The
#'     default is set to be \code{n.repeat = 1}. When \code{n.repeat > 1}, the
#'     variance matrix is not estimated by the sandwich estimator.
#'     User is suggested to use multiple subsamples (\code{n.repeat > 1}) only
#'     for the least-squares based semi-prarmetric AFT model whose slope matrix
#'     is estimated by a resampling method. And the corresponding
#'     sandwich estimator will overestimate the variance when the censoring
#'     rate is high censoring rate. Moreover, when \code{n.repeat > 1}, only the
#'     variance matrix with respect to the full sample estimator can be estimated.
#'     Moreover, \code{n.repeat} should be larger than the number of covariates \code{p}.
#' @param contrasts an optional list.
#' @param subset subset ot the observations to be used in the fit.
#' @param sspType the type of subsampling probabilities (SSPs). Three types of
#'     SSPs are provided:
#' \describe{
#'   \item{\code{uniform}}{the uniform SSPs.}
#'   \item{\code{optA}}{the A-optimal SSPs, always yields the least mean square errors.}
#'   \item{\code{optL}}{the L-optimal SSPs which is less time-consuming than the A-optimal SSPs.}
#' }
#' @param method method to fit the data:
#' \describe{
#'   \item{\code{weibull}}{the Weibull parametric AFT model.}
#'   \item{\code{ls}}{the least-square based semi-parametric AFT model.}
#'   \item{\code{rank}}{the rank based semi-parametric AFT model which requires much less subsample
#'                       size to get a converging results than the least-squares approach.}
#' }
#' @param combine methods to combine the results by the pilot sample and the subsample:
#' \describe{
#'   \item{\code{sample}}{combine the pilot sample and the subsample to derive the
#'                         final estimator by the combined subsample.}
#'   \item{\code{estimator}}{combine the pilot estimator and the subsample estimator}
#' }
#' @param control controls equation solver, maxiter, tolerance and
#'     resampling variance estimation.



#' @export
#'
#' @return \code{aftosmac} returns an object of class "\code{aftosmac}" representing the fit.
#' An object of class "\code{aftosmac}" is a list containing at least the following components:
#' \describe{
#'   \item{coefficients}{A vector of beta estimates}
#'   \item{covmat}{A matrix of covariance estimates}
#'   \item{convergence}{An integer code indicating type of convergence. If \code{n.repeat > 1},
#'        this object will be a frequency table of all possible convergences types.}
#'   \describe{
#'     \item{0}{indicates successful convergence.}
#'     \item{1}{indicates that the iteration limit \code{maxit} has been reached.}
#'     \item{2}{indicates failure due to stagnation.}
#'   }
#' }
#'
#' @references Yang, Z., Wang, H., Yan, J. (2023) Subsampling Approach for Least Squares
#' Fitting of Semi-parametric Accelerated Failure Time Models to Massive
#' Survival Data.
#' \emph{Research Square preprint}, https://www.researchsquare.com/article/rs-2866415/v1.
#' @references Yang Z., Wang, H., Yan, J. (2022) Optimal Subsampling for
#' Parametric Accelerated Failure Time Models with Massive
#' Survival Data.  2022; 27): 5421–5431.
#' \emph{Statistics in Medicine}, \bold{41}(27): 5421–5431.
#' @references Zeng, D. and Lin, D. Y. (2008)
#' Efficient Resampling Methods for Nonsmooth Estimating Functions.
#' \emph{Biostatistics}, \bold{9}, 355--363
#' @references Chiou, S., Kang, S. and Yan, J. (2014)
#' Fast Accelerated Failure Time Modeling for Case-Cohort Data. \emph{Statistics and Computing},
#' \bold{24}(4): 559--568.
#' @references Johnson, L. M. and Strawderman, R. L. (2009)
#' Induced Smoothing for the Semiparametric Accelerated Failure Time Model:
#' Asymptotic and Extensions to Clustered Data. \emph{Biometrika}, \bold{96}, 577 -- 590.
#'
#' @importFrom stats pnorm printCoefmat
#' @export
#' @example inst/examples/ex_aftosmac.R
#' @keywords aftosmac
#'

aftosmac <- function(formula, data, n.pilot, n.sub, n.repeat = 1,
                     contrasts = NULL, subset,
                     sspType = c("optA", "optL", "uniform"),
                     method = c("parametric", "least_squares", "rank"),
                     combine = c("estimator", "sample"),
                     control = list()) {
  method <- match.arg(method)
  sspType <- match.arg(sspType)
  scall <- match.call() # call function, including all inputs
  mnames <- c("", "formula", "data", "subset") # variable related to data
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
  DF <- as.data.frame(cbind(obj, model.matrix(mterms, m, contrasts))) # add covariates
  if (sspType == "uniform") {
    combine <- "sample"
  }
  mtd.name <- method
  if (method == "least_squares") {
    method <- "semi.ls"
  }else if (method == "rank") {
    # delete intercept for rank-based approach
    DF <- DF[,-which(colnames(DF) == "(Intercept)")]
    method <- "semi.rank.gehan.s"
  }else if (method == "parametric") {
    method <- "par.weibull"
  }
  # create engine
  engine.control <- control[names(control) %in% names(attr(getClass(method), "slots"))]
  engine <- do.call("new", c(list(Class = method), engine.control))
  if (engine@b0 == 0) {
    if (method == "par.weibull") {
      engine@b0 <- c(1, as.numeric(rep(0, ncol(DF)-2)))
    }else if (method == "semi.rank.gehan.s") {
      engine@b0 <- as.numeric(lsfit(DF[, -(1:2)], log(DF[, 1]))$coefficient)[-1]
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
  ptm1 <- proc.time()
  engine@r0 <- n.pilot
  engine@r <- n.sub
  engine@n <- nrow(DF)
  optSSPs <- aftosmac.ssps(DF, engine, sspType = sspType)
  if (optSSPs$converge != 0) {
    tol <- engine@tol
    engine@tol <- engine@tol_pilot
    optSSPs <- aftosmac.ssps(DF, engine, sspType = sspType)
    if (optSSPs$converge != 0) {
      stop(paste0("Fail to get a converging pilot estimator. The converging code is ",
                  optSSPs$converge))
    }
    engine@tol <- tol
  }
  DF$ssps <- optSSPs$ssp
  ptm2 <- proc.time()
  # get subsample estimator
  
  subest <- replicate(n.repeat,
                      list(onefit(DF, engine, optSSPs, combine, method, n.repeat)))
  ptm3 <- proc.time()
  engine@ind_sub <- seq_len(engine@r + engine@r0)
  if (n.repeat == 1) {
    covg.out <- subest[[1]]$covg
    if (covg.out != 0) {
      tm_tmp <- matrix(c(ptm2 - ptm1, rep(NA, 2 * length(ptm2))), 
                       ncol = length(ptm2), byrow = TRUE)
      rownames(tm_tmp) <- c("ssps", "coefficients", "vcov")
      out <- list(call = scall, vari.name = colnames(DF)[-c(1, 2, ncol(DF))],
                  coefficients = NA, covmat = NA,
                  convergence = c("pilot" = 0, "sec_step" = covg.out),
                  n.iteration = c("pilot" = optSSPs$iteration, "sec_step" = NA),
                  time = tm_tmp,
                  ssp.type = sspType, method = mtd.name,
                  combine = combine, n.repeat = n.repeat)
      out$x <- DF[-c(1, 2, ncol(DF))]
      out$y <- DF[, c(1, 2)]
      out$ssp <- DF[, ncol(DF)]
      class(out) <- "aftosmac"
      return(out)
    }
    coe.out <- subest[[1]]$coe
    itr.out <- subest[[1]]$iter
    DF.Samp <- subest[[1]]$DF.Samp
    if (method == "semi.rank.gehan.s") {
      engine@b <- coe.out[-1]
    }else {
      engine@b <- coe.out
    }
    covmat <- vcovm(DF.Samp, engine)
  }else {
    covg.out <- sapply(subest, function(t){t[["covg"]]})
    if (sum(covg.out) != 0) {
      out <- list(call = scall, vari.name = colnames(DF)[-c(1, 2, ncol(DF))],
                  coefficients = NA, covmat = NA, convergence = table(covg.out),
                  n.iteration = NA,
                  time = NA,
                  ssp.type = sspType, method = mtd.name,
                  combine = combine, n.repeat = n.repeat)
      out$x <- DF[-c(1, 2, ncol(DF))]
      out$y <- DF[, c(1, 2)]
      class(out) <- "aftosmac"
      return(out)
    }
    coe.mat <- sapply(subest, function(t){t[["coe"]]})
    coe.out <- rowMeans(coe.mat)
    covmat <- list("Full Data Estimates" = var(t(coe.mat[-1, ])) / n.repeat)
    itr.out <- sapply(subest, function(t){t[["iter"]]})
  }
  ptm4 <- proc.time()
  if (method == "par.weibull") {
    names(coe.out) <- c("Scale", colnames(DF)[-c(1, 2, ncol(DF))])
    covmat <- lapply(covmat, function(t){
      colnames(t) <- rownames(t) <- names(coe.out)
      return(t)
    })
  }else if (method == "semi.ls") {
    names(coe.out) <- colnames(DF)[-c(1, 2, ncol(DF))]
    covmat <- lapply(covmat, function(t){
      colnames(t) <- rownames(t) <- names(coe.out)[-1]
      return(t)
    })
  }else if (method == "semi.rank.gehan.s") {
    # add intercept term for the rank based estimator
    names(coe.out) <- c("(Intercept)", colnames(DF)[-c(1, 2, ncol(DF))])
    covmat <- lapply(covmat, function(t){
      colnames(t) <- rownames(t) <- names(coe.out)[-1]
      return(t)
    })
  }
  out.time <- rbind(ptm2 - ptm1, ptm3 - ptm2, ptm4 - ptm3)
  rownames(out.time) <- c("ssps", "coefficients", "vcov")
  convergence <- c(optSSPs$converge, covg.out)
  names(convergence) <- c("pilot", paste0("sec_step", seq_along(covg.out)))
  iteration <- c(optSSPs$iteration, itr.out)
  names(iteration) <- c("pilot", paste0("sec_step", seq_along(itr.out)))
  out <- list(call = scall, vari.name = names(coe.out),
              coefficients = coe.out, covmat = covmat, 
              convergence = convergence,
              n.iteration = iteration,
              time = out.time,
              ssp.type = sspType, method = mtd.name,
              combine = combine, n.repeat = n.repeat)
  out$x <- DF[-c(1, 2, ncol(DF))]
  out$y <- DF[, c(1, 2)]
  out$ssp <- DF[, ncol(DF)]
  class(out) <- "aftosmac"
  return(out)
}

aftest <- function(formula, data, contrasts = NULL, subset,
                   method = c("parametric", "least_squares", "rank"),
                   control = list()) {
  method <- match.arg(method)
  scall <- match.call() # call function, including all inputs
  mnames <- c("", "formula", "data", "subset") # variable related to data
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
  DF <- as.data.frame(cbind(obj, model.matrix(mterms, m, contrasts))) # add covariates
  mtd.name <- method
  if (method == "least_squares") {
    method <- "semi.ls"
  }else if (method == "rank") {
    # delete intercept for rank-based approach
    DF <- DF[,-which(colnames(DF) == "(Intercept)")]
    method <- "semi.rank.gehan.s"
  }else if (method == "parametric") {
    method <- "par.weibull"
  }
  # create engine
  engine.control <- control[names(control) %in% names(attr(getClass(method), "slots"))]
  engine <- do.call("new", c(list(Class = method), engine.control))
  if (engine@b0 == 0) {
    if (method == "par.weibull") {
      engine@b0 <- c(1, as.numeric(rep(0, ncol(DF)-2)))
    }else if (method == "semi.rank.gehan.s") {
      engine@b0 <- as.numeric(lsfit(DF[, -(1:2)], log(DF[, 1]))$coefficient)[-1]
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
  
  # get full sample estimator
  ptm1 <- proc.time()
  
  engine@n <- nrow(DF)
  DF$ssps <- rep(1/engine@n, engine@n)
  est.full <- aftosmac.fit(DF, engine)
  
  coe.out <- est.full$coe
  covg.out <- est.full$converge
  
  ptm2 <- proc.time()
  if (covg.out != 0) {
    out <- list(call = scall, vari.name = colnames(DF)[-c(1, 2, ncol(DF))],
                coefficients = NA, covmat = NA, convergence = covg.out,
                n.iteration = NA, time = NA, method = mtd.name)
    out$x <- DF[-c(1, 2, ncol(DF))]
    out$y <- DF[, c(1, 2)]
    
    class(out) <- "aftest"
    return(out)
  }

  engine@b <- coe.out
  engine@ind_sub <- seq_len(engine@n)
  g <- aftosmac.est(DF, engine)
  vc <- crossprod(g) / engine@n
  m_inv <- solve(aftosmac.slope(DF, engine))
  covmat <- m_inv %*% vc %*% t(m_inv) / engine@n
  
  if (method == c("semi.rank.gehan.s")) {
    coe.out <- intcp(DF, engine)
  }
  
  if (method == "par.weibull") {
    names(coe.out) <- c("Scale", colnames(DF)[-c(1, 2, ncol(DF))])
    colnames(covmat) <- rownames(covmat) <- names(coe.out)
  }else if (method == "semi.ls") {
    names(coe.out) <- colnames(DF)[-c(1, 2, ncol(DF))]
    colnames(covmat) <- rownames(covmat) <- names(coe.out)[-1]
  }else if (method == "semi.rank.gehan.s") {
    # add intercept term for the rank based estimator
    names(coe.out) <- c("(Intercept)", colnames(DF)[-c(1, 2, ncol(DF))])
    colnames(covmat) <- rownames(covmat) <- names(coe.out)[-1]
  }
  ptm3 <- proc.time()
  
  out.time <- rbind(ptm2 - ptm1, ptm3 - ptm2)
  rownames(out.time) <- c("coefficients", "vcov")
  out <- list(call = scall, vari.name = names(coe.out),
              coefficients = coe.out, covmat = covmat, 
              convergence = covg.out,
              n.iteration = est.full$iter,
              time = out.time,
              method = mtd.name)
  out$x <- DF[-c(1, 2, ncol(DF))]
  out$y <- DF[, c(1, 2)]
  class(out) <- "aftest"
  return(out)
}


