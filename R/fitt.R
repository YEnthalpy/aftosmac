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
    # delete intercept in rank-based approach
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
  engine@n <- nrow(DF)
  optSSPs <- aftosmac.ssps(DF, engine, fitMtd = fitMtd,
                           rankWt = rankWt, sspType = sspType)

}

semi_ls_fit <- function(x, y, delta, r0, r, ssp_type, method, se = TRUE, b = 20,
                        itr = 5, alpha = 0.2) {
  n <- nrow(x)
  ssps <- semi_ls_ssp(x, y, delta, r0, ssp_type, b, alpha)
  pi <- ssps$ssp
  if (ssps$converge != 0) {
    stop(paste0("Fail to get a converging pilot estimator. The converging code is ", ssps$converge))
  }
  if (method == "one") {
    ind_pt <- ssps$index.pilot
    ind_r <- sample(n, r, prob = pi, replace = TRUE)
    ind <- c(ind_r, ind_pt)
    sec_ssp <- c(pi[ind_r], rep(1 / n, r0))
    est <- semi_ls_est(x[ind, ], y[ind], delta[ind], sec_ssp, n)
    if (est$converge != 0) {
      stop(paste0("Fail to get a converging second-step estimator. The converging code is ", est$converge))
    }
    coe <- as.vector(est$coefficient)
    ite <- est$ite
    if (se) {
      r <- r + r0
      z <- matrix(rexp(b * r), nrow = r, byrow = FALSE)
      g <- uls(x[ind, ], y[ind], delta[ind], coe, sec_ssp, n, seq_len(r))
      vc <- var(crossprod(z, g)) / r
      m_inv <- solve(resp(x[ind, ], y[ind], delta[ind], coe, sec_ssp, n, b))
      vx <- m_inv %*% tcrossprod(vc, m_inv) / r
      std <- c(NA, sqrt(diag(vx)))
    } else {
      std <- NA
    }
  } else if (method == "multiple") {
    out <- replicate(itr, expr = {
      ind <- sample(n, r, prob = pi, replace = TRUE)
      est <- list(semi_ls_est(x[ind, ], y[ind], delta[ind], pi[ind], n))
    })
    conv <- sapply(out, function(i) {i$converge})
    if (sum(conv) != 0) {
      stop(paste0("At least one second-step estimators do not converge. ",
                  length(which(conv == 1)), "subsample estimators have the converging code 1. ",
                  length(which(conv == 2)), "subsample estimators have the converging code 2."))
    }
    coe_itr <- sapply(out, function(t){t$coefficient})
    ite <- mean(sapply(out, function(t){t$ite}))
    coe <- rowMeans(coe_itr)
    if (se) {
      std <- sqrt(diag(var(t(coe_itr)) / itr))
    } else {
      std <- NA
    }
  }
  if (is.null(colnames(x))) {
    names(coe) <- names(std) <- c("Intercept", paste0("Beta", seq(1, ncol(x)-1, 1)))
  }else {
    names(coe) <- names(std) <- colnames(x)
  }
  return(list(coefficient = coe, std = std, converge = 0, ite = ite))
}
