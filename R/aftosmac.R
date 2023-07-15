"%^%" <- function(x, n) {
  with(eigen(x), vectors %*% (values^n * t(vectors)))
}
# least-square approach
## 1. KM estimate and the integration part in yhat
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


## 2. Estimating function for least square approach
uls <- function(x, y, delta, beta, pi, n, ind_km) {
  xdif <- center(x[, -1], pi, n)
  py <- x %*% beta
  tmp <- eres((y - py), delta, pi, ind_km)[[1]]
  hy <- delta * y + (1 - delta) * (tmp + py)
  ydif <- hy - mean(hy / pi) / n
  beta <- beta[-1]
  out <- xdif * as.vector(ydif - xdif %*% beta) / pi / n
  return(out)
}

## 3. resampling approach to get slope matrix
resp <- function(x, y, delta, beta, pi, n, b) {
  p <- ncol(x) - 1
  r <- nrow(x)
  z <- matrix(rexp(b * r), nrow = r, byrow = FALSE)
  g <- uls(x, y, delta, beta, pi, n, seq_along(y))
  v <- var(crossprod(z, g)) / r
  zbs <- matrix(ifelse(rbinom(b * p, 1, 0.5) == 1, 1, -1),
    ncol = p, byrow = TRUE
  )
  zb <- zbs %*% (v %^% (-0.5))
  beta1 <- zb / sqrt(r) + matrix(rep(beta[-1], b),
    ncol = p,
    byrow = TRUE
  )
  beta1 <- cbind(rep(beta[1], b), beta1)
  response <- matrix(NA, nrow = b, ncol = p)
  for (i in seq_len(b)) {
    response[i, ] <- colSums(uls(
      x, y, delta, beta1[i, ],
      pi, n, seq_along(y)
    )) / sqrt(r)
  }
  return(t(lsfit(zb, response, intercept = FALSE)$coefficients))
}


## 4. Iterative procedure to get the estimator
semi_ls_est <- function(x, y, delta, pi, n, control = list(
                          init = rep(0, ncol(x)),
                          tol = 1e-5, maxit = 200
                        )) {
  beta <- control$init
  xdif <- center(x[, -1], pi, n)
  for (i in 1:control$maxit) {
    py <- x %*% beta
    tmp <- eres((y - py), delta, pi, seq_along(y))[[1]]
    hy <- delta * y + (1 - delta) * (tmp + py)
    ydif <- as.vector(hy - mean(hy / pi) / n)
    beta_new <- tryCatch(
      solve(
        crossprod(xdif, xdif / pi),
        colSums(xdif * ydif / pi)
      ),
      error = function(e) NA,
      warning = function(w) NA
    )
    if (is.na(beta_new[1])) {
      return(list(coefficient = rep(0, ncol(x)), converge = 1, ite = NA))
    }
    beta_new <- c(NA, as.vector(beta_new))
    beta_new[1] <- max(eres(
      (y - x[, -1] %*% beta_new[-1]),
      delta, pi, seq_along(y)
    )[[2]])
    e <- sqrt(sum((beta_new - beta)^2))
    if (e < control$tol) {
      return(list(coefficient = beta_new, converge = 0, ite = i))
    } else {
      beta <- beta_new
    }
    if (i == control$maxit) {
      return(list(coefficient = beta_new, converge = 2, ite = i))
    }
  }
}

## 5. get the optimal ssps
semi_ls_ssp <- function(x, y, delta, r0, ssp_type, b, alpha) {
  n <- nrow(x)
  if (ssp_type == "uniform") {
    ssp <- rep(1 / n, n)
    return(list(ssp = ssp, index.pilot = sample(n, r0, TRUE), converge = 0))
  } else {
    ssp_pt <- rep(1 / n, r0)
    ind_pt <- sample(n, r0, TRUE)
    x_pt <- x[ind_pt, ]
    y_pt <- y[ind_pt]
    dt_pt <- delta[ind_pt]
    mle_pt <- semi_ls_est(x_pt, y_pt, dt_pt, ssp_pt, n)
    if (mle_pt$converge != 0) {
      return(list(bt_pt = NA, ssp = NA, index.pilot = NA,
                  converge = mle_pt$converge))
    }
    bt_pt <- mle_pt$coefficient
    g <- uls(x, y, delta, bt_pt, rep(1 / n, n), n, ind_pt)
    if (ssp_type == "optL") {
      g_nm <- sqrt(rowSums(g^2))
      ssp <- g_nm / sum(g_nm) * (1 - alpha) + alpha / n
    } else if (ssp_type == "optA") {
      m_inv <- solve(resp(x_pt, y_pt, dt_pt, bt_pt, ssp_pt, n, b))
      m_mse <- sqrt(colSums((tcrossprod(m_inv, g))^2))
      ssp <- m_mse / sum(m_mse) * (1 - alpha) + alpha / n
    }
    return(list(ssp = ssp, index.pilot = ind_pt, converge = 0))
  }
}


## 6. Get estiamted coefficient and standard error
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

## 2. Estimating function for least square approach
uls <- function(x, y, delta, beta, pi, n, ind_km) {
  xdif <- center(x[, -1], pi, n)
  py <- x %*% beta
  tmp <- eres((y - py), delta, pi, ind_km)[[1]]
  hy <- delta * y + (1 - delta) * (tmp + py)
  ydif <- hy - mean(hy / pi) / n
  beta <- beta[-1]
  out <- xdif * as.vector(ydif - xdif %*% beta) / pi / n
  return(out)
}

## 3. resampling approach to get slope matrix
resp <- function(x, y, delta, beta, pi, n, b) {
  p <- ncol(x) - 1
  r <- nrow(x)
  z <- matrix(rexp(b * r), nrow = r, byrow = FALSE)
  g <- uls(x, y, delta, beta, pi, n, seq_along(y))
  v <- var(crossprod(z, g)) / r
  zbs <- matrix(ifelse(rbinom(b * p, 1, 0.5) == 1, 1, -1),
                ncol = p, byrow = TRUE
  )
  zb <- zbs %*% (v %^% (-0.5))
  beta1 <- zb / sqrt(r) + matrix(rep(beta[-1], b),
                                 ncol = p,
                                 byrow = TRUE
  )
  beta1 <- cbind(rep(beta[1], b), beta1)
  response <- matrix(NA, nrow = b, ncol = p)
  for (i in seq_len(b)) {
    response[i, ] <- colSums(uls(
      x, y, delta, beta1[i, ],
      pi, n, seq_along(y)
    )) / sqrt(r)
  }
  return(t(lsfit(zb, response, intercept = FALSE)$coefficients))
}


## 4. Iterative procedure to get the estimator
semi_ls_est <- function(x, y, delta, pi, n, control = list(
  init = rep(0, ncol(x)),
  tol = 1e-5, maxit = 200
)) {
  beta <- control$init
  xdif <- center(x[, -1], pi, n)
  for (i in 1:control$maxit) {
    py <- x %*% beta
    tmp <- eres((y - py), delta, pi, seq_along(y))[[1]]
    hy <- delta * y + (1 - delta) * (tmp + py)
    ydif <- as.vector(hy - mean(hy / pi) / n)
    beta_new <- tryCatch(
      solve(
        crossprod(xdif, xdif / pi),
        colSums(xdif * ydif / pi)
      ),
      error = function(e) NA,
      warning = function(w) NA
    )
    if (is.na(beta_new[1])) {
      return(list(coefficient = rep(0, ncol(x)), converge = 1, ite = NA))
    }
    beta_new <- c(NA, as.vector(beta_new))
    beta_new[1] <- max(eres(
      (y - x[, -1] %*% beta_new[-1]),
      delta, pi, seq_along(y)
    )[[2]])
    e <- sqrt(sum((beta_new - beta)^2))
    if (e < control$tol) {
      return(list(coefficient = beta_new, converge = 0, ite = i))
    } else {
      beta <- beta_new
    }
    if (i == control$maxit) {
      return(list(coefficient = beta_new, converge = 2, ite = i))
    }
  }
}

## 5. get the optimal ssps
semi_ls_ssp <- function(x, y, delta, r0, ssp_type, b, alpha) {
  n <- nrow(x)
  if (ssp_type == "uniform") {
    ssp <- rep(1 / n, n)
    return(list(ssp = ssp, index.pilot = sample(n, r0, TRUE), converge = 0))
  } else {
    ssp_pt <- rep(1 / n, r0)
    ind_pt <- sample(n, r0, TRUE)
    x_pt <- x[ind_pt, ]
    y_pt <- y[ind_pt]
    dt_pt <- delta[ind_pt]
    mle_pt <- semi_ls_est(x_pt, y_pt, dt_pt, ssp_pt, n)
    if (mle_pt$converge != 0) {
      return(list(bt_pt = NA, ssp = NA, index.pilot = NA,
                  converge = mle_pt$converge))
    }
    bt_pt <- mle_pt$coefficient
    g <- uls(x, y, delta, bt_pt, rep(1 / n, n), n, ind_pt)
    if (ssp_type == "optL") {
      g_nm <- sqrt(rowSums(g^2))
      ssp <- g_nm / sum(g_nm) * (1 - alpha) + alpha / n
    } else if (ssp_type == "optA") {
      m_inv <- solve(resp(x_pt, y_pt, dt_pt, bt_pt, ssp_pt, n, b))
      m_mse <- sqrt(colSums((tcrossprod(m_inv, g))^2))
      ssp <- m_mse / sum(m_mse) * (1 - alpha) + alpha / n
    }
    return(list(ssp = ssp, index.pilot = ind_pt, converge = 0))
  }
}


## 6. Get estiamted coefficient and standard error
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
