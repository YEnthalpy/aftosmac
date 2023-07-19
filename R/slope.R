lsSlp <- function(DF, engine) {
  estSlp(DF, engine, fitMtd = "ls")
}

rankSlp.gehan.ns <- function(DF, engine) {
  estSlp(DF, engine, fitMtd = "rank", rankWt = "gehan")
}

rankSlp.gehan.s <- function(DF, engine) {
  xmat <- as.matrix(DF$covaraites)
  y <- log(DF$time)
  gehan_s_jaco(xmat, y, DF$status, DF$ssps, engine@b, engine@n)
}

parSlp.weibull <- function(DF, engine) {
  xmat <- as.matrix(DF$covaraites)
  y <- log(DF$time)
  sigma <- engine@b[1]
  er <- drop((y - xmat %*% engine@b[-1]) / sigma)
  exper <- exp(er)
  d2.beta <- -t(xmat * exper  / DF$ssps) %*% xmat / (sigma^2)
  d2.sig <- sum((DF$status + 2 * er * (DF$status - exper) - (er^2)
                 * exper) / DF$ssps) / (sigma^2)
  d2.sigbeta <- colSums((xmat * (DF$status - exper * (1 + er)))) / (sigma ^ 2)
  cbind(c(d2.sig, d2.sigbeta), rbind(t(d2.sigbeta), d2.beta))
}

estSlp <- function(DF, engine, fitMtd = c("rank", "ls"),
                   rankWt = c("gehan")) {
  xmat <- as.matrix(DF$covaraites)
  y <- log(DF$time)
  ssps <- DF$ssps
  p <- ncol(xmat) - 1
  r <- length(y)
  fitMtd <- match.arg(fitMtd)
  rankWt <- match.arg(rankWt)
  if (fitMtd == "ls") {
    method <- "semi.ls"
  }else if (fitMtd == "rank") {
    method <- paste("semi", fitMtd, rankWt, "s", sep = ".")
  }
  engine <- do.call("new", list(Class = method))
  z <- matrix(rexp(engine@B * r), nrow = r, byrow = FALSE)
  g <- aftosmac.est(DF = DF, engine = engine)
  v <- var(crossprod(z, g)) / r
  zbs <- matrix(ifelse(rbinom(engine@B * p, 1, 0.5) == 1, 1, -1),
                ncol = p, byrow = TRUE)
  zb <- zbs %*% with(eigen(v), vectors %*% (values^(-0.5) * t(vectors)))
  # need to define b, B and ind_sub
  beta1 <- zb / sqrt(r) + matrix(rep(engine@b, engine@B), ncol = p, byrow = TRUE)
  beta1 <- cbind(rep(beta[1], engine@B), beta1)
  response <- matrix(NA, nrow = engine@B, ncol = p)
  for (i in seq_len(engine@B)) {
    engine@b <- beta1[i, ]
    response[i, ] <- colSums(aftosmac.est(DF = DF, engine = engine)) / sqrt(r)
  }
  t(lsfit(zb, response, intercept = FALSE)$coefficients)
}

## Generic function -- slopes of estimating functions
setGeneric("aftosmac.slope", function(DF, engine) standardGeneric("aftosmac.slope"))

setMethod("aftosmac.slope", signature(engine = "semi.ls"), lsSlp)
setMethod("aftosmac.slope", signature(engine = "semi.rank.gehan.s"), rankSlp.gehan.s)
setMethod("aftosmac.slope", signature(engine = "semi.rank.gehan.ns"), rankSlp.gehan.ns)
setMethod("aftosmac.slope", signature(engine = "par.weibull"), parSlp.weibull)
