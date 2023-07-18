## Point estimator
setClass("Engine",
         representation(tol = "numeric", b0 = "numeric", b = "numeric",
                        maxit = "numeric", n = "numeric", ind_sub = "numberic",
                        B = "numeric", r0 = "numeric", alpha = "numeric"),
         prototype(tol = 1e-5, b0 = 0, b = 0, maxit = 1000,
                   n = 0, ind_sub = 0, B = 20, r0 = 0, alpha = 0.2),
         contains = "VIRTUAL")

setClass("semi.ls", contains = "Engine")
setClass("semi.rank.gehan.s", contains = "Engine")
setClass("semi.rank.gehan.ns", contains = "Engine")
setClass("par.weibull", contains = "Engine")

## Generic function -- subsample estimates
setGeneric("aftosmac.fit", function(DF, engine) standardGeneric("aftosmac.fit"))

setMethod("aftosmac.fit", signature(engine = "semi.ls"), lsFit)
setMethod("aftosmac.fit", signature(engine = "semi.rank.gehan.s"), rankFit.gehan.s)
setMethod("aftosmac.fit", signature(engine = "semi.rank.gehan.ns"), rankFit.gehan.ns)
setMethod("aftosmac.fit", signature(engine = "par.weibull"), parFit.weibull)

## Generic function -- estimating functions
setGeneric("aftosmac.est", function(DF, engine) standardGeneric("aftosmac.est"))

setMethod("aftosmac.est", signature(engine = "semi.ls"), lsEst)
setMethod("aftosmac.est", signature(engine = "semi.rank.gehan.s"), rankEst.gehan.s)
setMethod("aftosmac.est", signature(engine = "semi.rank.gehan.ns"), rankEst.gehan.ns)
setMethod("aftosmac.est", signature(engine = "par.weibull"), parEst.weibull)

## Generic function -- slopes of estimating functions
setGeneric("aftosmac.slope", function(DF, engine) standardGeneric("aftosmac.slope"))

setMethod("aftosmac.slope", signature(engine = "semi.ls"), lsSlp)
setMethod("aftosmac.slope", signature(engine = "semi.rank.gehan.s"), rankSlp.gehan.s)
setMethod("aftosmac.slope", signature(engine = "semi.rank.gehan.ns"), rankSlp.gehan.ns)
setMethod("aftosmac.slope", signature(engine = "par.weibull"), parSlp.weibull)


