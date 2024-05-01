## Point estimator
setClass("Engine",
         representation(tol = "numeric", b0 = "numeric", b = "numeric",
                        maxit = "numeric", n = "numeric", ind_sub = "numeric",
                        B = "numeric", r0 = "numeric", r = "numeric",
                        alpha = "numeric", tol_pilot = "numeric"),
         prototype(tol = 1e-5, b0 = 0, b = 0, maxit = 200,
                   n = 0, ind_sub = 0, B = 20, r0 = 0, r = 0, alpha = 0.2,
                   tol_pilot = 1e-5),
         contains = "VIRTUAL")

setClass("semi.ls", contains = "Engine")
setClass("semi.rank.gehan.s", contains = "Engine")
setClass("semi.rank.gehan.ns", contains = "Engine")
setClass("par.weibull", contains = "Engine")
