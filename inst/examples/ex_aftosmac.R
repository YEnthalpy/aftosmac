## Prepare the dataset
data(lymp, package = "aftosmac")

## center and scale continuous covariates
new_lymp <- lymp
new_lymp[, c("Age", "Diagnostic_year")] <- 
  lapply(lymp[, c("Age", "Diagnostic_year")], scale)

## Fit the semi-prametric AFT model based on the rank-based approach
set.seed(100)
fit <- aftosmac(Surv(Survtime, Status) ~ Age + Male + 
                  Nonwhite + Diagnostic_year,
                n.pilot = 500, n.sub = 4000,
                sspType = "optA", method = "parametric",
                combine = "estimator", data = new_lymp)
summary(fit)
