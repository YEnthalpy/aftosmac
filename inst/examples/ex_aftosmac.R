## Prepare the dataset

data(lymp, package = "aftosmac")

# center and scale continuous covariates
new_lymp <- lymp
new_lymp$Age <- scale(new_lymp$Age)
new_lymp$Diagnosed_year <- scale(new_lymp$Diagnostic_year)

## Fit the semi-prametric AFT model based on the rank-based approach
fit <- aftosmac(Surv(Survtime, Status) ~ Age + Male + Nonwhite + Diagnostic_year,
                size.pilot = 500, size.subsample = 2000,
                sspType = "optA", model = "weibull", se = "parTrue",
                combine = "estimator", data = new_lymp)
summary(fit)
