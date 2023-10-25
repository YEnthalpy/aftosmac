## Prepare the dataset

data(lymp, package = "aftosmac")
x <- lymp$covariates
y <- lymp$response
delta <- lymp$delta

dat <- as.data.frame(cbind(x, y, delta))

## Fit the semi-prametric AFT model based on the rank-based approach
fit <- aftosmac(Surv(y, delta) ~ Age + Male + Nonwhite + Year,
                size.pilot = 500, size.subsample = 2000,
                sspType = "optA", model = "weibull", se = "parTrue",
                combine = "estimator", data = dat)
summary(fit)
