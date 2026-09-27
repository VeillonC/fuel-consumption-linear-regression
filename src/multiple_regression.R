library(plotly)

# Multiple linear regression, from scratch (normal equations), validated
# against R's lm(). Fits mpg as a linear function of displacement and
# weight, with a 3D visualization and the usual validity checks (R²,
# F-test, normality of residuals, confidence intervals on the coefficients).

# Build the working dataset
A <- mtcars
F1 <- as.factor(A[, 8])
A[, 8] <- F1
set.seed(4)
set.seed(4 * floor(100 * runif(1, 0, 3)))
row_to_drop <- sample(1:32, 1)
B <- A[-row_to_drop, ]
Y <- B[, 1]
all_cols <- 1:11
candidate_cols <- all_cols[-c(1, 8, 9)]
selected_cols <- c(8, sample(candidate_cols, 6, replace = FALSE))
X <- B[, selected_cols]

X <- B[, c("disp", "wt")]
X_reg <- cbind(1, X$disp, X$wt)  # 1 for the intercept

# Manual linear regression (normal equations)
regression <- function(X, Y) {
  Xt <- t(X)
  beta <- solve(Xt %*% X) %*% Xt %*% Y
  list(coefficients = as.vector(beta))
}

result <- regression(X_reg, Y)
coef <- result$coefficients

# Cross-check with lm()
lm_model <- lm(mpg ~ disp + wt, data = B)

# 3D visualization
disp_seq <- seq(min(B$disp), max(B$disp), length.out = 20)
wt_seq <- seq(min(B$wt), max(B$wt), length.out = 20)

mpg_pred <- outer(disp_seq, wt_seq, function(d, w) coef[1] + coef[2] * d + coef[3] * w)

plot_ly(B) %>%
  add_markers(x = ~disp, y = ~wt, z = ~mpg, marker = list(size = 5, color = 'black')) %>%
  add_surface(x = disp_seq, y = wt_seq, z = t(mpg_pred), colorscale = "Rainbow", opacity = 0.3) %>%
  layout(scene = list(
    xaxis = list(title = "disp"),
    yaxis = list(title = "wt"),
    zaxis = list(title = "mpg")
  ))


############################## VALIDITY CHECKS ##############################

# Predictions
y_hat <- X_reg %*% coef
residuals <- Y - y_hat

# Manual R-squared
tss <- sum((Y - mean(Y))^2)
rss <- sum(residuals^2)
ess <- sum((y_hat - mean(Y))^2)
r2 <- ess / tss

n <- nrow(X_reg)
p <- ncol(X_reg) - 1  # without the intercept
r2_adj <- 1 - (1 - r2) * (n - 1) / (n - p - 1)

cat(sprintf("R-squared (manual) = %.4f\n", r2))
cat(sprintf("Adjusted R-squared (manual) = %.4f\n", r2_adj))

# Global F-test
f_stat <- ((tss - rss) / p) / (rss / (n - p - 1))
p_value_f <- 1 - pf(f_stat, p, n - p - 1)
cat("\nF-test: F =", f_stat, "with p-value =", p_value_f, "\n")
if (p_value_f < 0.05) {
  cat("=> Significant model (H0 rejected)\n")
} else {
  cat("=> Model not significant\n")
}

# Normality of residuals: KS test
leverage <- diag(X_reg %*% solve(t(X_reg) %*% X_reg) %*% t(X_reg))
resid_se <- sqrt(sum(residuals^2) / (n - p - 1))
studentized_residuals <- residuals / (resid_se * sqrt(1 - leverage))

ks_test <- ks.test(studentized_residuals, "pnorm")
cat("\nKS normality test:\n")
cat("KS p-value =", ks_test$p.value, "\n")
if (ks_test$p.value > 0.05) {
  cat("=> Residuals are Gaussian (normality holds)\n")
} else {
  cat("=> Residuals are not Gaussian\n")
}


# Comparison with lm()'s own values
summary(lm_model)$r.squared        # R-squared from lm
summary(lm_model)$adj.r.squared    # adjusted R-squared from lm


# 95% confidence intervals
n <- nrow(X_reg)
p <- ncol(X_reg) - 1
sigma2_hat <- sum(residuals^2) / (n - p - 1)

# Coefficient variance-covariance matrix
xtx_inv <- solve(t(X_reg) %*% X_reg)
beta_var <- sigma2_hat * xtx_inv

# Standard errors of the coefficients
beta_se <- sqrt(diag(beta_var))

alpha <- 0.05
t_crit <- qt(1 - alpha / 2, df = n - p - 1)

lower_bounds <- coef - t_crit * beta_se
upper_bounds <- coef + t_crit * beta_se

cat("\n95% confidence interval for each coefficient:\n")
for (i in 1:length(coef)) {
  cat(sprintf("B%d : [%.4f ; %.4f]\n", i - 1, lower_bounds[i], upper_bounds[i]))
}
