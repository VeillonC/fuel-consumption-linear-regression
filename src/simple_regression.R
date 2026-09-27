
# Simple linear regression, from scratch (normal equations), validated
# against R's lm(). Fits fuel consumption (1/mpg, in L/100km) as a linear
# function of engine displacement, on a random subset of the mtcars dataset.

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

# Simple linear regression function (least squares, closed form)
linear_regression <- function(X, Y) {
  n <- length(Y)
  sum_x <- sum(X)
  sum_y <- sum(Y)
  sum_xy <- sum(X * Y)
  sum_x2 <- sum(X^2)

  a <- (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x^2)
  b <- (sum_y - a * sum_x) / n

  list(a = a, b = b)
}

# R-squared
r_squared <- function(X, Y, a, b) {
  y_pred <- a * X + b
  ss_tot <- sum((Y - mean(Y))^2)
  ss_res <- sum((y_pred - mean(Y))^2)
  ss_res / ss_tot
}

X_explanatory <- X[, 7]  # explanatory variable: qsec (quarter mile time, in seconds)
X_s <- sort(X_explanatory)
Y_s <- sort(Y)  # response variable: mpg

result <- linear_regression(X_s, Y_s)
a <- result$a
b <- result$b
n <- length(Y_s)
alpha <- 0.05

plot(X_s, Y_s, main = "Simple linear regression: mpg vs qsec", ylab = "Fuel efficiency (mpg)", xlab = "Quarter mile time, qsec (s)", col = "blue", pch = 10)
abline(a = result$b, b = result$a, col = "red", lwd = 1)

r2 <- r_squared(X_s, Y_s, a = result$a, b = result$b)
if (r2 > 0.5) {
  cat("R-squared =", r2, "> 0.5, so the linear regression is valid.\n")
} else {
  cat("R-squared is too low, the linear regression is not a good fit.\n")
}

# Predictions and residuals
y_hat <- a * X_s + b
residuals <- Y_s - y_hat

# Residual standard error
resid_se <- sqrt(sum(residuals^2) / (n - 2))

# Leverage (hat values)
X_design <- cbind(1, X_s)
leverage <- diag(X_design %*% solve(t(X_design) %*% X_design) %*% t(X_design))
# Standardized residuals
std_residuals <- residuals / (resid_se * sqrt(1 - leverage))
unique_std_residuals <- unique(std_residuals)
# Kolmogorov-Smirnov test
ks <- ks.test(unique_std_residuals, "pnorm")
ks_pvalue <- ks$p.value
if (ks_pvalue > alpha) {
  cat("KS p-value =", ks_pvalue, "> 5%, the residuals are consistent with a Gaussian distribution.\n")
} else {
  cat("The residuals are not Gaussian.\n")
}

t_a <- a / resid_se / sqrt(sum(X_s - mean(X_s))^2)
p_value <- 1 - pt(abs(t_a), df = n - 2)
if (p_value < alpha) {
  cat("p-value for the slope test =", p_value, "< 5%, we reject H0: a = 0.\n")
} else {
  cat("We do not reject H0: no significant evidence that a differs from 0.\n")
}

# Prediction interval for a new x value (here, the sample mean)
x_new <- mean(X_s)
x_bar <- mean(X_s)
s_xx <- sum((X_s - x_bar)^2)
t_crit <- qt(1 - alpha / 2, df = n - 2)

y_pred <- a * x_new + b
margin <- t_crit * resid_se * sqrt(1 + (1 / n) + ((x_new - x_bar)^2 / s_xx))
ci_lower <- y_pred - margin
ci_upper <- y_pred + margin

cat("95% confidence interval for x_new =", x_new, "\n")
cat("[", ci_lower, ",", ci_upper, "]\n")

plot(Y_s, a * X_s + b, main = "Observed vs predicted Y", xlab = "Observed Y", ylab = "Predicted Y", pch = 10)
abline(0, 1, col = "red", lwd = 2)                 # Y = X reference line
abline(-margin, 1, col = "blue", lty = 2)          # lower confidence bound
abline(margin, 1, col = "darkgreen", lty = 2)      # upper confidence bound
legend("bottomright",
       legend = c("Predicted Y", "Y = X", "CI lower", "CI upper"),
       col = c("black", "red", "blue", "darkgreen"),
       pch = c(10, NA, NA, NA),
       lty = c(NA, 1, 2, 2),
       lwd = c(NA, 2, 1, 1))

qqnorm(unique_std_residuals,
       main = "QQ-plot of the standardized residuals",
       xlab = "Theoretical quantiles",
       ylab = "Empirical quantiles (standardized)",
       pch = 19,
       col = "springgreen")
qqline(unique_std_residuals, lwd = 3, col = "tomato")


### Cross-check against R's built-in lm() ###
cat("Cross-check with lm():\n")
lm_model <- lm(Y_s ~ X_s)
studentized_residuals <- unique(rstudent(lm_model))
lm_ks_test <- ks.test(studentized_residuals, "pt", df = lm_model$df.residual)
print(lm_ks_test)
cat("R-squared =", summary(lm_model)$r.squared, "\n")
cat("p-value for H0 (slope = 0):", summary(lm_model)$coefficients[2, 4], "\n", "We conclude H1: a != 0")
qqnorm(studentized_residuals,
       main = "QQ-plot of the studentized residuals (via lm)",
       xlab = "Theoretical quantiles",
       ylab = "Empirical quantiles (studentized)",
       pch = 19,
       col = "springgreen")
qqline(studentized_residuals, lwd = 3, col = "tomato")

# Empirical CDF of the studentized residuals vs the theoretical Student's t CDF
empirical_cdf <- ecdf(studentized_residuals)
sorted_vals <- sort(studentized_residuals)
x_vals <- seq(min(sorted_vals), max(sorted_vals), 0.01)
cumulative_prob <- empirical_cdf(x_vals)
student_theoretical <- pt(x_vals, df = n - 3)
plot(x_vals, cumulative_prob, type = "s", col = "red", main = "Empirical vs theoretical CDF of the residuals", xlab = "x", ylab = "Cumulative probability")
lines(x_vals, student_theoretical, type = "l", col = "blue")
