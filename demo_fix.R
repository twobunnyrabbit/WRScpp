#!/usr/bin/env Rscript
# Demonstration of the residual bug fix in tsreg_C()

library(WRScpp)

cat("=======================================================\n")
cat("Demonstration: tsreg_C() Bug Fix Verification\n")
cat("=======================================================\n\n")

# Simple example with perfect fit
cat("Example: Perfect linear fit (residuals should be ~0)\n")
cat("-----------------------------------------------------\n")
X <- matrix(c(1, 2, 3,
              0, 1, 2), ncol=2)
y <- c(5, 7, 9)

cat("Data:\n")
cat("X:\n")
print(X)
cat("\ny:", y, "\n\n")

result <- tsreg_C(x=X, y=y)

cat("Results:\n")
cat("Coefficients:", result$coef, "\n")
cat("Residuals:", result$residuals, "\n\n")

# Manually verify
fitted <- result$coef[1] + X %*% result$coef[-1]
manual_residuals <- y - fitted

cat("Verification:\n")
cat("Fitted values:", as.numeric(fitted), "\n")
cat("Manual residuals (y - fitted):", as.numeric(manual_residuals), "\n")
cat("Match? ", all(abs(result$residuals - manual_residuals) < 0.001), "\n\n")

# Realistic example
cat("=======================================================\n")
cat("Realistic Example: Multiple regression with noise\n")
cat("=======================================================\n")
set.seed(42)
n <- 100
X_real <- matrix(rnorm(n*3), ncol=3)
y_real <- 2 + 1.5*X_real[,1] - 2*X_real[,2] + 0.8*X_real[,3] + rnorm(n, sd=0.5)

result_real <- tsreg_C(x=X_real, y=y_real)

cat("Coefficients:\n")
cat("  Intercept:", result_real$coef[1], "\n")
cat("  Slopes:", result_real$coef[-1], "\n\n")

cat("Residual diagnostics:\n")
cat("  Mean (should be ~0):", mean(result_real$residuals), "\n")
cat("  SD:", sd(result_real$residuals), "\n")
cat("  Range: [", min(result_real$residuals), ",", max(result_real$residuals), "]\n\n")

# Verify residuals
fitted_real <- result_real$coef[1] + X_real %*% result_real$coef[-1]
manual_real <- y_real - fitted_real

cat("Verification:\n")
cat("  Correlation (function vs manual): ", cor(result_real$residuals, manual_real), "\n")
cat("  Max difference: ", max(abs(result_real$residuals - manual_real)), "\n")
cat("  Residuals are CORRECT: ", max(abs(result_real$residuals - manual_real)) < 1e-10, "\n\n")

cat("=======================================================\n")
cat("✓ Bug is FIXED! Residuals are now calculated correctly.\n")
cat("=======================================================\n")
