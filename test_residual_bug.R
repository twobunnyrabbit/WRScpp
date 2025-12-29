#!/usr/bin/env Rscript
# Test script to verify the potential bug in tsreg_C residual calculation
# Bug location: src/robustmethods_CPP.cpp line 546
# Current code: res = y = x*temp - coef(0)
# Expected code: res = y - x*temp - coef(0)

library(WRScpp)

cat("=======================================================\n")
cat("Testing for Residual Bug in tsreg_C()\n")
cat("=======================================================\n\n")

# Create a simple test case with known structure
set.seed(123)
n <- 50
p <- 3

# Generate clean data with known relationship
X <- matrix(rnorm(n*p), ncol=p)
true_intercept <- 2.5
true_slopes <- c(1.5, -2.0, 0.8)
true_coef <- c(true_intercept, true_slopes)

# Generate y with small noise
y_original <- true_intercept + X %*% true_slopes + rnorm(n, sd=0.1)

# Make a copy to preserve original
y_test <- as.numeric(y_original)
X_test <- X

cat("TEST 1: Single Predictor (Should work correctly)\n")
cat("================================================\n")
x1 <- X[,1]
y1 <- as.numeric(y_original)

result1 <- tsreg_C(x=x1, y=y1)

# Manually compute residuals
fitted1 <- result1$coef[1] + result1$coef[2] * x1
manual_residuals1 <- y1 - fitted1

cat("Coefficients:", result1$coef, "\n")
cat("\nResidual diagnostics:\n")
cat("  Function residuals - mean:", mean(result1$residuals), "sd:", sd(result1$residuals), "\n")
cat("  Manual residuals   - mean:", mean(manual_residuals1), "sd:", sd(manual_residuals1), "\n")
cat("  Difference (should be ~0):", mean(abs(result1$residuals - manual_residuals1)), "\n")

if (mean(abs(result1$residuals - manual_residuals1)) < 0.001) {
  cat("✓ Single predictor residuals are CORRECT\n\n")
} else {
  cat("✗ Single predictor residuals are INCORRECT\n\n")
}

cat("\n=======================================================\n")
cat("TEST 2: Multiple Predictors (Bug location)\n")
cat("=======================================================\n")

# Store original values
y_before <- y_test
X_before <- X_test

result2 <- tsreg_C(x=X_test, y=y_test)

cat("Coefficients:", result2$coef, "\n\n")

# Manually compute what residuals SHOULD be
# fitted = intercept + X %*% slopes
fitted2_manual <- result2$coef[1] + X_before %*% result2$coef[-1]
manual_residuals2 <- y_before - fitted2_manual

cat("Residual Analysis:\n")
cat("------------------\n")
cat("Returned residuals:\n")
cat("  Mean:  ", mean(result2$residuals), "\n")
cat("  SD:    ", sd(result2$residuals), "\n")
cat("  Range: [", min(result2$residuals), ",", max(result2$residuals), "]\n")
cat("  Sum:   ", sum(result2$residuals), "\n\n")

cat("Manual residuals (y - fitted):\n")
cat("  Mean:  ", mean(manual_residuals2), "\n")
cat("  SD:    ", sd(manual_residuals2), "\n")
cat("  Range: [", min(manual_residuals2), ",", max(manual_residuals2), "]\n")
cat("  Sum:   ", sum(manual_residuals2), "\n\n")

cat("Difference between returned and manual residuals:\n")
cat("  Max absolute diff:", max(abs(result2$residuals - manual_residuals2)), "\n")
cat("  Mean absolute diff:", mean(abs(result2$residuals - manual_residuals2)), "\n")
cat("  Correlation:", cor(result2$residuals, manual_residuals2), "\n\n")

# Check if the bug manifests as the chained assignment
# If bug exists: res = y = x*temp - coef(0) would set res to fitted values (wrong sign)
# Expected: res = y - x*temp - coef(0) (correct residuals)

# The bug would make residuals equal to -(fitted - mean(y)) approximately
# Let's check various hypotheses

cat("Hypothesis Testing:\n")
cat("-------------------\n")

# Hypothesis 1: Residuals are correct (y - fitted)
hyp1_match <- mean(abs(result2$residuals - manual_residuals2)) < 0.001
cat("H1: Residuals = y - fitted values:          ", hyp1_match, "\n")

# Hypothesis 2: Residuals are fitted values (completely wrong)
hyp2_residuals <- as.numeric(fitted2_manual)
hyp2_match <- mean(abs(result2$residuals - hyp2_residuals)) < 0.001
cat("H2: Residuals = fitted values:              ", hyp2_match, "\n")

# Hypothesis 3: Residuals are negative fitted (from bug)
hyp3_residuals <- -as.numeric(fitted2_manual)
hyp3_match <- mean(abs(result2$residuals - hyp3_residuals)) < 0.001
cat("H3: Residuals = -fitted values:             ", hyp3_match, "\n")

# Hypothesis 4: Some variation of the bug
# res = y = x*temp - coef(0) means res gets x*temp - coef(0)
# This would be fitted_without_intercept - intercept
hyp4_residuals <- as.numeric(X_before %*% result2$coef[-1] - result2$coef[1])
hyp4_match <- mean(abs(result2$residuals - hyp4_residuals)) < 0.001
cat("H4: Residuals = X*slopes - intercept:       ", hyp4_match, "\n")

# Hypothesis 5: Correct but with sign flipped
hyp5_residuals <- -manual_residuals2
hyp5_match <- mean(abs(result2$residuals - hyp5_residuals)) < 0.001
cat("H5: Residuals = -(y - fitted):              ", hyp5_match, "\n\n")

cat("=======================================================\n")
cat("CONCLUSION:\n")
cat("=======================================================\n")

if (hyp1_match) {
  cat("✓ Residuals are CORRECT despite the suspicious code!\n")
  cat("  The line 'res = y = x*temp - coef(0)' might be correct\n")
  cat("  in the context of Armadillo matrix operations.\n")
} else if (hyp2_match || hyp3_match || hyp4_match || hyp5_match) {
  cat("✗ BUG CONFIRMED! Residuals are computed INCORRECTLY.\n")
  if (hyp4_match) {
    cat("  The bug manifests as: residuals = X*slopes - intercept\n")
    cat("  This confirms the suspicious line: res = y = x*temp - coef(0)\n")
  }
} else {
  cat("? Residuals don't match expected patterns.\n")
  cat("  Further investigation needed.\n")
  cat("  Showing first 10 values for manual inspection:\n\n")
  comparison <- data.frame(
    returned = result2$residuals[1:10],
    manual = manual_residuals2[1:10],
    fitted = fitted2_manual[1:10],
    y_orig = y_before[1:10]
  )
  print(comparison)
}

cat("\n=======================================================\n")
cat("Additional Diagnostic: Check if y was corrupted\n")
cat("=======================================================\n")
cat("This tests if the chained assignment 'res = y = ...' \n")
cat("corrupted the original y vector in memory.\n\n")

# Note: In R, the function receives a copy, so y_test outside
# shouldn't be modified, but let's check the internal behavior
cat("y values before function call - first 5:", y_before[1:5], "\n")
cat("y values after function call  - first 5:", y_test[1:5], "\n")
cat("Are they identical?", identical(y_before, y_test), "\n\n")

cat("=======================================================\n")
cat("TEST 3: Verify with Simple Hand-Calculable Example\n")
cat("=======================================================\n")

# Create super simple data we can verify by hand
X_simple <- matrix(c(1, 2, 3,
                      0, 1, 2), ncol=2)
y_simple <- c(5, 7, 9)  # Perfect fit: y = 3 + 2*x1 + 0*x2

cat("Simple dataset (should have perfect/near-perfect fit):\n")
cat("X:\n")
print(X_simple)
cat("\ny:", y_simple, "\n\n")

result_simple <- tsreg_C(x=X_simple, y=y_simple)

cat("Estimated coefficients:", result_simple$coef, "\n")
cat("Returned residuals:", result_simple$residuals, "\n\n")

# Compute manual residuals
fitted_simple <- result_simple$coef[1] + X_simple %*% result_simple$coef[-1]
manual_res_simple <- y_simple - fitted_simple

cat("Manual residuals (y - fitted):", manual_res_simple, "\n")
cat("Fitted values:", fitted_simple, "\n\n")

cat("Match check (difference < 0.01):",
    all(abs(result_simple$residuals - manual_res_simple) < 0.01), "\n")
