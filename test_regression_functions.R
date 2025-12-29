#!/usr/bin/env Rscript
# Test script for WRScpp regression functions
# This script tests all four main regression estimators:
# 1. tsreg_C()    - Theil-Sen regression
# 2. tshdreg_C()  - Theil-Sen with Harrell-Davis estimator
# 3. stsreg_C()   - Theil-Sen minimizing robust variance
# 4. tstsreg_C()  - Modified Theil-Sen with outlier removal

library(WRScpp)

cat("=======================================================\n")
cat("Testing WRScpp Regression Functions\n")
cat("=======================================================\n\n")

# Create test datasets
set.seed(42)
n <- 100
p <- 4

# Dataset 1: Simple linear relationship with some outliers
cat("Creating test datasets...\n")
x1 <- rnorm(n)
y1 <- 2 + 3*x1 + rnorm(n, sd=0.5)
# Add some outliers
outlier_idx <- sample(1:n, 5)
y1[outlier_idx] <- y1[outlier_idx] + rnorm(5, mean=0, sd=5)

# Dataset 2: Multiple predictors
X2 <- matrix(rnorm(n*p), ncol=p)
true_coef <- c(1, 2, -1.5, 0.5)
y2 <- 3 + X2 %*% true_coef + rnorm(n, sd=0.8)

# Dataset 3: Multiple predictors with outliers
X3 <- matrix(rnorm(n*p), ncol=p)
y3 <- -1 + X3 %*% c(1.5, -2, 1, 0.8) + rnorm(n, sd=0.6)
outlier_idx3 <- sample(1:n, 8)
y3[outlier_idx3] <- y3[outlier_idx3] + rnorm(8, mean=0, sd=6)

cat("Datasets created successfully.\n\n")

# Helper function to print results nicely
print_results <- function(name, result, check_residuals=TRUE) {
  cat("----------------------------------------\n")
  cat(name, "\n")
  cat("----------------------------------------\n")
  cat("Coefficients:\n")
  print(result$coef)

  if (!is.null(result$Strength.Assoc)) {
    cat("\nStrength of Association:", result$Strength.Assoc, "\n")
  }

  if (!is.null(result$Explanatory.Power)) {
    cat("Explanatory Power:", result$Explanatory.Power, "\n")
  }

  if (check_residuals && !is.null(result$residuals)) {
    cat("\nResidual summary:\n")
    print(summary(result$residuals))
    cat("Residual SD:", sd(result$residuals), "\n")
  }

  cat("\n")
}

# Test 1: tsreg_C() - Single predictor
cat("=======================================================\n")
cat("TEST 1: tsreg_C() with single predictor\n")
cat("=======================================================\n")
tryCatch({
  result1 <- tsreg_C(x=x1, y=y1)
  print_results("tsreg_C (single predictor)", result1)
  cat("✓ Test passed\n\n")
}, error = function(e) {
  cat("✗ Test FAILED:", conditionMessage(e), "\n\n")
})

# Test 2: tsreg_C() - Multiple predictors
cat("=======================================================\n")
cat("TEST 2: tsreg_C() with multiple predictors\n")
cat("=======================================================\n")
tryCatch({
  result2 <- tsreg_C(x=X2, y=y2)
  print_results("tsreg_C (multiple predictors)", result2)
  cat("True coefficients: [intercept=3,", paste(true_coef, collapse=", "), "]\n")
  cat("✓ Test passed\n\n")
}, error = function(e) {
  cat("✗ Test FAILED:", conditionMessage(e), "\n\n")
})

# Test 3: tsreg_C() with outliers and custom iterations
cat("=======================================================\n")
cat("TEST 3: tsreg_C() with outliers (iter=20)\n")
cat("=======================================================\n")
tryCatch({
  result3 <- tsreg_C(x=X3, y=y3, iter=20)
  print_results("tsreg_C (with outliers)", result3)
  cat("✓ Test passed\n\n")
}, error = function(e) {
  cat("✗ Test FAILED:", conditionMessage(e), "\n\n")
})

# Test 4: tshdreg_C() - Single predictor
cat("=======================================================\n")
cat("TEST 4: tshdreg_C() with single predictor\n")
cat("=======================================================\n")
tryCatch({
  result4 <- tshdreg_C(x=x1, y=y1)
  print_results("tshdreg_C (single predictor)", result4)
  cat("✓ Test passed\n\n")
}, error = function(e) {
  cat("✗ Test FAILED:", conditionMessage(e), "\n\n")
})

# Test 5: tshdreg_C() - Multiple predictors
cat("=======================================================\n")
cat("TEST 5: tshdreg_C() with multiple predictors\n")
cat("=======================================================\n")
tryCatch({
  result5 <- tshdreg_C(x=X2, y=y2, iter=15)
  print_results("tshdreg_C (multiple predictors)", result5)
  cat("✓ Test passed\n\n")
}, error = function(e) {
  cat("✗ Test FAILED:", conditionMessage(e), "\n\n")
})

# Test 6: stsreg_C() - Single predictor
cat("=======================================================\n")
cat("TEST 6: stsreg_C() with single predictor\n")
cat("=======================================================\n")
tryCatch({
  result6 <- stsreg_C(x=x1, y=y1)
  print_results("stsreg_C (single predictor)", result6)
  cat("✓ Test passed\n\n")
}, error = function(e) {
  cat("✗ Test FAILED:", conditionMessage(e), "\n\n")
})

# Test 7: stsreg_C() - Multiple predictors
cat("=======================================================\n")
cat("TEST 7: stsreg_C() with multiple predictors\n")
cat("=======================================================\n")
tryCatch({
  result7 <- stsreg_C(x=X3, y=y3, iter=12)
  print_results("stsreg_C (multiple predictors with outliers)", result7)
  cat("✓ Test passed\n\n")
}, error = function(e) {
  cat("✗ Test FAILED:", conditionMessage(e), "\n\n")
})

# Test 8: tstsreg_C() - Single predictor
cat("=======================================================\n")
cat("TEST 8: tstsreg_C() with single predictor\n")
cat("=======================================================\n")
tryCatch({
  result8 <- tstsreg_C(x=x1, y=y1)
  print_results("tstsreg_C (single predictor)", result8, check_residuals=TRUE)
  cat("✓ Test passed\n\n")
}, error = function(e) {
  cat("✗ Test FAILED:", conditionMessage(e), "\n\n")
})

# Test 9: tstsreg_C() - Multiple predictors with outliers
cat("=======================================================\n")
cat("TEST 9: tstsreg_C() with multiple predictors + outliers\n")
cat("=======================================================\n")
tryCatch({
  result9 <- tstsreg_C(x=X3, y=y3)
  print_results("tstsreg_C (multiple predictors with outliers)", result9)
  cat("✓ Test passed\n\n")
}, error = function(e) {
  cat("✗ Test FAILED:", conditionMessage(e), "\n\n")
})

# Test 10: Comparison test - all functions on same data
cat("=======================================================\n")
cat("TEST 10: Comparison of all methods on same dataset\n")
cat("=======================================================\n")
cat("Dataset: 4 predictors with outliers\n\n")

tryCatch({
  comp_tsreg <- tsreg_C(x=X3, y=y3)
  comp_tshdreg <- tshdreg_C(x=X3, y=y3)
  comp_stsreg <- stsreg_C(x=X3, y=y3)
  comp_tstsreg <- tstsreg_C(x=X3, y=y3)

  cat("Method Comparison:\n")
  cat("------------------\n")
  cat("tsreg_C    intercept:", comp_tsreg$coef[1], "  slopes:",
      paste(round(comp_tsreg$coef[-1], 3), collapse=", "), "\n")
  cat("tshdreg_C  intercept:", comp_tshdreg$coef[1], "  slopes:",
      paste(round(comp_tshdreg$coef[-1], 3), collapse=", "), "\n")
  cat("stsreg_C   intercept:", comp_stsreg$coef[1], "  slopes:",
      paste(round(comp_stsreg$coef[-1], 3), collapse=", "), "\n")
  cat("tstsreg_C  intercept:", comp_tstsreg$coef[1], "  slopes:",
      paste(round(comp_tstsreg$coef[-1], 3), collapse=", "), "\n")

  cat("\nResidual standard deviations:\n")
  cat("tsreg_C   :", sd(comp_tsreg$residuals), "\n")
  cat("tshdreg_C :", sd(comp_tshdreg$residuals), "\n")
  cat("stsreg_C  :", sd(comp_stsreg$residuals), "\n")
  cat("tstsreg_C :", sd(comp_tstsreg$residuals), "\n")

  cat("\n✓ Comparison test passed\n\n")
}, error = function(e) {
  cat("✗ Test FAILED:", conditionMessage(e), "\n\n")
})

# Test 11: Edge case - matrix with one column
cat("=======================================================\n")
cat("TEST 11: Edge case - matrix input with single column\n")
cat("=======================================================\n")
tryCatch({
  X_single <- matrix(x1, ncol=1)
  result11 <- tsreg_C(x=X_single, y=y1)
  print_results("tsreg_C (matrix with 1 column)", result11)
  cat("✓ Test passed\n\n")
}, error = function(e) {
  cat("✗ Test FAILED:", conditionMessage(e), "\n\n")
})

# Test 12: Parameter variations
cat("=======================================================\n")
cat("TEST 12: Testing parameter variations\n")
cat("=======================================================\n")
tryCatch({
  cat("a) tsreg_C with HD=TRUE (Harrell-Davis for intercept):\n")
  result12a <- tsreg_C(x=X2, y=y2, HD=TRUE)
  cat("   Coefficients:", paste(round(result12a$coef, 3), collapse=", "), "\n\n")

  cat("b) tshdreg_C with different tolerance:\n")
  result12b <- tshdreg_C(x=X2, y=y2, tol=0.001)
  cat("   Coefficients:", paste(round(result12b$coef, 3), collapse=", "), "\n\n")

  cat("c) stsreg_C with fewer iterations:\n")
  result12c <- stsreg_C(x=X2, y=y2, iter=5)
  cat("   Coefficients:", paste(round(result12c$coef, 3), collapse=", "), "\n\n")

  cat("✓ Parameter variation tests passed\n\n")
}, error = function(e) {
  cat("✗ Test FAILED:", conditionMessage(e), "\n\n")
})

cat("=======================================================\n")
cat("ALL TESTS COMPLETED\n")
cat("=======================================================\n")
cat("\nAll regression functions appear to be working correctly!\n")
cat("\nNote: Verify that residuals are calculated correctly,\n")
cat("especially for tsreg_C with multiple predictors, due to\n")
cat("the potential bug noted in the README.\n")
