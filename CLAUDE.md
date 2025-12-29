# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

WRScpp is a pre-compiled R package that provides C++ sub-routines for the R package WRS (Robust Statistics) by Dr. Rand Wilcox. This repository contains a **binary distribution** of the package, not the source code.

**Important**: This repository contains the compiled package (libs/WRScpp.so) but does NOT contain the C++ source files (.cpp/.h). The actual source code is maintained in a separate repository: https://github.com/mrxiaohe/robustmethods_cplusplus

## Package Structure

This is an installed R package with the following structure:
- `libs/WRScpp.so` - Pre-compiled shared library (macOS x86_64)
- `R/WRScpp` - Lazy-loaded R functions database
- `DESCRIPTION` - Package metadata
- `NAMESPACE` - Package namespace configuration
- `Meta/` - Package metadata (RDS files)
- `help/` - Help documentation database
- `html/` - HTML documentation

## Dependencies

Required R packages (from DESCRIPTION):
- Rcpp (>= 0.10.2)
- RcppArmadillo (>= 0.3.900.0)
- MASS (>= 7.3-27)
- scatterplot3d (>= 0.3-33)

## Installing and Using the Package

To install from this repository:
```r
library("devtools")
install_github(repo="WRScpp", username="twobunnyrabbit")
```

Or to install dependencies and the WRS package:
```r
install.packages("WRS", repos="http://R-Forge.R-project.org", type="source")
install.packages(c("RcppArmadillo", "devtools"))
```

## Main Functions

The package provides C++-accelerated versions of robust regression estimators:

1. **tsreg_C()** - Theil-Sen regression estimator using Gauss-Seidel algorithm for multiple predictors
2. **tshdreg_C()** - Theil-Sen variant using Harrell-Davis estimator instead of sample median
3. **stsreg_C()** - Theil-Sen variant that minimizes robust variance of residuals (percentage bend midvariance)
4. **tstsreg_C()** - Modified Theil-Sen that uses stsreg_C() for initial estimates, detects outliers, then runs regular Theil-Sen on cleaned data

All functions provide significant performance improvements over pure R implementations (typically 5-30x faster).

## Platform Notes

- This repository contains a macOS x86_64 binary (libs/WRScpp.so is Mach-O format)
- For 64-bit Linux: https://github.com/JoeJohnston/WRScppLin64
- For 64-bit Windows: http://github.com/mrxiaohe/WRScppWIN
- The compiled library is platform-specific and cannot be used across different OS/architectures

## Development Workflow

Since this is a binary distribution repository:

1. **To modify the package**: Work in the source code repository at https://github.com/mrxiaohe/robustmethods_cplusplus
2. **This repository is for distribution**: Contains pre-built binaries for macOS users to install via devtools::install_github()
3. **No compilation happens here**: The .so file is pre-compiled and committed to version control
4. **To update the binary**: Rebuild from source repository and commit the new .so file

## Testing

Example usage patterns are documented in README.md, including:
- Creating test datasets
- Comparing performance with pure R implementations using `rbenchmark` or `system.time()`
- Visualizing regression lines with base R plotting

No automated test suite is present in this repository.
