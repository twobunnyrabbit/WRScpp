# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

WRScpp is an R package that provides C++ sub-routines for the R package WRS (Robust Statistics) by Dr. Rand Wilcox. This repository contains the **source code** and can be compiled on Linux and macOS (both Intel and Apple Silicon).

## Package Structure

This is a standard R source package with the following structure:
- `src/robustmethods_CPP.cpp` - C++ source code (main implementation)
- `src/Makevars` - Build configuration for Linux/macOS
- `src/WRScpp.so` - Compiled shared library (created during build)
- `R/` - R wrapper functions
- `DESCRIPTION` - Package metadata
- `NAMESPACE` - Package namespace configuration
- `man/` - Documentation files

## Dependencies

Required R packages (from DESCRIPTION):
- Rcpp (>= 0.10.2)
- RcppArmadillo (>= 0.3.900.0)
- MASS (>= 7.3-27)
- scatterplot3d (>= 0.3-33)

## Installing and Using the Package

### Prerequisites

**For macOS M1/M2/M3 (Apple Silicon):**
You must install gfortran for ARM64 before installing this package:
1. Download from: https://mac.r-project.org/tools/
2. Get "gfortran for arm64 (Apple Silicon)" and install the .pkg file
3. Or use Homebrew: `brew install gcc`

**For Intel Macs and Linux:**
Standard development tools should be sufficient.

### Installation

To install from GitHub:
```r
library("devtools")
install_github(repo="WRScpp", username="twobunnyrabbit")
```

To install from local directory:
```r
devtools::install_local(".")
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

- **macOS (Intel & Apple Silicon)**: Requires gfortran (see Prerequisites above)
- **Linux**: Should compile with standard development tools (gcc, g++, gfortran)
- **Windows**: May require additional configuration; see http://github.com/mrxiaohe/WRScppWIN for Windows-specific builds
- The compiled `.so` file is platform-specific and will be rebuilt during installation

## Development Workflow

This is a standard R source package:

1. **Source code**: Main C++ implementation is in `src/robustmethods_CPP.cpp`
2. **Build system**: Uses standard R package build tools with Rcpp/RcppArmadillo
3. **Compilation**: The `.so` library is compiled during package installation
4. **Build configuration**: `src/Makevars` controls compilation flags and linking

### Making Changes

1. Modify C++ code in `src/robustmethods_CPP.cpp`
2. Update R wrappers in `R/` if needed
3. Rebuild and test: `devtools::load_all()` or `devtools::install_local(".")`
4. The package will automatically recompile when installed

## Testing

Example usage patterns are documented in README.md, including:
- Creating test datasets
- Comparing performance with pure R implementations using `rbenchmark` or `system.time()`
- Visualizing regression lines with base R plotting

No automated test suite is present in this repository.
