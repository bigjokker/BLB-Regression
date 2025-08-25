# BLB Regression

## Overview

This R script (`blb_regression.R`) implements the Bag of Little Bootstraps (BLB) method for linear regression. BLB is designed for massive datasets where fitting models or bootstrapping on the full data is computationally infeasible. It provides bagged point estimates and percentile-t confidence intervals that approximate full-data results, using only small subsample fits.

The core function is `pooled_t_blb_multi`, which supports robust (HC1) standard errors, factor variables, and ridge regularization for stability.

## Requirements

- R (base installation sufficient)
- Optional: `MASS` package for generalized inverse in singular cases (falls back to eigen decomposition if unavailable)

## Installation

Simply source the script in your R session:

```r
source("blb_regression.R")
```

## Usage

The function signature:

```r
pooled_t_blb_multi(
  formula,          # Model formula (e.g., y ~ X1 + X2)
  data,             # Data frame with the full large dataset
  b = NULL,         # Subsample size (default: max(1000, N^0.7))
  s = 20,           # Number of outer subsamples
  r = 50,           # Number of inner bootstraps per subsample
  robust = TRUE,    # Use HC1 robust SEs (FALSE for classical)
  ridge = 1e-10,    # Small ridge for matrix inversion stability
  seed = NULL,      # Optional seed for reproducibility
  contrasts.arg = NULL  # Optional contrasts for factors
)
```

Returns a list with:
- `est`: Bagged coefficient estimates
- `ci`: 95% confidence intervals (matrix with lower/upper rows)
- `pooled_q`: Pooled t-quantiles
- `se_global`: Global SE estimates
- `settings`: Input parameters

### Example

```r
set.seed(123)
N <- 1e6
X <- rnorm(N)
y <- 2 + 3 * X + rnorm(N)
df <- data.frame(y, X)

res <- pooled_t_blb_multi(y ~ X, data = df, b = 5000, s = 20, r = 50, robust = FALSE, seed = 123)
res$est  # Bagged estimates
res$ci   # 95% CIs (approximates full-data tightness)
```

For heteroskedastic errors, set `robust = TRUE`.

## Notes

- Subsamples are drawn without replacement; inner bootstraps use multinomial weights summing to N.
- Handles NA values via `na.omit` per subsample.
- For very large N, ensure data fits in memory (only indexing, no full matrix ops).
- Inspired by the Bag of Little Bootstraps method for scalable inference.
