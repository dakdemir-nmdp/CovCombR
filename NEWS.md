# CovCombR 1.6.0

## Highlights

* **Factor-analytic GRM combination is now the primary documented workflow.**
  Documentation, README, and DESCRIPTION have been rewritten to foreground the
  FA model (Σ = ΛΛ⊤ + Ψ) as the recommended approach for combining incomplete
  genomic relationship matrices across multi-platform genotyping scenarios.

## New Vignettes

* `combining-grms-factor-model`: Comprehensive vignette demonstrating FA-model
  GRM combination with the BGLR wheat dataset, including:
  - Chain-overlap cohort design with unobserved pairs
  - Model comparison by BIC (2-factor, 3-factor, 5-factor, free)
  - Recovery metrics separated by observed vs. unobserved pairs
  - Downstream genomic prediction (GBLUP) with the combined GRM

## Documentation

* README rewritten to stress FA methods for GRM combination as the key
  differentiator — particularly the ability to predict relatedness for
  individual pairs never jointly observed.
* DESCRIPTION updated: title and description now emphasize factor-analytic
  models and genomic relationship matrices.
* Added BGLR, ggplot2, reshape2, gridExtra to Suggests for the new vignette.

---

# CovCombR 1.5.0

## Breaking Changes

* **Package renamed from WishartEM to CovCombR**
* **Main function renamed**: `fit_wishart_em()` → `fit_covcomb()`
* **S3 class renamed**: `wishart_em` → `covcomb`
* All S3 methods updated accordingly: `print.covcomb()`, `summary.covcomb()`, `coef.covcomb()`, `fitted.covcomb()`

## Migration Guide

### Old code (WishartEM):
```r
library(WishartEM)
result <- fit_wishart_em(S_list, nu, se_method = "plugin")
class(result)  # "wishart_em"
```

### New code (CovCombR):
```r
library(CovCombR)
result <- fit_covcomb(S_list, nu, se_method = "plugin")
class(result)  # "covcomb"
```

All other functionality remains unchanged.

---

# Previous Versions (as WishartEM)

## WishartEM 1.4.0

* Enhanced bootstrap standard error computation
* Improved convergence diagnostics
* Added vignettes for statistical methods

## WishartEM 1.3.0

* Initial CRAN release
* Core EM algorithm implementation
* Support for heterogeneous scaling
