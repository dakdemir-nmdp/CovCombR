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
