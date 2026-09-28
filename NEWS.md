# CovCombR 1.7.0

## Identifiability reporting

* **New `identifiability_report()`.** Determines, from the observation design
  alone and without fitting, which entries of the combined covariance matrix
  the data can determine. A variable pair never jointly observed in any sample
  carries no information about its covariance: under the free-Sigma model the
  observed-data log-likelihood is exactly flat in that coordinate, so the value
  the EM algorithm returns there is a deterministic function of `init_sigma`
  rather than of the data.

* **`fit_covcomb()` now warns** when a free-Sigma fit leaves entries
  unidentified, naming the count, the proportion, and example pairs. Under a
  factor model it emits a message instead, since the factor structure is a
  genuine identifying assumption; the point there is that such entries rest on
  the assumption rather than on direct observation, and should be reported that
  way.

* **Results carry an `$identifiability` component** and `summary()` reports it.

* This closes a real gap. The previously documented requirement that the
  free-Sigma model "is only identifiable when every variable pair is jointly
  observed in at least one study" was never checked at runtime. The existing
  graph-connectivity check is strictly weaker: a chain design observing
  variables 1-8, 6-15 and 13-20 is fully connected, yet 95 of its 190
  off-diagonal entries (50%) are never jointly observed. Such fits previously
  returned initialization-determined values silently.

* One boundary case is documented but not detected: if either conditional
  variance entering the feasible interval is exactly zero, the entry is
  identified after all, forced by positive-definiteness. That case is not
  generic.

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

# CovCombR 1.4.0

(Released as WishartEM 1.4.0.)

* Enhanced bootstrap standard error computation
* Improved convergence diagnostics
* Added vignettes for statistical methods

# CovCombR 1.3.0

(Released as WishartEM 1.3.0; initial CRAN release.)

* Core EM algorithm implementation
* Support for heterogeneous scaling
