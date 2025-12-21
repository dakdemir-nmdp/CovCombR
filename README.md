# CovCombR: Combining Incomplete Covariance Matrices with EM

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-1.5.0-blue)](https://github.com/dakdemir-nmdp/CovCombR)
[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue)](https://www.r-project.org/)

## Overview

**CovCombR** combines incomplete covariance matrices from multiple studies into a single, complete covariance estimate using the Expectation-Maximization (EM) algorithm. The method handles arbitrary missing data patterns and returns a positive-definite combined covariance matrix ready for downstream analysis.

**Typical use case:** You have covariance/correlation matrices from K studies, each measuring overlapping but distinct variable sets, and need a unified covariance estimate.

## When to Use CovCombR

**✓ Use CovCombR if you have:**
- Multiple studies with overlapping variable sets
- Covariance or correlation matrices from each study
- Missing data patterns (different variables per study)
- Need for a positive-definite combined estimate

**Minimum requirements:**

| Requirement | Value | Why |
|-------------|-------|-----|
| Studies (K) | ≥ 2 | Need multiple sources to combine |
| Variables (p) | ≥ 3 | Covariance estimation requires variance |
| Connectivity | All variables connected | Ensures identifiability |
| Sample size (ν) | ≥ p_k per study | Statistical requirement for Wishart |
| Overlap | > 50% recommended | Better stability and coverage |

**✗ When NOT to use CovCombR:**
- You have complete data from one study → use standard methods
- Variables are completely disjoint across studies → cannot combine
- Single study with random missingness → use imputation methods
- You need sparse precision matrix estimation → use graphical lasso

See `vignette("study-design-guide")` for detailed design guidance.

## Installation

```r
# Install from CRAN (when available)
install.packages("CovCombR")

# Or install development version from GitHub
# install.packages("devtools")
devtools::install_github("dakdemir-nmdp/CovCombR", build_vignettes = TRUE)
```

## Quick Start

```r
library(CovCombR)

# Simulate 3 studies with overlapping variables
set.seed(2025)
p <- 8
true_Sigma <- diag(p)
true_Sigma[1:4, 1:4] <- 0.6
true_Sigma[5:8, 5:8] <- 0.5
diag(true_Sigma) <- 1

# Study 1: Variables 1-5
S1 <- rWishart(1, 60, true_Sigma[1:5, 1:5])[,,1] / 60

# Study 2: Variables 3-8
S2 <- rWishart(1, 70, true_Sigma[3:8, 3:8])[,,1] / 70

# Study 3: Variables 1,3,5,7
S3 <- rWishart(1, 55, true_Sigma[c(1,3,5,7), c(1,3,5,7)])[,,1] / 55

# Combine with plugin SEs (fast)
fit <- fit_covcomb(
  S_list = list(Study1 = S1, Study2 = S2, Study3 = S3),
  nu = c(Study1 = 60, Study2 = 70, Study3 = 55),
  se_method = "plugin"
)

# Extract combined covariance
Sigma_combined <- coef(fit)
```

## Key Features

- **Handles arbitrary missing patterns** - No restrictions on overlap structure
- **Guaranteed convergence** - Monotonic likelihood increase
- **Standard errors** - Plugin (fast), bootstrap (reliable), or SEM (experimental)
- **Scale estimation** - Optional heterogeneous scale correction
- **Numerically stable** - Cholesky decomposition, ridge regularization

## Main Function Reference

### fit_covcomb()

```r
fit <- fit_covcomb(
  S_list,                  # List of sample covariances
  nu,                      # Degrees of freedom
  scale_method = "none",   # "none" or "estimate"
  se_method = "plugin",    # "none", "plugin", "bootstrap", "sem"
  control = list()         # Advanced options
)
```

### Essential Parameters

**S_list** - Named list of sample covariance matrices
- Each matrix: symmetric, positive semi-definite
- Row/column names identify variables
- Can have different dimensions (incomplete data)
- Example: `list(Study1 = S1, Study2 = S2)`

**nu** - Named vector of degrees of freedom
- Typically: sample size - 1
- Or: effective number of markers (for genomic relationship matrices)
- Must match names in S_list
- Example: `c(Study1 = 100, Study2 = 95)`

**scale_method** - Scale estimation
- `"none"` (default): assume comparable scales across studies
- `"estimate"`: jointly estimate scale factors α_k for each study

**se_method** - Standard error computation

| Method | Speed | Use For | Caveat |
|--------|-------|---------|--------|
| `"none"` | Instant | Exploration only | No SEs computed |
| `"plugin"` | Fast | Quick diagnostics | Anti-conservative |
| `"bootstrap"` | Slow | **Publication** | Most reliable |
| `"sem"` | Medium | Large problems | Experimental, validate |

See `vignette("standard-errors")` for detailed selection guide.

### Return Value

```r
fit$Sigma_hat      # PRIMARY OUTPUT - combined covariance
fit$Sigma_se       # Standard errors (if requested)
fit$S_hat          # Alternative scale (usually identical)
fit$alpha_hat      # Scale factors (if estimated)
fit$convergence    # Convergence diagnostics

coef(fit)          # Extract Sigma_hat
fitted(fit)        # Extract S_hat
```

**Which output to use:** For most applications, use `Sigma_hat` via `coef(fit)`. See `?fit_covcomb` for complete documentation.

## Complete Workflow Example

```r
library(CovCombR)

# Realistic scenario: 3 genomic studies with partial overlap
set.seed(2025)
p <- 12
var_names <- paste0("SNP", 1:p)

# True genetic relationship matrix
true_GRM <- diag(p)
true_GRM[1:6, 1:6] <- 0.4  # Related individuals block 1
true_GRM[7:12, 7:12] <- 0.35  # Related individuals block 2
diag(true_GRM) <- 1
dimnames(true_GRM) <- list(var_names, var_names)

# Study 1: SNPs 1-8 (1000 markers)
GRM1 <- rWishart(1, 1000, true_GRM[1:8, 1:8])[,,1] / 1000
dimnames(GRM1) <- list(var_names[1:8], var_names[1:8])

# Study 2: SNPs 5-12 (1200 markers)
GRM2 <- rWishart(1, 1200, true_GRM[5:12, 5:12])[,,1] / 1200
dimnames(GRM2) <- list(var_names[5:12], var_names[5:12])

# Study 3: SNPs 1-4, 9-12 (900 markers)
idx3 <- c(1:4, 9:12)
GRM3 <- rWishart(1, 900, true_GRM[idx3, idx3])[,,1] / 900
dimnames(GRM3) <- list(var_names[idx3], var_names[idx3])

# Combine with bootstrap SEs for publication
fit <- fit_covcomb(
  S_list = list(Study1 = GRM1, Study2 = GRM2, Study3 = GRM3),
  nu = c(Study1 = 1000, Study2 = 1200, Study3 = 900),
  se_method = "bootstrap",
  control = list(
    bootstrap = list(B = 200, seed = 2025, progress = TRUE)
  )
)

# Extract results
GRM_combined <- coef(fit)
GRM_se <- fit$Sigma_se

# Hypothesis test: Is correlation between SNP1-SNP10 significant?
# (These were NEVER observed together in any study!)
cor_estimate <- GRM_combined[1, 10]
cor_se <- GRM_se[1, 10]
z_stat <- cor_estimate / cor_se
p_value <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)

cat(sprintf("Correlation: %.3f ± %.3f\n", cor_estimate, cor_se))
cat(sprintf("p-value: %.4f\n", p_value))
```

## Documentation

Comprehensive vignettes cover all aspects:

```r
# Getting started
vignette("CovCombR-vignette")      # Quick introduction
vignette("iris-example")           # Beginner-friendly example

# Choosing methods
vignette("standard-errors")        # SE method selection guide
vignette("study-design-guide")     # Design recommendations

# Advanced topics
vignette("advanced-configuration") # Scale parameters, control options
vignette("combining-grms")         # Genomic applications
vignette("statistical-methods")    # Mathematical foundations

# Reference
?fit_covcomb                       # Complete function documentation
```

## Common Workflows

### Exploratory Analysis
```r
# Quick fit without SEs
fit <- fit_covcomb(S_list, nu, se_method = "none")
print(fit)
```

### Publication-Quality Inference
```r
# Bootstrap SEs (recommended)
fit <- fit_covcomb(S_list, nu, se_method = "bootstrap",
                  control = list(bootstrap = list(B = 200)))
```

### Heterogeneous Scales (e.g., multi-platform GRMs)
```r
# Estimate scale factors
fit <- fit_covcomb(S_list, nu, scale_method = "estimate")
print(fit$alpha_hat)  # Relative scales
```

## Statistical Background

Given K studies observing overlapping variable subsets, the model assumes:

$$\mathbf{S}_k \sim \text{Wishart}(\nu_k, \alpha_k \boldsymbol{\Sigma}_{\mathcal{O}_k})$$

where:
- $\mathbf{S}_k$ = sample covariance from study k
- $\nu_k$ = degrees of freedom (sample size - 1)
- $\alpha_k$ = optional scale factor
- $\boldsymbol{\Sigma}$ = common underlying covariance (estimated)
- $\mathcal{O}_k$ = observed variable indices in study k

The EM algorithm iterates between:
- **E-step**: Compute conditional expectations of missing blocks
- **M-step**: Update covariance (and optionally scale factors)

See `vignette("statistical-methods")` for complete mathematical details.

## Performance

**Computational complexity:** O(K × p³) per iteration

**Typical convergence:** 20-100 iterations

**For large problems (p > 100):**
- Start with `se_method = "none"` to assess convergence
- Use `init_sigma = "identity"` for faster initialization
- Consider `se_method = "sem"` (validate against bootstrap)

See `vignette("optimization")` for large-scale strategies.

## Troubleshooting

**Convergence issues?**
- Increase `max_iter` or relax `tol`
- Try `init_sigma = "identity"`
- Increase ridge: `control = list(ridge = 1e-6)`

**Non-positive definite estimates?**
- Check overlap: all variables connected?
- Increase ridge regularization
- Use `scale_method = "estimate"` if heterogeneous scaling

**High condition number warning?**
- Indicates weak identifiability
- Check connectivity with `vignette("study-design-guide")`
- May need more overlap or studies

See `vignette("troubleshooting")` for comprehensive diagnostics.

## Contributing

Contributions are welcome! Please submit issues or pull requests at:
https://github.com/dakdemir-nmdp/CovCombR

## License

MIT License. See [LICENSE](LICENSE) for details.

## Citation

If you use CovCombR in your research, please cite:

```
Akdemir, D. (2025). CovCombR: Combining Incomplete Covariance Matrices with EM.
R package version 1.5.0. https://github.com/dakdemir-nmdp/CovCombR
```

## Contact

- **Author:** Deniz Akdemir
- **Email:** dakdemir@nmdp.org
- **GitHub:** https://github.com/dakdemir-nmdp/CovCombR
- **Issues:** https://github.com/dakdemir-nmdp/CovCombR/issues

---

**Quick Links:**
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Documentation](#documentation)
- [Vignettes](https://dakdemir-nmdp.github.io/CovCombR/)
