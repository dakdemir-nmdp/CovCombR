# CovCombR: Combining Incomplete Genomic Relationship Matrices via Factor-Analytic EM

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-1.6.0-blue)](https://github.com/dakdemir-nmdp/CovCombR)
[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue)](https://www.r-project.org/)

## Overview

**CovCombR** combines incomplete covariance matrices — especially **Genomic Relationship Matrices (GRMs)** — from multiple studies into a single, complete estimate using the Expectation-Maximization (EM) algorithm. Its core strength is a **factor-analytic (FA) model** that can predict pairwise relatedness *even for individuals that never appeared in the same study*.

**Key capability:** When cohorts are genotyped on different SNP platforms with only partial individual overlap, the FA model (Σ = ΛΛ⊤ + Ψ) leverages shared anchor individuals to infer relationships between individuals that were never jointly observed — something no unconstrained model can do.

### Why Factor-Analytic Models for GRMs?

GRMs are inherently low-rank: population structure is driven by a small number of principal axes. A factor model exploits this by estimating per-individual loadings (Λ) on k latent factors, so that:

- If individual *i*'s loading is determined via shared anchors in Cohort 1, and
- individual *j*'s loading is determined via shared anchors in Cohort 2,
- then their relatedness σᵢⱼ = Λ[i,·] · Λ[j,·] is **identified** — even though *i* and *j* never appeared in the same GRM.

This makes CovCombR uniquely suited for:

- **Multi-platform genomic prediction** across cohorts genotyped on different SNP chips
- **Meta-analysis of breeding values** from independent trials with partial individual overlap
- **Imputation-free workflows** where raw genotypes are unavailable but per-platform GRMs are

## Installation

```r
# Install development version from GitHub
# install.packages("devtools")
devtools::install_github("dakdemir-nmdp/CovCombR", build_vignettes = TRUE)
```

## Quick Start: Combining GRMs with a Factor Model

```r
library(CovCombR)

# Three cohorts genotyped on different SNP platforms
# Cohort 1: individuals 1-20 on Chip A
# Cohort 2: individuals 16-35 on Chip B (overlap: 16-20)
# Cohort 3: individuals 31-40 + 1-5 on Chip C (overlap: 31-35, 1-5)
# → Many pairs (e.g., individual 10 & 38) are NEVER jointly observed

# Combine with a 3-factor model
fit <- fit_covcomb(
  S_list = list(GRM1 = G1, GRM2 = G2, GRM3 = G3),
  nu = c(GRM1 = 799, GRM2 = 878, GRM3 = 878),
  n_factors = 3,
  se_method = "plugin"
)

# Combined GRM — including predicted entries for unobserved pairs
GRM_combined <- coef(fit)

# Factor loadings reveal population structure
fit$Lambda_hat   # p x 3 loadings matrix
fit$Psi_hat      # per-individual unique variances
```

### Choosing the Number of Factors

```r
# Compare models by BIC
fit_2f <- fit_covcomb(S_list, nu, n_factors = 2, se_method = "none")
fit_3f <- fit_covcomb(S_list, nu, n_factors = 3, se_method = "none")
fit_5f <- fit_covcomb(S_list, nu, n_factors = 5, se_method = "none")

# Or use automatic selection (best for small problems)
fit_auto <- fit_covcomb(S_list, nu, n_factors = "auto")
```

**Rule of thumb:** Set k equal to the number of major population subgroups. For large GRMs, manually compare 2–3 values of k via BIC rather than using `n_factors = "auto"`.

## When to Use CovCombR

**✓ Use CovCombR if you have:**
- Multiple GRMs (or covariance matrices) with overlapping but non-identical individuals/variables
- Cohorts genotyped on different SNP platforms
- Need to predict relatedness for pairs never jointly observed (**→ use FA model**)
- Need a unified positive-definite GRM for downstream genomic prediction

**Minimum requirements:**

| Requirement | Value | Why |
|-------------|-------|-----|
| Studies (K) | ≥ 2 | Need multiple sources to combine |
| Variables (p) | ≥ 3 | Covariance estimation requires variance |
| Connectivity | All variables connected | Ensures identifiability |
| Sample size (ν) | ≥ p_k per study | Statistical requirement for Wishart |
| Overlap | > 50% recommended | Better stability and coverage |

**✗ When NOT to use CovCombR:**
- You have a complete GRM from one study → use standard methods
- Individuals are completely disjoint across studies → cannot combine (no anchors)
- Single study with random missingness → use imputation methods
- You need sparse precision matrix estimation → use graphical lasso

## Two Model Types

CovCombR offers two approaches for combining matrices:

### 1. Factor-Analytic Model (Recommended for GRMs)

```r
# Σ = ΛΛ⊤ + Ψ  (k latent factors, diagonal unique variances)
fit <- fit_covcomb(S_list, nu, n_factors = 3)
```

**Advantages:**
- ✅ Predicts unobserved pairwise relationships via shared factor loadings
- ✅ Parsimonious: fewer parameters than free model
- ✅ Biologically interpretable (factors ≈ population structure axes)
- ✅ Regularized by construction (low-rank + diagonal)

**When to use:** Multi-platform GRM combination, breeding programs, any setting with unobserved pairs.

### 2. Free (Unconstrained) Model

```r
# Σ estimated freely (no structural constraint)
fit <- fit_covcomb(S_list, nu, n_factors = NULL)
```

**Advantages:**
- ✅ No assumption about covariance structure
- ✅ Appropriate when all pairs are observed in at least one study

**Limitation:** Cannot reliably estimate entries for pairs never jointly observed.

## Main Function Reference

### fit_covcomb()

```r
fit <- fit_covcomb(
  S_list,                  # List of sample covariances / GRMs
  nu,                      # Degrees of freedom (n_markers - 1 for GRMs)
  n_factors = "auto",      # "auto" (default), integer k, or NULL/"free"
  scale_method = "none",   # "none" or "estimate"
  se_method = "plugin",    # "none", "plugin", "bootstrap", "sem"
  control = list()         # Advanced options
)
```

### Key Parameters

**n_factors** — Model selection (NEW)
- `"auto"` (default): fits k = 1, 2, … plus the free model and selects the best by BIC
- Integer `k`: k-factor model Σ = ΛΛ⊤ + Ψ
- `NULL` or `"free"`: unconstrained (free) model

**S_list** — Named list of sample covariance matrices / GRMs
- Row/column names identify individuals (or variables)
- Can have different dimensions (partial overlap)

**nu** — Named vector of degrees of freedom
- For GRMs: number of markers - 1
- For covariance matrices: sample size - 1

**se_method** — Standard error computation

| Method | Speed | Use For | Caveat |
|--------|-------|---------|--------|
| `"none"` | Instant | Exploration only | No SEs computed |
| `"plugin"` | Fast | Quick diagnostics | Anti-conservative |
| `"bootstrap"` | Slow | **Publication** | Most reliable |
| `"sem"` | Medium | Large problems | Experimental, validate |

### Return Value

```r
# Common outputs
fit$Sigma_hat      # Combined covariance / GRM
fit$Sigma_se       # Standard errors (if requested)
fit$convergence    # Convergence diagnostics

# Factor model outputs
fit$Lambda_hat     # p × k factor loadings (varimax-rotated)
fit$Psi_hat        # Per-variable unique variances
fit$n_factors      # Number of factors used
fit$model          # "factor" or "free"

# Scale estimation
fit$alpha_hat      # Scale factors (if scale_method = "estimate")

coef(fit)          # Extract Sigma_hat
fitted(fit)        # Extract S_hat
```

## Complete Workflow: Multi-Platform GRM Combination

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

# Combine with 3-factor model and bootstrap SEs
fit <- fit_covcomb(
  S_list = list(Study1 = GRM1, Study2 = GRM2, Study3 = GRM3),
  nu = c(Study1 = 1000, Study2 = 1200, Study3 = 900),
  n_factors = 3,
  se_method = "bootstrap",
  control = list(bootstrap = list(B = 200, seed = 2025, progress = TRUE))
)

# Combined GRM — includes predicted entries for unobserved pairs
GRM_combined <- coef(fit)
GRM_se <- fit$Sigma_se

# Factor loadings reveal population structure
print(round(fit$Lambda_hat, 3))

# Test: Is relatedness between SNP1-SNP10 significant?
# (These were NEVER observed together in any study!)
est <- GRM_combined[1, 10]
se <- GRM_se[1, 10]
z <- est / se
p_val <- 2 * pnorm(abs(z), lower.tail = FALSE)

cat(sprintf("Predicted relatedness: %.3f ± %.3f\n", est, se))
cat(sprintf("p-value: %.4f\n", p_val))
```

## Documentation

Comprehensive vignettes cover all aspects:

```r
# Getting started
vignette("CovCombR-vignette")                # Quick introduction
vignette("iris-example")                     # Beginner-friendly example

# Factor-analytic GRM combination
vignette("combining-grms-factor-model")      # ★ FA models for GRMs with unobserved pairs

# Choosing methods
vignette("standard-errors")                  # SE method selection guide
vignette("study-design-guide")               # Design recommendations

# Advanced topics
vignette("advanced-configuration")           # Scale parameters, control options
vignette("combining-grms")                   # General genomic applications
vignette("statistical-methods")              # Mathematical foundations

# Reference
?fit_covcomb                                 # Complete function documentation
```

## Common Workflows

### Factor-Analytic GRM Combination (Primary Use Case)
```r
# Combine partial GRMs with a k-factor model
fit <- fit_covcomb(S_list, nu, n_factors = 3, se_method = "plugin")
GRM_combined <- coef(fit)
Lambda <- fit$Lambda_hat  # Inspect population structure
```

### Exploratory Analysis
```r
# Quick fit without SEs
fit <- fit_covcomb(S_list, nu, se_method = "none")
print(fit)
```

### Publication-Quality Inference
```r
# Bootstrap SEs (recommended)
fit <- fit_covcomb(S_list, nu, n_factors = 3, se_method = "bootstrap",
                  control = list(bootstrap = list(B = 200)))
```

### Heterogeneous Scales (e.g., Multi-Platform GRMs)
```r
# Estimate scale factors when SNP chips differ in marker density
fit <- fit_covcomb(S_list, nu, n_factors = 3, scale_method = "estimate")
print(fit$alpha_hat)  # Relative scales per platform
```

## Statistical Background

### Free Model

Given K studies observing overlapping variable subsets:

$$\mathbf{S}_k \sim \text{Wishart}(\nu_k, \alpha_k \boldsymbol{\Sigma}_{\mathcal{O}_k})$$

The EM algorithm iterates between:
- **E-step**: Compute conditional expectations of missing blocks
- **M-step**: Update covariance (and optionally scale factors)

### Factor-Analytic Model

The FA model constrains the covariance as:

$$\boldsymbol{\Sigma} = \boldsymbol{\Lambda}\boldsymbol{\Lambda}^\top + \boldsymbol{\Psi}$$

where:
- **Λ** (p × k) = factor loadings (population structure axes)
- **Ψ** (diagonal) = individual-specific unique variances
- **k** = number of latent factors (population subgroups)

The M-step uses **Rubin & Thayer (1982)** inner iterations to update Λ and Ψ, followed by varimax rotation for interpretability. This low-rank structure enables:
1. **Prediction of unobserved pairs** via the inner product Λ[i,·] · Λ[j,·]
2. **Regularization** through the rank-k + diagonal constraint
3. **Biological interpretability** of factors as population structure axes

See `vignette("statistical-methods")` and `vignette("combining-grms-factor-model")` for complete mathematical details.

## Performance

**Computational complexity:** O(K × p³) per iteration

**Typical convergence:** 20-100 iterations (FA model); 50-200 (free model)

**For large GRMs (p > 100):**
- Use factor model (`n_factors = k`) — fewer parameters, faster convergence
- Start with `se_method = "none"` to assess convergence
- Set `n_factors` manually rather than `"auto"`
- Consider `se_method = "sem"` (validate against bootstrap)

## Troubleshooting

**Convergence issues?**
- Increase `max_iter` or relax `tol`
- Try `init_sigma = "identity"`
- Increase ridge: `control = list(ridge = 1e-6)`

**Non-positive definite estimates?**
- Check overlap: all individuals connected through shared anchors?
- Increase ridge regularization
- Use `scale_method = "estimate"` if heterogeneous scaling
- Try factor model — it's regularized by construction

**Factor model: how many factors?**
- Start with k = number of known population subgroups
- Compare 2–3 values by BIC
- Use `n_factors = "auto"` for small problems only

See `vignette("study-design-guide")` for comprehensive diagnostics.

## Contributing

Contributions are welcome! Please submit issues or pull requests at:
https://github.com/dakdemir-nmdp/CovCombR

## License

MIT License. See [LICENSE](LICENSE) for details.

## Citation

If you use CovCombR in your research, please cite:

```
Akdemir, D. (2025). CovCombR: Combining Incomplete Covariance Matrices
with Factor-Analytic EM. R package version 1.6.0.
https://github.com/dakdemir-nmdp/CovCombR
```

## References

- Rubin, D. B. & Thayer, D. T. (1982). EM algorithms for ML factor analysis. *Psychometrika*, 47(1), 69–76.
- VanRaden, P. M. (2008). Efficient methods to compute genomic predictions. *Journal of Dairy Science*, 91(11), 4414–4423.
- Perez, P. & de los Campos, G. (2014). Genome-wide regression and prediction with the BGLR statistical package. *Genetics*, 198(2), 483–495.

## Contact

- **Author:** Deniz Akdemir
- **Email:** dakdemir@nmdp.org
- **GitHub:** https://github.com/dakdemir-nmdp/CovCombR
- **Issues:** https://github.com/dakdemir-nmdp/CovCombR/issues

---

**Quick Links:**
- [Installation](#installation)
- [Quick Start](#quick-start-combining-grms-with-a-factor-model)
- [Factor-Analytic vs. Free Models](#two-model-types)
- [Documentation](#documentation)
- [Vignettes](https://dakdemir-nmdp.github.io/CovCombR/)
