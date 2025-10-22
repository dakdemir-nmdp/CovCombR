# CovCombR: Combining Incomplete Covariance Matrices with EM

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-1.5.0-blue)](https://github.com/dakdemir-nmdp/CovCombR)

## Overview

**CovCombR** implements the Expectation-Maximization (EM) algorithm to combine incomplete covariance (or genetic relationship) matrices into a single data-scale covariance estimate. Internally the method leverages Wishart theory, but the primary deliverable is the combined covariance that respects the original measurement scale. The package is designed for situations where:

- Multiple experiments measure overlapping but distinct variable sets
- Missing data patterns prevent direct pooling of covariance estimates
- Heterogeneous scaling across experiments requires correction

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

# Simulate incomplete covariance samples (Wishart draws)
set.seed(2025)
p <- 8
var_names <- paste0("Var", seq_len(p))
true_Sigma <- diag(p)
true_Sigma[1:4, 1:4] <- 0.4
true_Sigma[5:8, 5:8] <- 0.35
true_Sigma[1:4, 5:8] <- 0.1
true_Sigma[5:8, 1:4] <- 0.1
diag(true_Sigma) <- 1
dimnames(true_Sigma) <- list(var_names, var_names)

# Three samples with different observation patterns
nu <- c(Sample_A = 60, Sample_B = 70, Sample_C = 55)
index_list <- list(
  Sample_A = 1:5,
  Sample_B = 3:8,
  Sample_C = c(1, 3, 5, 7)
)
S_list <- lapply(names(index_list), function(name) {
  idx <- index_list[[name]]
  S <- rWishart(1, nu[name], true_Sigma[idx, idx])[,,1] / nu[name]
  dimnames(S) <- list(var_names[idx], var_names[idx])
  S
})
names(S_list) <- names(index_list)

# Fit model (no scaling by default: scale_method = "none")
fit <- fit_covcomb(
  S_list,
  nu,
  scale_method = "none",
  se_method = "plugin",
  control = list(max_iter = 1500, tol = 1e-5, ridge = 1e-4)
)

# Combined covariance on the data scale
S_hat <- fitted(fit)

# Parameter-scale matrix (useful for inference/SEs)
Sigma_hat <- coef(fit)

# Standard errors on both scales
Sigma_se <- fit$Sigma_se
S_hat_se <- fit$S_hat_se
```

## Key Features

### Robust EM Algorithm
- Guaranteed monotonic increase in likelihood
- Automatic convergence detection
- Multiple initialization strategies

### Missing Data Handling
- Arbitrary missing data patterns
- No restrictions on overlap structure
- Efficient conditional expectations

### Covariance-Focused Output
- Returns data-scale covariance (`fitted()`) ready for downstream use
- Provides parameter-scale `Sigma_hat` for standard errors and diagnostics
- Optional experiment-specific scale estimation when required

### Standard Errors
- Plugin method: Fast but anti-conservative
- Bootstrap method: Parametric simulation from the fitted model
- Outputs standard errors on both parameter and data scales

### Numerical Stability
- Cholesky decomposition for inversions
- Positive definiteness projection
- Ridge regularization options
- Memoization for performance

## User Interface

### Main Function: `fit_covcomb()`

```r
fit <- fit_covcomb(
  S_list,                  # List of covariance/Wishart matrices
  nu,                      # Degrees of freedom vector
  init_sigma = "eigendecomp",
  se_method = "plugin",
  control = list(
    max_iter = 1000,
    tol = 1e-7,
    ridge = 1e-8
  )
)
```

#### Input Arguments

**`S_list`** (required, list of matrices)
- List of sample covariance matrices, one per experiment/study
- Each matrix must be square, symmetric, and positive semi-definite
- Matrices can have different dimensions (incomplete data)
- Names optional but recommended for identification
- Example: `list(Study1 = S1, Study2 = S2, Study3 = S3)`

**`nu`** (required, numeric vector)
- Degrees of freedom for each sample
- Must equal the length of `S_list`
- Typically: number of observations minus 1 (for sample covariance)
- Or: effective number of markers (for genomic relationships)
- Can be named to match `S_list` entries
- Example: `c(Study1 = 100, Study2 = 95, Study3 = 110)`

**`scale_method`** (optional, character, default: `"none"`)
- `"none"`: No scaling; matrices assumed to be on comparable scales
  - Use when input matrices are already standardized
  - Fastest option, recommended default
- `"estimate"`: Jointly estimate scale factors $\alpha_k$ for each sample
  - When input matrices may have heterogeneous scaling
  - Returns `fit$alpha_hat` with estimated scale factors
  - Useful for combining GRMs from different marker panels

**`alpha_normalization`** (optional, character, default: `"geometric"`)
- Only applies when `scale_method = "estimate"`
- `"geometric"`: Constrains geometric mean of $\alpha_k$ to 1
  - Preferred for multiplicative scale effects
- `"arithmetic"`: Constrains arithmetic mean of $\alpha_k$ to 1
  - Alternative normalization scheme

**`init_sigma`** (optional, character, default: `"eigendecomp"`)
- `"eigendecomp"`: Initialize from eigendecomposition of stitched matrix
  - Most stable, recommended for most applications
- `"identity"`: Start with identity matrix
  - Faster but may require more iterations
- Numeric matrix: Provide custom initialization

**`se_method`** (optional, character, default: `"plugin"`)
- `"none"`: No standard error computation (fastest)
- `"plugin"`: Fisher information-based standard errors (fast, may underestimate)
- `"bootstrap"`: Parametric bootstrap standard errors (slower, more reliable)
  - Additional arguments: `boot_samples` (default: 100)

**`control`** (optional, list of parameters)
- `max_iter`: Maximum EM iterations (default: 1000)
- `tol`: Convergence tolerance for relative change (default: 1e-7)
- `ridge`: Ridge regularization for Cholesky stability (default: 1e-8)
- `min_eigen`: Minimum eigenvalue floor (default: 1e-10)

#### Return Value: Fitted Object

A list of class `"covcomb"` containing:

**Primary Covariance Estimate:**
- `$Sigma_hat`: Combined covariance matrix (dimension: $p \times p$)
  - The main output: estimated covariance on the same scale as input matrices
  - Use with `coef(fit)` for extraction
  - Ready for downstream analysis, genomic prediction, etc.
- `$Sigma_se`: Standard errors for $\Sigma_{hat}$ (dimension: $p \times p$)
  - Only when `se_method != "none"`
  - Plugin SEs may underestimate uncertainty for imputed entries
  - Bootstrap SEs recommended for formal inference

**Data-Scale Presentation (when `scale_method = "estimate"`):**
- $\hat{S}$: Rescaled version of $\hat{\Sigma}$ (dimension: $p \times p$)
  - When `scale_method = "none"` or `alpha_normalization = "geometric"`: equals $\hat{\Sigma}$
  - When `scale_method = "estimate"` with `alpha_normalization = "arithmetic"`: equals $c \times \hat{\Sigma}$ where $c$ is the scale factor
  - Use with `fitted(fit)`
- Scale factor $c$: Rescaling factor (scalar)
  - Equals 1 when `scale_method = "none"` or `alpha_normalization = "geometric"`
  - Equals mean($\alpha_k$) when `scale_method = "estimate"` with `alpha_normalization = "arithmetic"`
- $\text{SE}(\hat{S})$: Standard errors for $\hat{S}$ (dimension: $p \times p$)
  - Only when `se_method != "none"`

**Scale Estimation Results:**
- `$alpha_hat`: Vector of estimated scale factors (when `scale_method = "estimate"`)
  - Length: number of samples in `S_list`
  - Named if `S_list` was named
  - Interpretation: relative scaling of each matrix

**Convergence Diagnostics:**
- `$convergence$converged`: Logical, TRUE if converged
- `$convergence$iterations`: Number of EM iterations performed
- `$convergence$rel_change`: Final relative change in objective

**Algorithm Details:**
- `$log_lik`: Final log-likelihood value
- `$call`: The function call

### S3 Methods

| Method | Description | Output |
|--------|-------------|--------|
| `print(fit)` | Summary of fit with covariance stats | Console output |
| `summary(fit)` | Detailed diagnostics & eigenspectra | Console output |
| `coef(fit)` | Extract combined covariance $\hat{\Sigma}$ | $p \times p$ matrix |
| `fitted(fit)` | Extract rescaled presentation $\hat{S}$ | $p \times p$ matrix |

#### Using the Fitted Object

```r
# Extract main results
Sigma_hat <- coef(fit)            # Primary output: combined covariance estimate
S_hat <- fitted(fit)              # Rescaled version (often identical to Sigma_hat)

# Standard errors (if computed)
Sigma_se <- fit$Sigma_se          # SE for Sigma_hat
S_se <- fit$S_hat_se              # SE for S_hat (rescaled)

# Check convergence
print(fit$convergence)

# Quick summary
print(fit)                         # Compact output
summary(fit)                       # Detailed diagnostics

# Downstream use for genomic prediction
combined_cov <- coef(fit)  # or fitted(fit), often identical
# Use combined_cov in GWAS or prediction models
```

### Main Functions Summary

| Function | Purpose | Returns |
|----------|---------|---------|
| `fit_covcomb()` | Main EM fitting | List with class `"covcomb"` |
| `print()` | Display fit summary | Console output |
| `summary()` | Detailed diagnostics | Console output + eigen spectra |
| `coef()` | Extract combined covariance $\hat{\Sigma}$ | $p \times p$ matrix |
| `fitted()` | Extract rescaled covariance $\hat{S}$ | $p \times p$ matrix |

## Documentation

Comprehensive documentation is available:

```r
# Main function help
?fit_covcomb

# Package overview
?CovCombR

# Full tutorial vignette
vignette("CovCombR-vignette")

# Iris example
vignette("iris-example")

# GRM combination vignette
vignette("combining-grms")

# Statistical methods vignette
vignette("statistical-methods")
```

## Interpreting Results

### The Two Output Scales

CovCombR provides estimates on **two related scales**:

**Primary Output (`coef(fit)` / $\hat{\Sigma}$):**
- The combined covariance estimate on the same scale as your input matrices
- Appropriate for genomic prediction, GWAS, downstream analysis
- Use this for most applications

**Data-Scale Presentation (`fitted(fit)` / $\hat{S}$):**
- When `scale_method = "none"`: identical to $\hat{\Sigma}$
- When `scale_method = "estimate"` with `alpha_normalization = "geometric"`: identical to $\hat{\Sigma}$
- When `scale_method = "estimate"` with `alpha_normalization = "arithmetic"`: rescaled version $\hat{S} = c \times \hat{\Sigma}$ where $c$ is the scale factor
- Primarily for compatibility; use $\hat{\Sigma}$ for most applications

### Understanding Scale Parameters

**Experiment-Specific Scales $\hat{\alpha}_k$** (when `scale_method = "estimate"`):
- Relative scaling factors for each experiment
- Stored in `fit$alpha_hat`
- Normalized so geometric mean = 1 (default) or arithmetic mean = 1
- Useful for diagnostics when inputs differ in scale
- When `scale_method = "none"`: all equal to 1

### Standard Errors

**Plugin Standard Errors** (`se_method = "plugin"`):
- Fast to compute (no simulation needed)
- Based on observed Fisher information
- Anti-conservative (tend to be too small)
- Good for quick diagnostics

**Bootstrap Standard Errors** (`se_method = "bootstrap"`):
- Parametric bootstrap from fitted model
- More realistic uncertainty quantification
- Recommended for publication and inference
- Slower but more reliable (especially for small samples)

Access results:
```r
Sigma_se <- fit$Sigma_se        # SE on parameter scale
S_hat_se <- fit$S_hat_se        # SE on data scale

# Compute 95% confidence intervals
lower_95 <- coef(fit) - 1.96 * fit$Sigma_se
upper_95 <- coef(fit) + 1.96 * fit$Sigma_se
```

## Example Workflows

The following helper constructs a reproducible synthetic data set used across
the next examples. Run it once to make the examples fully self-contained.

```r
simulate_multi_site_example <- function(seed = 456) {
  if (!is.null(seed)) set.seed(seed)

  p <- 12
  var_names <- paste0("Var", seq_len(p))

  true_cov <- diag(p)
  true_cov[1:6, 1:6] <- 0.35
  true_cov[7:12, 7:12] <- 0.25
  true_cov[1:6, 7:12] <- 0.05
  true_cov[7:12, 1:6] <- 0.05
  diag(true_cov) <- 1
  dimnames(true_cov) <- list(var_names, var_names)

  sites <- list(
    Site_A = 1:8,
    Site_B = 5:12,
    Site_C = c(1:4, 9:12)
  )
  nu <- c(Site_A = 260, Site_B = 200, Site_C = 220)

  S_list <- lapply(names(sites), function(site) {
    idx <- sites[[site]]
    S <- rWishart(1, nu[site], true_cov[idx, idx])[,,1] / nu[site]
    dimnames(S) <- list(var_names[idx], var_names[idx])
    S
  })
  names(S_list) <- names(sites)

  list(S_list = S_list, nu = nu, true_cov = true_cov)
}
```

### Example 1: Multi-Site Study (Data Scale Focus)

```r
example <- simulate_multi_site_example(seed = 456)
S_list <- example$S_list
nu <- example$nu

fit <- fit_covcomb(
  S_list,
  nu,
  scale_method = "none",
  se_method = "plugin",
  control = list(max_iter = 1500, tol = 1e-5, ridge = 1e-4)
)

combined_cov <- coef(fit)
print(fit)

# Visualize the combined covariance matrix
heatmap(combined_cov)
```

### Example 2: Multi-GRM in Genomics (Scale Estimation)

```r
example <- simulate_multi_site_example(seed = 456)
scale_factors <- c(Site_A = 1.0, Site_B = 1.15, Site_C = 0.9)

scaled_inputs <- Map(function(S, alpha) alpha * S, example$S_list, scale_factors)

fit <- fit_covcomb(
  S_list = scaled_inputs,
  nu = example$nu,
  scale_method = "estimate",
  se_method = "none",
  control = list(max_iter = 1500, tol = 1e-5, ridge = 1e-4)
)

print(fit$alpha_hat)
genomic_relationship <- coef(fit)
```

### Example 3: Statistical Inference (Parameter Scale Focus)

```r
example <- simulate_multi_site_example(seed = 456)

fit <- fit_covcomb(
  S_list = example$S_list,
  nu = example$nu,
  se_method = "bootstrap",
  control = list(
    max_iter = 1500,
    tol = 1e-5,
    ridge = 1e-4,
    bootstrap = list(B = 8, seed = 99, progress = FALSE)
  )
)

Sigma_est <- coef(fit)
Sigma_se <- fit$Sigma_se

idx1 <- which(rownames(Sigma_est) == "Var1")
idx2 <- which(rownames(Sigma_est) == "Var2")
corr_se <- Sigma_se[idx1, idx2]
z_stat <- Sigma_est[idx1, idx2] / corr_se
p_value <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)

cat("Correlation:", round(Sigma_est[idx1, idx2], 3), "\n")
cat("SE:", round(corr_se, 3), "\n")
cat("p-value:", signif(p_value, 3), "\n")
```

### Example 4: Multi-Site Study Simulation

```r
set.seed(789)
n_genes <- 40
gene_names <- paste0("Gene", seq_len(n_genes))
true_cov <- matrix(0.5, n_genes, n_genes)
diag(true_cov) <- 1
dimnames(true_cov) <- list(gene_names, gene_names)

site_genes <- list(
  Site1 = 1:28,
  Site2 = 10:35,
  Site3 = 18:40,
  Site4 = c(1:12, 24:36)
)
nu_sites <- c(110, 130, 120, 100)

S_sites <- vector("list", length(site_genes))
for (i in seq_along(site_genes)) {
  genes <- site_genes[[i]]
  S <- rWishart(1, nu_sites[i], true_cov[genes, genes])[,,1] / nu_sites[i]
  dimnames(S) <- list(gene_names[genes], gene_names[genes])
  S_sites[[i]] <- S
}
names(S_sites) <- names(site_genes)
nu_sites <- setNames(nu_sites, names(site_genes))

fit <- fit_covcomb(
  S_list = S_sites,
  nu = nu_sites,
  scale_method = "estimate",
  control = list(max_iter = 1500, tol = 1e-5, ridge = 1e-4)
)

fro_error <- norm(coef(fit) - true_cov, "F") / norm(true_cov, "F")
print(paste("Relative error:", signif(fro_error, 3)))
```

## Statistical Background

### Model

Given $K$ samples where sample $k$ observes variables $\mathcal{O}_k$:

$$
\mathbf{S}_k = \mathbf{W}_{k,\mathcal{O}_k\mathcal{O}_k} \sim \text{Wishart}(\nu_k, \alpha_k \boldsymbol{\Sigma}_{\mathcal{O}_k\mathcal{O}_k})
$$

### EM Algorithm

**E-step:** Compute conditional expectations using Wishart properties

**M-step:** Update covariance and (optionally) scale factors:

$$
\boldsymbol{\Sigma}^{(t+1)} = \frac{1}{\sum_k \nu_k} \sum_{k=1}^K \frac{1}{\alpha_k} E[\mathbf{W}_k | \mathbf{S}_k, \boldsymbol{\Sigma}^{(t)}]
$$

### References

- McLachlan, G. J., & Krishnan, T. (2008). *The EM Algorithm and Extensions* (2nd ed.). Wiley.
- Anderson, T. W. (2003). *An Introduction to Multivariate Statistical Analysis* (3rd ed.). Wiley.

## Performance

Computational complexity per iteration: $O(Kp^3)$

### Tips for Large Problems

1. Start with `se_method = "none"` to assess convergence
2. Use `init_sigma = "identity"` for faster initialization
3. Implement parallel bootstrap for standard errors

## Advanced Usage

### Bootstrap Standard Errors

CovCombR now includes built-in parametric bootstrap for more reliable standard errors:

```r
# Built-in bootstrap (recommended for publication-quality inference)
example <- simulate_multi_site_example(seed = 888)

fit_boot <- fit_covcomb(
  example$S_list,
  example$nu,
  se_method = "bootstrap",
  control = list(
    max_iter = 1500,
    tol = 1e-5,
    ridge = 1e-4,
    bootstrap = list(
      B = 200,              # Number of bootstrap replicates
      seed = 2025,          # For reproducibility
      progress = TRUE,      # Show progress bar
      verbose = FALSE       # Suppress replicate warnings
    )
  )
)

# Access bootstrap results
Sigma_se <- fit_boot$Sigma_se          # SEs on parameter scale
S_hat_se <- fit_boot$S_hat_se          # SEs on data scale
boot_meta <- fit_boot$bootstrap        # Bootstrap metadata (successes, failures, etc.)
```

The bootstrap accounts for both sampling variance and EM imputation uncertainty, providing valid inference even for entries never jointly observed in any study.

### Custom Control Parameters

```r
fit <- fit_covcomb(
  S_list, nu,
  control = list(
    max_iter = 1000,    # Maximum iterations
    tol = 1e-8,         # Convergence tolerance
    ridge = 1e-8,       # Ridge for stability
    min_eigen = 1e-10   # Eigenvalue floor
  )
)
```

## Troubleshooting

### Issue: Poor Convergence or Non-Convergence

**Symptoms:** Algorithm stops before iteration limit with `convergence$converged = FALSE`

**Solutions:**
```r
# 1. Relax convergence tolerance
fit <- fit_covcomb(S_list, nu, control = list(tol = 1e-6))

# 2. Increase maximum iterations
fit <- fit_covcomb(S_list, nu, control = list(max_iter = 5000))

# 3. Use simpler initialization
fit <- fit_covcomb(S_list, nu, init_sigma = "identity")

# 4. Increase ridge regularization for stability
fit <- fit_covcomb(S_list, nu, control = list(ridge = 1e-6))
```

### Issue: Non-Positive Definite Estimates

**Symptoms:** Error about Cholesky decomposition failure or negative eigenvalues

**Solutions:**
```r
# 1. Increase ridge regularization
fit <- fit_covcomb(S_list, nu, control = list(ridge = 1e-6))

# 2. Project to PSD after fitting
S_fit <- fitted(fit)
eig <- eigen(S_fit, symmetric = TRUE)
S_psd <- eig$vectors %*% diag(pmax(eig$values, 0)) %*% t(eig$vectors)

# 3. Use pre-standardization to stabilize
```

### Issue: Unexpected Scale Factors

**Symptoms:** `$alpha_hat` has unexpected values (very large/small)

**Solutions:**
```r
# 1. Pre-scale inputs to reasonable ranges
S_list_scaled <- lapply(S_list, function(S) S / mean(diag(S)))
fit <- fit_covcomb(S_list_scaled, nu, scale_method = "none")

# 2. Use scale estimation to handle heterogeneity
fit <- fit_covcomb(S_list, nu, scale_method = "estimate")

# 3. Provide external scaling if you have domain knowledge
fit <- fit_covcomb(S_list, nu, scale_method = "none")
```

### Issue: Standard Error Computation Fails

**Symptoms:** Error during `se_method = "bootstrap"` or `"plugin"`

**Solutions:**
```r
# 1. Start without SEs to verify fit
fit <- fit_covcomb(S_list, nu, se_method = "none")

# 2. Use plugin method (faster, usually works)
fit <- fit_covcomb(S_list, nu, se_method = "plugin")

# 3. Reduce bootstrap samples if memory-limited
fit <- fit_covcomb(S_list, nu, se_method = "bootstrap", 
                      boot_samples = 50)

# 4. Increase ridge regularization
fit <- fit_covcomb(S_list, nu, se_method = "bootstrap",
                      control = list(ridge = 1e-6))
```

## Performance Tips

### For Large Problems

1. **Reduce SE computation:**
   ```r
   # Assess convergence first without SEs
   fit <- fit_covcomb(S_list, nu, se_method = "none")
   
   # Then compute SEs if needed
   fit <- fit_covcomb(S_list, nu, se_method = "plugin")
   ```

2. **Pre-scale inputs to unit diagonal:**
   ```r
   # Normalize each matrix to have unit diagonal
   S_list_scaled <- lapply(S_list, function(S) {
     cov2cor(S)  # Convert to correlation matrix
   })
   fit <- fit_covcomb(S_list_scaled, nu, scale_method = "none")
   ```

3. **Coarsen convergence tolerance:**
   ```r
   # Sufficient for many applications
   fit <- fit_covcomb(S_list, nu, 
                        control = list(tol = 1e-5, max_iter = 500))
   ```

### Memory Management

- For very large $p$ (e.g., $p > 10000$):
  - Consider block-wise analysis
  - Use single precision matrices if applicable
  - Pre-allocate bootstrap arrays

### Computational Complexity

- **Per-iteration cost:** $O(Kp^3)$ dominated by matrix inversions
- **Total cost:** $\approx$ iterations $\times K \times p^3$
- Typical convergence: 20-50 iterations for well-conditioned problems
- With bootstrap: Multiply by `boot_samples` factor

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This package is licensed under the MIT License. See [LICENSE](LICENSE) file for details.

## Citation

If you use this package in your research, please cite:

```
Akdemir, D. (2025). CovCombR: Combining Incomplete Covariance Matrices with EM.
R package version 1.5.0.
```

## Contact

- **Author:** Deniz Akdemir
- **GitHub:** https://github.com/dakdemir-nmdp/CovCombR
- **Issues:** https://github.com/dakdemir-nmdp/CovCombR/issues

## Acknowledgments

This package implements statistical methods described in McLachlan & Krishnan (2008) with enhancements for numerical stability and missing data handling.
