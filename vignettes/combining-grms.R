## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width = 8,
  fig.height = 6,
  fig.align = "center"
)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
# Suppress lint warnings for ggplot aesthetic variables
utils::globalVariables(c("Col", "Row", "value_plot", "is_observed"))
library(CovCombR)

## ----simulation-setup---------------------------------------------------------
# Simulation parameters
n_total <- 1000 # Total genotypes
n_markers_study <- 5000 # Total markers per study
min_genotypes_per_study <- 200 # Minimum genotypes sampled per study
max_genotypes_per_study <- 400 # Maximum genotypes sampled per study
n_markers_per_study <- n_markers_study # Alias for inline code
n_studies <- 10

# We'll generate a smaller set of base markers to save time/code for the vignette
n_markers_base_sim <- 10000
n_markers_total <- n_markers_base_sim # Alias for total markers used in inline code

# Generate base population structure (3 subpopulations)
subpop_sizes <- c(400, 350, 250) # Scaled up to sum to 1000
subpop_labels <- rep(1:3, times = subpop_sizes)
genotype_ids <- paste0("G", sprintf("%03d", 1:n_total))

# Generate true correlation matrix with block structure
true_cor <- matrix(0.1, n_total, n_total)
for (i in 1:3) {
  idx <- which(subpop_labels == i)
  true_cor[idx, idx] <- 0.6
}
diag(true_cor) <- 1
# Ensure PD
eig <- eigen(true_cor, symmetric = TRUE)
true_cor <- eig$vectors %*% diag(pmax(eig$values, 0.01)) %*% t(eig$vectors)
true_cor <- cov2cor(true_cor)
dimnames(true_cor) <- list(genotype_ids, genotype_ids)

# Generate marker data consistent with this structure
L <- chol(true_cor)
marker_mat <- matrix(0, n_total, n_markers_base_sim)
marker_freqs <- runif(n_markers_base_sim, 0.1, 0.9)

# Helper for VanRaden GRM
calc_grm <- function(M, ids = NULL) {
  p <- colMeans(M) / 2
  Z <- M - matrix(2 * p, nrow(M), ncol(M), byrow = TRUE)
  G <- tcrossprod(Z) / (sum(diag(tcrossprod(Z))) / nrow(M))
  if (!is.null(ids)) dimnames(G) <- list(ids, ids)
  G
}

# Generate markers (simplified for vignette)
# Using a subset of markers for speed in this example
for (m in 1:n_markers_base_sim) {
  z <- crossprod(L, rnorm(n_total))
  p <- marker_freqs[m]
  thresh1 <- qnorm((1 - p)^2)
  thresh2 <- qnorm((1 - p)^2 + 2 * p * (1 - p))
  marker_mat[, m] <- ifelse(z < thresh1, 0, ifelse(z < thresh2, 1, 2))
}
rownames(marker_mat) <- genotype_ids
colnames(marker_mat) <- paste0("M", 1:n_markers_base_sim)

true_grm <- calc_grm(marker_mat, genotype_ids)

# Generate Study GRMs with overlap
study_grms <- list()
study_genotypes <- list()
study_markers <- list()

set.seed(2025)
for (k in 1:n_studies) {
  # Sample genotypes with overlap
  # Higher probability for already sampled individuals to ensure connectivity
  if (k == 1) {
    prob <- rep(1, n_total)
  } else {
    sampled <- unique(unlist(study_genotypes))
    prob <- rep(1, n_total)
    prob[which(genotype_ids %in% genotype_ids[sampled])] <- 3 # 3x more likely to be resampled
  }

  n_s <- sample(200:400, 1)
  s_idx <- sort(sample(1:n_total, n_s, prob = prob))
  study_genotypes[[k]] <- s_idx

  # Sample markers
  m_idx <- sort(sample(1:n_markers_base_sim, n_markers_study))
  study_markers[[k]] <- m_idx

  # Calculate GRM
  M_sub <- marker_mat[s_idx, m_idx]
  study_grms[[k]] <- calc_grm(M_sub, genotype_ids[s_idx])
}
names(study_grms) <- paste0("Study", 1:n_studies)

# Track marker counts for each study (needed for degrees of freedom)
study_marker_counts <- sapply(study_markers, length)

## ----study-overlap------------------------------------------------------------
# Calculate pairwise overlap between studies (genotypes)
genotype_overlap_matrix <- matrix(0, n_studies, n_studies)
for (i in 1:n_studies) {
  for (j in 1:n_studies) {
    genotype_overlap_matrix[i, j] <- length(intersect(
      study_genotypes[[i]],
      study_genotypes[[j]]
    ))
  }
}

colnames(genotype_overlap_matrix) <- paste0("S", 1:n_studies)
rownames(genotype_overlap_matrix) <- paste0("S", 1:n_studies)

cat("\nGenotype overlap between studies:\n")
print(genotype_overlap_matrix)

# Calculate pairwise marker overlap between studies
marker_overlap_matrix <- matrix(0, n_studies, n_studies)
for (i in 1:n_studies) {
  for (j in 1:n_studies) {
    marker_overlap_matrix[i, j] <- length(intersect(
      study_markers[[i]],
      study_markers[[j]]
    ))
  }
}

colnames(marker_overlap_matrix) <- paste0("S", 1:n_studies)
rownames(marker_overlap_matrix) <- paste0("S", 1:n_studies)

cat("\nMarker overlap between studies:\n")
print(marker_overlap_matrix)

cat("\nStudy coverage statistics:\n")
all_genotypes_covered <- unique(unlist(study_genotypes))
cat(sprintf(
  "  Genotypes: %d / %d covered (%.1f%%)\n",
  length(all_genotypes_covered), n_total,
  100 * length(all_genotypes_covered) / n_total
))

genotype_frequency <- table(unlist(study_genotypes))
cat(sprintf(
  "  Mean studies per genotype: %.1f (range: %d-%d)\n",
  mean(genotype_frequency),
  min(genotype_frequency), max(genotype_frequency)
))

all_markers_used <- unique(unlist(study_markers))
cat(sprintf(
  "  Markers: %d / %d used (%.1f%%)\n",
  length(all_markers_used), n_markers_base_sim,
  100 * length(all_markers_used) / n_markers_base_sim
))

marker_frequency <- table(unlist(study_markers))
cat(sprintf(
  "  Mean studies per marker: %.1f (range: %d-%d)\n",
  mean(marker_frequency),
  min(marker_frequency), max(marker_frequency)
))

## ----compare-relationships, fig.width=7, fig.height=5, out.width="90%"--------
# Find genotypes that appear in at least 3 studies
# Convert study_genotypes to genotype IDs
all_study_genotype_ids <- lapply(study_genotypes, function(idx) genotype_ids[idx])
genotype_id_counts <- table(unlist(all_study_genotype_ids))
well_sampled_genotypes <- names(genotype_id_counts)[genotype_id_counts >= 3]

if (length(well_sampled_genotypes) >= 10) {
  # Select 10 well-sampled genotypes
  selected_genotypes <- sample(
    well_sampled_genotypes,
    min(10, length(well_sampled_genotypes))
  )

  # For each pair, collect estimates from different studies
  comparison_data <- list()

  for (i in 1:(length(selected_genotypes) - 1)) {
    for (j in (i + 1):length(selected_genotypes)) {
      g1 <- selected_genotypes[i]
      g2 <- selected_genotypes[j]

      true_value <- true_grm[g1, g2]

      # Find studies that have both genotypes
      studies_with_pair <- which(sapply(all_study_genotype_ids, function(study_ids) {
        (g1 %in% study_ids) && (g2 %in% study_ids)
      }))

      if (length(studies_with_pair) > 0) {
        for (study in studies_with_pair) {
          study_value <- study_grms[[study]][g1, g2]
          comparison_data[[length(comparison_data) + 1]] <- data.frame(
            Genotype1 = g1,
            Genotype2 = g2,
            Study = study,
            True_Relationship = true_value,
            Estimated_Relationship = study_value,
            Error = study_value - true_value
          )
        }
      }
    }
  }

  comparison_df <- do.call(rbind, comparison_data)

  # Plot: True vs Estimated
  ggplot(comparison_df, aes(x = True_Relationship, y = Estimated_Relationship)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "Sample GRM vs True GRM",
      subtitle = "Relationships between well-sampled genotype pairs",
      x = "True Relationship (Population GRM)",
      y = "Estimated Relationship (Study GRM)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  cat("\nComparison statistics:\n")
  cat(sprintf(
    "  Number of genotype pairs compared: %d\n",
    nrow(comparison_df)
  ))
  cat(sprintf("  RMSE: %.4f\n", sqrt(mean(comparison_df$Error^2))))
  cat(sprintf("  MAE: %.4f\n", mean(abs(comparison_df$Error))))
  cat(sprintf(
    "  Correlation: %.4f\n",
    cor(
      comparison_df$True_Relationship,
      comparison_df$Estimated_Relationship
    )
  ))
}

## ----aligned-grm-viz, fig.width=9, fig.height=18, out.width="100%", fig.alt="Aligned heatmaps showing observed and missing entries for each study GRM alongside the true GRM."----
# Create full population matrices for each study, with NA for unobserved pairs
study_grms_full <- lapply(1:n_studies, function(k) {
  # Create a matrix for the full population
  full_grm <- matrix(NA, n_total, n_total)
  rownames(full_grm) <- genotype_ids
  colnames(full_grm) <- genotype_ids

  # Fill in observed values
  study_idx <- study_genotypes[[k]]
  study_ids <- genotype_ids[study_idx]
  full_grm[study_ids, study_ids] <- study_grms[[k]]

  return(full_grm)
})

# Order genotypes by hierarchical clustering of the true GRM
# This will group similar genotypes together
hc <- hclust(as.dist(1 - true_grm), method = "average")
ordered_genotype_ids <- genotype_ids[hc$order]

# Helper function to plot aligned GRM heatmap
plot_aligned_grm <- function(grm_matrix, title = "GRM", ordered_ids,
                             show_observed_only = FALSE,
                             text_size = 8) {
  n <- nrow(grm_matrix)

  # Reorder matrix
  grm_ordered <- grm_matrix[ordered_ids, ordered_ids]

  # Convert to long format for ggplot
  grm_long <- expand.grid(
    Row = 1:n,
    Col = 1:n
  )
  grm_long$value <- as.vector(grm_ordered)
  grm_long$is_observed <- !is.na(grm_long$value)

  # For visualization, set NA to a specific value (will be colored gray)
  grm_long$value_plot <- ifelse(is.na(grm_long$value), -999, grm_long$value)

  # Create plot
  p <- ggplot(grm_long, aes(x = Col, y = n - Row + 1)) +
    geom_tile(aes(fill = value_plot, alpha = is_observed)) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0.5,
      limits = c(-0.5, 2),
      name = "Relationship",
      na.value = "#808080"
    ) +
    scale_alpha_manual(
      values = c("TRUE" = 1, "FALSE" = 0.3),
      guide = "none"
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = text_size) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.key.height = unit(0.3, "cm"),
      legend.key.width = unit(0.15, "cm")
    ) +
    coord_equal()

  return(p)
}

# Create aligned heatmaps for all studies
aligned_study_plots <- lapply(1:n_studies, function(k) {
  n_obs <- length(study_genotypes[[k]])
  n_pairs_observed <- n_obs * (n_obs + 1) / 2
  coverage_pct <- 100 * n_pairs_observed / (n_total * (n_total + 1) / 2)

  plot_aligned_grm(
    study_grms_full[[k]],
    title = sprintf("Study %d\n(%.1f%% coverage)", k, coverage_pct),
    ordered_ids = ordered_genotype_ids,
    text_size = 8
  )
})

# Create aligned heatmap for true GRM
true_grm_aligned_plot <- plot_aligned_grm(
  true_grm,
  title = "True Population GRM\n(100% coverage)",
  ordered_ids = ordered_genotype_ids,
  text_size = 8
)

# Arrange all plots in a grid (1 true + n_studies)
all_aligned_plots <- c(list(true_grm_aligned_plot), aligned_study_plots)
n_plots_aligned <- length(all_aligned_plots)
n_cols_aligned <- 3
n_rows_aligned <- ceiling(n_plots_aligned / n_cols_aligned)
grid.arrange(grobs = all_aligned_plots, ncol = n_cols_aligned, nrow = n_rows_aligned)

## ----coverage-by-block, fig.width=9, fig.height=4, out.width="100%", fig.alt="Heatmaps summarizing within- and between-population coverage percentages for every study."----
# Calculate coverage by population structure blocks
# Define blocks based on subpopulation membership (after ordering)
subpop_labels_ordered <- subpop_labels[hc$order]

# For each study, calculate coverage within and between subpopulations
coverage_summary <- do.call(rbind, lapply(1:n_studies, function(k) {
  study_idx <- study_genotypes[[k]]
  study_subpops <- subpop_labels[study_idx]

  # Count coverage for each subpopulation pair
  coverage_data <- expand.grid(Pop1 = 1:3, Pop2 = 1:3)
  coverage_data$Coverage <- sapply(seq_len(nrow(coverage_data)), function(i) {
    pop1 <- coverage_data$Pop1[i]
    pop2 <- coverage_data$Pop2[i]

    # Genotypes in each subpopulation
    pop1_all <- which(subpop_labels == pop1)
    pop2_all <- which(subpop_labels == pop2)

    # Genotypes observed in study
    pop1_obs <- sum(study_idx %in% pop1_all)
    pop2_obs <- sum(study_idx %in% pop2_all)

    # Number of observed pairs
    if (pop1 == pop2) {
      n_pairs_total <- length(pop1_all) * (length(pop1_all) + 1) / 2
      n_pairs_obs <- pop1_obs * (pop1_obs + 1) / 2
    } else {
      n_pairs_total <- length(pop1_all) * length(pop2_all)
      n_pairs_obs <- pop1_obs * pop2_obs
    }

    100 * n_pairs_obs / n_pairs_total
  })

  coverage_data$Study <- paste0("Study ", k)
  return(coverage_data)
}))

# Plot coverage heatmap
ggplot(coverage_summary, aes(x = Pop2, y = Pop1, fill = Coverage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f%%", Coverage)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "#2166AC", name = "Coverage (%)") +
  facet_wrap(~Study, ncol = 5) +
  scale_x_continuous(breaks = 1:3, labels = paste0("Pop", 1:3)) +
  scale_y_continuous(breaks = 1:3, labels = paste0("Pop", 1:3)) +
  labs(
    title = "Coverage of Population Structure Blocks by Study",
    x = "Population 2",
    y = "Population 1"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  ) +
  coord_equal()

## ----prepare-wishart-data-----------------------------------------------------
# For each study, prepare the Wishart-scale matrix and degrees of freedom
# The GRM is already an estimate of the covariance, so we need to scale it
# appropriately and provide degrees of freedom

# Effective degrees of freedom can be approximated by the number of markers
# For Van Raden GRM: nu approximately equals n_markers
# We'll use a slightly conservative estimate

S_list_wishart <- lapply(1:n_studies, function(k) {
  study_grm <- study_grms[[k]]

  # The GRM is already a covariance estimate on the data scale
  # fit_covcomb expects covariance matrices on the data scale (not per-df)
  # So we pass the GRM directly without any scaling
  S_k <- study_grm

  return(S_k)
})

names(S_list_wishart) <- paste0("Study", 1:n_studies)

# Degrees of freedom for each study
nu_vec <- sapply(1:n_studies, function(k) {
  round(study_marker_counts[k] * 0.8)
})
names(nu_vec) <- paste0("Study", 1:n_studies)

cat("Prepared data for CovCombR:\n")
for (k in 1:n_studies) {
  cat(sprintf(
    "  Study %2d: %d genotypes, nu = %d\n",
    k, nrow(study_grms[[k]]), nu_vec[k]
  ))
}

## ----fit-wishart-em, fig.width=7, fig.height=4.5, out.width="90%", fig.alt="Line plot of CovCombR log-likelihood across EM iterations for the simulated GRM example."----
# Fit CovCombR with reduced iterations for testing
fit <- fit_covcomb(
  S_list = S_list_wishart,
  nu = nu_vec,
  init_sigma = "identity",
  scale_method = "none",
  control = list(
    max_iter = 50, # Small number for fast testing
    tol = 1e-5,
    ridge = 1e-6
  ),
  se_method = "none" # Skip SE computation for speed
)

print(fit)

# plot log-likelihood convergence
plot(fit$history$log_likelihood,
  type = "o", xlab = "Iteration", ylab = "Log-Likelihood",
  main = "CovCombR Log-Likelihood Convergence"
)

# Extract the combined GRM
combined_grm_raw <- fit$S_hat

# IMPORTANT: Rescale to match the true GRM scale
# Each study GRM was independently scaled to have mean diagonal = 1
# But the absolute variance scale differs due to different marker counts
# We need to rescale the combined GRM to match the true GRM's scale
#
# Strategy: Use the diagonal elements (which represent self-relationships)
# These should be close to 1.0 in a properly scaled GRM
# We'll rescale so that the mean diagonal matches the true GRM's mean diagonal
true_mean_diag <- mean(diag(true_grm))
combined_mean_diag <- mean(diag(combined_grm_raw))
scale_factor <- true_mean_diag / combined_mean_diag

cat("\nRescaling combined GRM:\n")
cat(sprintf("  True GRM mean diagonal: %.4f\n", true_mean_diag))
cat(sprintf("  Combined GRM mean diagonal (before rescaling): %.4f\n", combined_mean_diag))
cat(sprintf("  Scale factor applied: %.4f\n", scale_factor))

combined_grm_raw <- combined_grm_raw * scale_factor

cat(sprintf("  Combined GRM mean diagonal (after rescaling): %.4f\n", mean(diag(combined_grm_raw))))

# Check which genotypes are in the combined GRM
combined_genotype_ids <- rownames(combined_grm_raw)

cat("\nCombined GRM coverage:\n")
cat(sprintf(
  "  Genotypes in combined GRM: %d / %d\n",
  length(combined_genotype_ids), n_total
))

# Create full population matrix, filling in with NA for unobserved genotypes
combined_grm <- matrix(NA, n_total, n_total)
rownames(combined_grm) <- genotype_ids
colnames(combined_grm) <- genotype_ids

# Fill in observed genotypes
combined_grm[combined_genotype_ids, combined_genotype_ids] <- combined_grm_raw

cat("\nCombined GRM statistics (observed genotypes only):\n")
cat(sprintf("  Dimensions: %d x %d\n", nrow(combined_grm_raw), ncol(combined_grm_raw)))
cat(sprintf("  Mean diagonal: %.3f\n", mean(diag(combined_grm_raw))))
cat(sprintf("  Mean off-diagonal: %.3f\n", mean(combined_grm_raw[upper.tri(combined_grm_raw)])))
cat(sprintf("  SD off-diagonal: %.3f\n", sd(combined_grm_raw[upper.tri(combined_grm_raw)])))

## ----evaluate-recovery, fig.width=9, fig.height=8, out.width="100%", fig.alt="Heatmaps comparing the true GRM, CovCombR combined GRM, and residual differences across genotypes."----
# Create comparison visualization
true_grm_plot <- plot_aligned_grm(
  true_grm,
  title = "True Population GRM",
  ordered_ids = ordered_genotype_ids,
  text_size = 10
)

combined_grm_plot <- plot_aligned_grm(
  combined_grm,
  title = "CovCombR Combined GRM",
  ordered_ids = ordered_genotype_ids,
  text_size = 10
)

# Calculate residuals (only for observed genotypes)
residual_grm <- combined_grm - true_grm

residual_plot <- ggplot(data.frame(
  Row = rep(1:n_total, n_total),
  Col = rep(1:n_total, each = n_total),
  Residual = as.vector(residual_grm[ordered_genotype_ids, ordered_genotype_ids])
), aes(x = Col, y = n_total - Row + 1, fill = Residual)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-0.5, 0.5),
    name = "Residual",
    na.value = "#CCCCCC"
  ) +
  labs(title = "Residuals (Combined - True)", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  coord_equal()

# Scatter plot: True vs Combined with categorical coloring
# First, determine which pairs were observed in at least one study
observed_pairs <- matrix(FALSE, n_total, n_total)
rownames(observed_pairs) <- genotype_ids
colnames(observed_pairs) <- genotype_ids

for (k in 1:n_studies) {
  study_idx <- study_genotypes[[k]]
  study_ids <- genotype_ids[study_idx]
  observed_pairs[study_ids, study_ids] <- TRUE
}

# Create indices for classification
upper_tri_with_diag <- upper.tri(true_grm, diag = TRUE)
diag_idx <- matrix(FALSE, n_total, n_total)
diag(diag_idx) <- TRUE

# Classify each pair
pair_type <- rep(NA_character_, sum(upper_tri_with_diag))
true_vals <- true_grm[upper_tri_with_diag]
combined_vals <- combined_grm[upper_tri_with_diag]
diag_vals <- diag_idx[upper_tri_with_diag]
observed_vals <- observed_pairs[upper_tri_with_diag]

# Classify: diagonal, observed off-diagonal, or unobserved off-diagonal
pair_type[diag_vals] <- "Diagonal"
pair_type[!diag_vals & observed_vals] <- "Observed (off-diagonal)"
pair_type[!diag_vals & !observed_vals] <- "Unobserved (imputed)"

scatter_data <- data.frame(
  True = true_vals,
  Combined = combined_vals,
  Type = factor(pair_type, levels = c("Diagonal", "Observed (off-diagonal)", "Unobserved (imputed)"))
)

# Remove NA values (should only be unobserved genotypes)
scatter_data <- scatter_data[!is.na(scatter_data$Combined), ]

# Print counts for each category
cat("\nScatter plot point counts:\n")
cat(sprintf("  Diagonal: %d\n", sum(scatter_data$Type == "Diagonal", na.rm = TRUE)))
cat(sprintf("  Observed (off-diagonal): %d\n", sum(scatter_data$Type == "Observed (off-diagonal)", na.rm = TRUE)))
cat(sprintf("  Unobserved (imputed): %d\n", sum(scatter_data$Type == "Unobserved (imputed)", na.rm = TRUE)))

# Define colors for each category
type_colors <- c(
  "Diagonal" = "#E31A1C", # Red for diagonal
  "Observed (off-diagonal)" = "#1F78B4", # Blue for observed
  "Unobserved (imputed)" = "#33A02C" # Green for imputed
)

scatter_plot <- ggplot(scatter_data, aes(x = True, y = Combined, color = Type, alpha = Type)) +
  geom_point(size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(
    values = type_colors,
    name = "Pair Type",
    drop = FALSE
  ) +
  scale_alpha_manual(
    values = c("Diagonal" = 0.8, "Observed (off-diagonal)" = 0.6, "Unobserved (imputed)" = 0.4),
    name = "Pair Type",
    drop = FALSE
  ) +
  labs(
    title = "True vs Combined GRM",
    x = "True Relationship",
    y = "CovCombR Combined Relationship"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.1, "cm")
  )

# Arrange plots (2x2 grid)
grid.arrange(
  true_grm_plot, combined_grm_plot,
  residual_plot, scatter_plot,
  ncol = 2, nrow = 2
)

## ----recovery-stats-----------------------------------------------------------
# Calculate recovery statistics
upper_tri_idx <- upper.tri(true_grm, diag = TRUE)

# Overall statistics (only for genotypes present in combined GRM)
valid_idx <- upper_tri_idx & !is.na(combined_grm)
cat("Overall Recovery Statistics (genotypes in combined GRM):\n")
cat(sprintf(
  "  N pairs: %d (%.1f%% of total)\n",
  sum(valid_idx),
  100 * sum(valid_idx) / sum(upper_tri_idx)
))
cat(sprintf(
  "  RMSE: %.4f\n",
  sqrt(mean((combined_grm[valid_idx] - true_grm[valid_idx])^2, na.rm = TRUE))
))
cat(sprintf(
  "  MAE: %.4f\n",
  mean(abs(combined_grm[valid_idx] - true_grm[valid_idx]), na.rm = TRUE)
))
cat(sprintf(
  "  Correlation: %.4f\n",
  cor(combined_grm[valid_idx], true_grm[valid_idx], use = "complete.obs")
))

# Separate statistics for observed vs unobserved pairs
# Note: observed_pairs was already computed in the evaluate-recovery chunk
observed_idx <- observed_pairs & upper_tri_idx
unobserved_idx <- (!observed_pairs) & upper_tri_idx

# Additionally filter out NA values (genotypes not in any study)
# combined_grm may have NA for genotypes never observed
observed_valid <- observed_idx & !is.na(combined_grm)
unobserved_valid <- unobserved_idx & !is.na(combined_grm)

if (sum(observed_valid) > 0) {
  cat("\nObserved Pairs (in at least one study):\n")
  cat(sprintf(
    "  N pairs: %d (%.1f%%)\n",
    sum(observed_valid),
    100 * sum(observed_valid) / sum(upper_tri_idx)
  ))
  cat(sprintf(
    "  RMSE: %.4f\n",
    sqrt(mean((combined_grm[observed_valid] - true_grm[observed_valid])^2, na.rm = TRUE))
  ))
  cat(sprintf(
    "  MAE: %.4f\n",
    mean(abs(combined_grm[observed_valid] - true_grm[observed_valid]), na.rm = TRUE)
  ))
  cat(sprintf(
    "  Correlation: %.4f\n",
    cor(combined_grm[observed_valid], true_grm[observed_valid], use = "complete.obs")
  ))
}

if (sum(unobserved_valid) > 0) {
  cat("\nUnobserved Pairs (never jointly observed, but both genotypes in studies):\n")
  cat(sprintf(
    "  N pairs: %d (%.1f%%)\n",
    sum(unobserved_valid),
    100 * sum(unobserved_valid) / sum(upper_tri_idx)
  ))
  cat(sprintf(
    "  RMSE: %.4f\n",
    sqrt(mean((combined_grm[unobserved_valid] - true_grm[unobserved_valid])^2, na.rm = TRUE))
  ))
  cat(sprintf(
    "  MAE: %.4f\n",
    mean(abs(combined_grm[unobserved_valid] - true_grm[unobserved_valid]), na.rm = TRUE)
  ))
  cat(sprintf(
    "  Correlation: %.4f\n",
    cor(combined_grm[unobserved_valid], true_grm[unobserved_valid], use = "complete.obs")
  ))
}

# Statistics by subpopulation blocks
cat("\nRecovery by Population Structure Blocks:\n")
for (pop1 in 1:3) {
  for (pop2 in pop1:3) {
    idx1 <- which(subpop_labels == pop1)
    idx2 <- which(subpop_labels == pop2)

    block_mask <- matrix(FALSE, n_total, n_total)
    block_mask[idx1, idx2] <- TRUE
    block_mask[idx2, idx1] <- TRUE

    block_idx <- block_mask & upper_tri_idx
    block_valid <- block_idx & !is.na(combined_grm)

    if (sum(block_valid) > 0) {
      block_rmse <- sqrt(mean((combined_grm[block_valid] - true_grm[block_valid])^2, na.rm = TRUE))
      block_cor <- cor(combined_grm[block_valid], true_grm[block_valid], use = "complete.obs")

      cat(sprintf(
        "  Pop%d-Pop%d: N=%d, RMSE=%.4f, Cor=%.4f\n",
        pop1, pop2, sum(block_valid), block_rmse, block_cor
      ))
    }
  }
}

## ----session-info-------------------------------------------------------------
sessionInfo()

