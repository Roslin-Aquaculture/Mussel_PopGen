# ============================================================================
# TRIANGLE PLOTS FOR ENTROPY - K=2
# FIXED: Proper Q12 calculation from zprob
# ============================================================================

library(rhdf5)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

# --- Check packages ---
required_packages <- c("rhdf5", "ggplot2", "dplyr", "grid", "gridExtra")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if(length(missing_packages) > 0) {
  cat("Missing packages:", paste(missing_packages, collapse = ", "), "\n")
  stop("Install missing packages")
}

# ============================================================================
# 1. EXTRACT DATA - PROPER Q12 FROM ZPROB
# ============================================================================

k_value <- 2
n_chains <- 3
base_path <- "Full_entropy/"
output_dir <- "triangle_plots_corrected"
dir.create(output_dir, showWarnings = FALSE)

# Read individual list
individuals <- readLines("inds_mussel_all.txt")
individuals <- gsub('"', '', individuals)

# Extract q (this part was correct)
extract_q_data <- function(chain_num, individuals) {
  hdf5_file <- sprintf("%s/mussel_k2_chain%d.hdf5", base_path, chain_num)
  q_array <- h5read(hdf5_file, "q")
  q_mean <- apply(q_array, c(2, 3), mean)
  data.frame(individual = individuals, q = q_mean[2, ], chain = chain_num)
}

chain_q_list <- lapply(1:n_chains, function(ch) extract_q_data(ch, individuals))
combined_q <- chain_q_list[[1]]
for(i in 2:length(chain_q_list)) combined_q$q <- combined_q$q + chain_q_list[[i]]$q
combined_q$q <- combined_q$q / n_chains

# ============================================================================
# FIX: PROPER Q12 CALCULATION FROM ZPROB
# ============================================================================
# zprob dimensions: 2 x 306 x 6802
# zprob[1,,] = probability allele from species 1 (p0)
# zprob[2,,] = probability allele from species 1 (p1)
# Q12 = 2 * p0 * p1 (probability of interspecific heterozygote)

extract_Q12_data <- function(chain_num, individuals) {
  hdf5_file <- sprintf("%s/mussel_k2_chain%d.hdf5", base_path, chain_num)
  
  # Read zprob: species-specific allele probabilities
  zprob <- h5read(hdf5_file, "zprob")
  
  cat(sprintf("  zprob dimensions: %d x %d x %d\n", dim(zprob)[1], dim(zprob)[2], dim(zprob)[3]))
  
  # Calculate Q12 per locus: Q12 = 2 * p_species1 * p_species2
  # zprob[1,,] is p_species1, zprob[2,,] is p_species2
  Q12_per_locus <- 2 * zprob[1,,] * zprob[2,,]
  
  # Average across loci for each individual
  # Apply across columns (loci) to get per-individual mean
  Q12_values <- apply(Q12_per_locus, 1, mean, na.rm = TRUE)
  
  # Verify we have the right number of values
  if(length(Q12_values) != length(individuals)) {
    stop(sprintf("Dimension mismatch: Q12 has %d values but %d individuals", 
                 length(Q12_values), length(individuals)))
  }
  
  data.frame(individual = individuals, Q12 = Q12_values, chain = chain_num)
}

# Extract Q12 from each chain
chain_Q12_list <- lapply(1:n_chains, function(ch) extract_Q12_data(ch, individuals))

# Check that Q12 is different from q
cat("\nQ12 summary from chain 1:\n")
print(summary(chain_Q12_list[[1]]$Q12))

cat("\nq summary from chain 1:\n")
print(summary(chain_q_list[[1]]$q))

cat("\nCorrelation between q and Q12 (should be moderate, < 0.9):\n")
cat(sprintf("  r = %.3f\n", cor(chain_q_list[[1]]$q, chain_Q12_list[[1]]$Q12)))

# Average across chains
combined_Q12 <- chain_Q12_list[[1]]
for(i in 2:length(chain_Q12_list)) combined_Q12$Q12 <- combined_Q12$Q12 + chain_Q12_list[[i]]$Q12
combined_Q12$Q12 <- combined_Q12$Q12 / n_chains

# Merge q and Q12
plot_data <- merge(combined_q[, c("individual", "q")], 
                   combined_Q12[, c("individual", "Q12")], 
                   by = "individual")

# ============================================================================
# 2. VERIFY DATA IS CORRECT
# ============================================================================

cat("\n" , rep("=", 60), "\n", sep="")
cat("FINAL DATA VERIFICATION\n")
cat(rep("=", 60), "\n", sep="")

cat("\nFinal q summary:\n")
print(summary(plot_data$q))

cat("\nFinal Q12 summary:\n")
print(summary(plot_data$Q12))

cat("\nFinal correlation between q and Q12:\n")
cor_final <- cor.test(plot_data$q, plot_data$Q12)
cat(sprintf("  r = %.3f, p = %.2e\n", cor_final$estimate, cor_final$p.value))

if(cor_final$estimate > 0.85) {
  cat("\n*** WARNING: q and Q12 are still too correlated! ***\n")
  cat("Something is wrong with the Q12 extraction.\n")
} else {
  cat("\n? q and Q12 correlation looks reasonable.\n")
}

# Save raw data for inspection
write.csv(plot_data, file.path(output_dir, "raw_data_extracted.csv"), row.names = FALSE)

# ============================================================================
# 3. ASSIGN POPULATIONS
# ============================================================================

plot_data$population <- NA
plot_data$population[grepl("^Aberdeen_", plot_data$individual)] <- "Aberdeen"
plot_data$population[grepl("^Cromarty_", plot_data$individual)] <- "Cromarty"
plot_data$population[grepl("^Ireland_", plot_data$individual)] <- "Ireland"
plot_data$population[grepl("^M\\.edulis_Eng_", plot_data$individual)] <- "M.edulis_Eng"
plot_data$population[grepl("^M\\.gallo_Atl_", plot_data$individual)] <- "M.gallo_Atl"
plot_data$population[grepl("^Shetland_A_", plot_data$individual)] <- "Shetland_A"
plot_data$population[grepl("^Shetland_B_", plot_data$individual)] <- "Shetland_B"
plot_data$population[grepl("^Western_Isles_", plot_data$individual)] <- "Western_Isles"

plot_data <- plot_data[!is.na(plot_data$population), ]

cat("\nFinal population counts:\n")
print(table(plot_data$population))

# ============================================================================
# 4. CLASSIFY HYBRIDS
# ============================================================================

classify_hybrids <- function(q, Q12) {
  result <- rep("Admixed", length(q))
  
  # Parental species (strict thresholds)
  result[q < 0.05 & Q12 < 0.05] <- "p0"
  result[q > 0.95 & Q12 < 0.05] <- "p1"
  
  # F1 hybrids (high Q12, q near 0.5)
  result[abs(q - 0.5) < 0.15 & Q12 > 0.5] <- "F1"
  
  # F2 hybrids
  f2_idx <- which(abs(q - 0.5) < 0.25 & Q12 > 0.2 & Q12 < 0.5)
  if(length(f2_idx) > 0) result[f2_idx] <- "F2"
  
  # Backcrosses
  bc0_idx <- which(q > 0.05 & q < 0.4 & Q12 > 0.1 & Q12 < 0.3)
  if(length(bc0_idx) > 0) result[bc0_idx] <- "BC_p0"
  
  bc1_idx <- which(q > 0.6 & q < 0.95 & Q12 > 0.1 & Q12 < 0.3)
  if(length(bc1_idx) > 0) result[bc1_idx] <- "BC_p1"
  
  return(result)
}

plot_data$hybrid_class <- classify_hybrids(plot_data$q, plot_data$Q12)

cat("\nHybrid classification:\n")
print(table(plot_data$hybrid_class))

# ============================================================================
# 5. CREATE TRIANGLE PLOTS
# ============================================================================

class_colors <- c(
  "p0" = "#1f77b4", "p1" = "#ff7f0e", "F1" = "#2ca02c",
  "F2" = "#d62728", "BC_p0" = "#9467bd", "BC_p1" = "#8c564b", "Admixed" = "#7f7f7f"
)

create_triangle_plot <- function(pop_data, pop_name) {
  pop_data <- pop_data[!is.na(pop_data$q) & !is.na(pop_data$Q12), ]
  
  ggplot(pop_data, aes(x = q, y = Q12)) +
    stat_function(fun = function(x) 2 * x * (1 - x), 
                  color = "black", linetype = "dashed", linewidth = 0.7) +
    geom_point(aes(color = hybrid_class, fill = hybrid_class), 
               size = 3, alpha = 0.7, shape = 21, stroke = 0.5) +
    geom_rug(aes(color = hybrid_class), sides = "b", alpha = 0.3) +
    geom_rug(aes(color = hybrid_class), sides = "l", alpha = 0.3) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), name = "Admixture proportion (q)") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), name = "Interspecific heterozygosity (Q12)") +
    scale_color_manual(values = class_colors, name = "Hybrid class") +
    scale_fill_manual(values = class_colors) +
    ggtitle(sprintf("%s (n=%d)", pop_name, nrow(pop_data))) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      legend.position = "right",
      aspect.ratio = 1
    )
}

# Create plots
populations <- c("Aberdeen", "Cromarty", "Ireland", "M.edulis_Eng", 
                 "M.gallo_Atl", "Shetland_A", "Shetland_B", "Western_Isles")

plots_list <- list()
for(pop in populations) {
  pop_subset <- plot_data[plot_data$population == pop, ]
  if(nrow(pop_subset) > 0) {
    cat(sprintf("Plot: %s (%d individuals)\n", pop, nrow(pop_subset)))
    plots_list[[pop]] <- create_triangle_plot(pop_subset, pop)
    ggsave(file.path(output_dir, sprintf("triangle_%s.png", pop)), 
           plots_list[[pop]], width = 7, height = 6, dpi = 300)
  }
}

# Combined plot
if(length(plots_list) > 0) {
  cat("Creating combined plot...\n")
  n_pops <- length(plots_list)
  n_cols <- 2
  n_rows <- ceiling(n_pops / n_cols)
  
  combined_plot <- grid.arrange(
    grobs = plots_list,
    ncol = n_cols,
    top = textGrob(sprintf("Triangle Plots - All Populations (K = %d)", k_value),
                   gp = gpar(fontsize = 16, fontface = "bold"))
  )
  
  ggsave(file.path(output_dir, "all_populations.png"), 
         combined_plot, width = 14, height = 6 * n_rows, dpi = 300)
}

# ============================================================================
# 6. SAVE RESULTS
# ============================================================================

write.csv(plot_data, file.path(output_dir, "plot_data_final.csv"), row.names = FALSE)

summary_table <- plot_data %>%
  group_by(population, hybrid_class) %>%
  summarise(count = n(), mean_q = mean(q), mean_Q12 = mean(Q12))

write.csv(summary_table, file.path(output_dir, "hybrid_summary.csv"), row.names = FALSE)

cat("\n", rep("=", 60), "\n", sep="")
cat("? COMPLETE! Check output in:", output_dir, "\n")
cat(rep("=", 60), "\n", sep="")
