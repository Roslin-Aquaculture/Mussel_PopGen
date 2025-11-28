# ============================================================================
# TRIANGLE PLOTS FOR ENTROPY - CORRECTED SCALE
# Q12 scaled 0-1 to match Mytilus literature
# Produces only Rplots.pdf
# ============================================================================

library(rhdf5)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

# --- Parameters ---
k_value <- 2
n_chains <- 3
base_path <- "Full_entropy/"

# Read individual list
individuals <- readLines("inds_mussel_all.txt")
individuals <- gsub('"', '', individuals)

# --- Extract q ---
extract_q <- function(chain_num) {
  hdf5_file <- sprintf("%s/mussel_k2_chain%d.hdf5", base_path, chain_num)
  q_raw <- h5read(hdf5_file, "q")
  q_mean <- apply(q_raw, c(2, 3), mean)
  data.frame(individual = individuals, q = q_mean[2, ], chain = chain_num)
}
q_chain_list <- lapply(1:n_chains, extract_q)
q_final <- q_chain_list[[1]]
for(i in 2:length(q_chain_list)) q_final$q <- q_final$q + q_chain_list[[i]]$q
q_final$q <- q_final$q / n_chains

# --- Extract & Calculate Q12 ---
extract_Q12 <- function(chain_num) {
  hdf5_file <- sprintf("%s/mussel_k2_chain%d.hdf5", base_path, chain_num)
  zprob <- h5read(hdf5_file, "zprob")
  Q12_per_locus <- 2 * zprob[1,,] * zprob[2,,]
  Q12_values <- apply(Q12_per_locus, 1, mean, na.rm = TRUE)
  
  # Scale Q12 to 0-1 range (multiply by 2)
  Q12_scaled <- Q12_values * 2
  
  data.frame(individual = individuals, Q12 = Q12_scaled, chain = chain_num)
}
Q12_chain_list <- lapply(1:n_chains, extract_Q12)
Q12_final <- Q12_chain_list[[1]]
for(i in 2:length(Q12_chain_list)) Q12_final$Q12 <- Q12_final$Q12 + Q12_chain_list[[i]]$Q12
Q12_final$Q12 <- Q12_final$Q12 / n_chains

# --- Combine data ---
plot_data <- merge(q_final[, c("individual", "q")], 
                   Q12_final[, c("individual", "Q12")], 
                   by = "individual")

# --- Assign Populations ---
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

# --- Classify Hybrids ---
classify_hybrids <- function(q, Q12) {
  result <- rep("Admixed", length(q))
  result[q < 0.05 & Q12 < 0.1] <- "p0"
  result[q > 0.95 & Q12 < 0.1] <- "p1"
  result[abs(q - 0.5) < 0.15 & Q12 > 0.9] <- "F1"
  f2_idx <- which(abs(q - 0.5) < 0.25 & Q12 > 0.4 & Q12 < 0.7 & result == "Admixed")
  if(length(f2_idx) > 0) result[f2_idx] <- "F2"
  bc0_idx <- which(q > 0.05 & q < 0.4 & Q12 > 0.2 & Q12 < 0.6 & result == "Admixed")
  if(length(bc0_idx) > 0) result[bc0_idx] <- "BC_p0"
  bc1_idx <- which(q > 0.6 & q < 0.95 & Q12 > 0.2 & Q12 < 0.6 & result == "Admixed")
  if(length(bc1_idx) > 0) result[bc1_idx] <- "BC_p1"
  return(result)
}
plot_data$hybrid_class <- classify_hybrids(plot_data$q, plot_data$Q12)

# --- Color palette ---
class_colors <- c(
  "p0" = "#1f77b4", "p1" = "#ff7f0e", "F1" = "#2ca02c",
  "F2" = "#d62728", "BC_p0" = "#9467bd", "BC_p1" = "#8c564b", "Admixed" = "#7f7f7f"
)

# --- Plot function ---
create_triangle_plot <- function(pop_data, pop_name) {
  pop_data <- pop_data[!is.na(pop_data$q) & !is.na(pop_data$Q12), ]
  
  ggplot(pop_data, aes(x = q, y = Q12, color = hybrid_class)) +
    geom_segment(aes(x = 0, y = 0, xend = 0.5, yend = 1), 
                 color = "black", linewidth = 0.4) +
    geom_segment(aes(x = 0.5, y = 1, xend = 1, yend = 0), 
                 color = "black", linewidth = 0.4) +
    geom_point(size = 3, alpha = 0.7) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), name = "Admixture proportion (q)") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), 
                       name = expression("Interspecific heterozygosity ("*Q[12]*")")) +
    scale_color_manual(values = class_colors, name = "Hybrid class") +
    ggtitle(sprintf("%s (n=%d)", pop_name, nrow(pop_data))) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.position = "right",
      aspect.ratio = 1,
      panel.grid.minor = element_blank(),
      legend.key = element_blank()
    )
}

# --- Create plots ---
populations <- c("Aberdeen", "Cromarty", "Ireland", "M.edulis_Eng", 
                 "M.gallo_Atl", "Shetland_A", "Shetland_B", "Western_Isles")

plots_list <- list()
for(pop in populations) {
  pop_subset <- plot_data[plot_data$population == pop, ]
  if(nrow(pop_subset) > 0) {
    cat(sprintf("Creating plot for %s (%d individuals)\n", pop, nrow(pop_subset)))
    plots_list[[pop]] <- create_triangle_plot(pop_subset, pop)
  }
}

# --- Create combined plot ---
if(length(plots_list) > 0) {
  n_pops <- length(plots_list)
  n_cols <- 2
  n_rows <- ceiling(n_pops / n_cols)
  
  pdf(file = "Rplots.pdf", width = 8, height = 3.6 * n_rows)
  
  combined_plot <- grid.arrange(
    grobs = plots_list,
    ncol = n_cols,
    top = textGrob(sprintf("Triangle Plots - All Populations (K = %d)", k_value),
                   gp = gpar(fontsize = 16, fontface = "bold"))
  )
  
  print(combined_plot)
  dev.off()
  
  cat("\nSUCCESS! Combined plot saved to: Rplots.pdf\n")
}

cat(sprintf("\nProcessed %d individuals across %d populations\n", nrow(plot_data), length(plots_list)))
