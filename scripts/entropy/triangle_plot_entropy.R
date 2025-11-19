#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rhdf5)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) stop("Usage: Rscript script.R out.pdf inds.txt chain1.hdf5 [chain2.hdf5 ...]")

outfile  <- args[1]
indsfile <- args[2]
chains   <- args[3:length(args)]

inds <- read.table(indsfile, stringsAsFactors = FALSE, header = FALSE)[,1]
n_inds <- length(inds)

# ------------------------------------------------------------------
#  Population extractor: take up to TWO underscore-separated fields
# ------------------------------------------------------------------
extract_pop <- function(x) {
  parts <- strsplit(x, "_")[[1]]

  # M.edulis_Eng_* → keep first 2 fields
  # M.gallo_Atl_*  → keep first 2 fields
  if(grepl("^M\\.edulis$", parts[1]) || grepl("^M\\.gallo$", parts[1])) {
    return(paste(parts[1], parts[2], sep = "_"))
  }

  # Western_Isles_* → keep first 2 fields
  if(parts[1] == "Western" && parts[2] == "Isles") {
    return("Western_Isles")
  }

  # Shetland_A_* → keep first 2 fields
  # Shetland_B_* → keep first 2 fields
  if(parts[1] == "Shetland") {
    return(paste(parts[1], parts[2], sep = "_"))
  }

  # Aberdeen_Aber* → first field only
  # Ireland_ATU*   → first field only
  # Cromarty_C*    → first field only
  return(parts[1])
}


# ------------------------------------------------------------------
# Readers
# ------------------------------------------------------------------

read_q <- function(f) {
  qraw <- h5read(f, "q")
  d <- dim(qraw)

  k_dim   <- which(d == 2)
  if(length(k_dim) == 0) k_dim <- which.min(d)[1] else k_dim <- k_dim[1]

  ind_dim <- which(d == n_inds)
  if(length(ind_dim) == 0) {
    other <- setdiff(1:3, k_dim)
    ind_dim <- other[which.min(d[other])]
  } else ind_dim <- ind_dim[1]

  iter_dim <- setdiff(1:3, c(k_dim, ind_dim))

  qmean <- apply(qraw, c(k_dim, ind_dim), mean)
  t(qmean)
}

read_Q12 <- function(f) {
  z <- h5read(f, "zprob")
  d <- dim(z)

  k_dim <- which(d == 2)
  if(length(k_dim) == 0) k_dim <- 1

  ind_dim <- which(d == n_inds)
  if(length(ind_dim) == 0) {
    other <- setdiff(1:3, k_dim)
    ind_dim <- other[which.min(d[other])]
  } else ind_dim <- ind_dim[1]

  loci_dim <- setdiff(1:3, c(k_dim, ind_dim))
  z_perm <- aperm(z, c(k_dim, ind_dim, loci_dim))

  K      <- dim(z_perm)[1]
  n_ind  <- dim(z_perm)[2]
  n_loci <- dim(z_perm)[3]

  Q12 <- numeric(n_ind)
  for(i in seq_len(n_ind)) {
    p1 <- z_perm[1, i, ]
    p2 <- z_perm[2, i, ]
    Q12[i] <- mean(p1 * p2, na.rm = TRUE)
  }
  Q12
}

classify <- function(q, Q12) {
  if(is.na(q) || is.na(Q12)) return("Other")
  if(Q12 > 0.75 && q > 0.45 && q < 0.55) return("F1")
  if(Q12 > 0.25 && q > 0.35 && q < 0.65) return("F2/F3")
  if(Q12 < 0.25 && q > 0.1  && q < 0.9 ) return("Backcross")
  if(q <= 0.05) return("Parental_A")
  if(q >= 0.95) return("Parental_B")
  "Other"
}

# ------------------------------------------------------------------
# Read all chains
# ------------------------------------------------------------------

all_rows <- list()

for(f in chains) {
  message("Reading ", f)
  qmat <- read_q(f)
  Q12  <- read_Q12(f)

  if(nrow(qmat) != n_inds) {
    warning("q matrix individuals != inds file. Adjusting...")
    if(nrow(qmat) > n_inds) qmat <- qmat[1:n_inds, , drop=FALSE]
    else qmat <- rbind(qmat, matrix(NA, n_inds - nrow(qmat), ncol(qmat)))
  }

  df <- data.frame(
    ind    = inds,
    q      = qmat[,1],
    Q12    = Q12,
    file   = basename(f),
    stringsAsFactors = FALSE
  )

  # use new population extractor
  df$pop <- vapply(df$ind, extract_pop, FUN.VALUE = character(1))
  df$hybridclass <- mapply(classify, df$q, df$Q12)

  all_rows[[f]] <- df
}

df_all <- do.call(rbind, all_rows)

# ------------------------------------------------------------------
# Final population levels (stable order)
# ------------------------------------------------------------------
pop_from_inds <- vapply(inds, extract_pop, FUN.VALUE = character(1))
df_all$pop <- factor(df_all$pop, levels = unique(pop_from_inds))

# ------------------------------------------------------------------
# User-defined population colours (must match extracted pop names)
# ------------------------------------------------------------------

mussel_colors <- c(
  "Shetland_A"   = "#FF69B4",
  "Shetland_B"   = "#8B00FF",
  "Ireland"      = "#FFA500",
  "Cromarty"     = "#1B9E77",
  "M.gallo_Atl"  = "#D95F02",
  "M.edulis_Eng" = "#7570B3",
  "Western_Isles"= "#66A61E",
  "Aberdeen"     = "#1E7878"
)

# ------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------

# ----- Draw triangle boundaries -----
triangle_df <- data.frame(
  x = c(0, 1, 0.5, 0),
  y = c(0, 0, 1, 0)
)

p <- ggplot() +
  geom_path(data = triangle_df, aes(x=x, y=y), linewidth=1.1, colour="black") +
  geom_point(data = df_all,
             aes(x=q, y=Q12, colour=pop),
             size = 2.2, alpha = 0.9,
             position = position_jitter(width=0.02, height=0.02)) +
  scale_colour_manual(values = mussel_colors) +
  xlab("Ancestry proportion (q)") +
  ylab("Inter-ancestry heterozygosity (Q12)") +
  xlim(0,1) + ylim(0,1) +
  theme_bw(base_size=14)


ggsave(outfile, p, width=10, height=10)

csvout <- sub("\\.pdf$", ".triangle_summary.csv", outfile)
write.csv(df_all, csvout, row.names = FALSE)

cat("Saved triangle plot to:", outfile, "\n")
cat("Saved CSV to:", csvout, "\n")
