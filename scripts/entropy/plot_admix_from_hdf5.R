#!/usr/bin/env Rscript
# plot_admix_from_hdf5.R
# Usage:
# Rscript plot_admix_from_hdf5.R inds.txt out_prefix K chain1.hdf5 chain2.hdf5 [chain3.hdf5 ...]

suppressPackageStartupMessages({
  library(rhdf5)
  library(ggplot2)
  library(reshape2)
  library(dplyr)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 4) stop("Usage: Rscript plot_admix_from_hdf5.R inds.txt out_prefix K chain1.hdf5 [chain2 ...]")

indsfile <- args[1]
out_prefix <- args[2]
Kval <- as.integer(args[3])
chains <- args[4:length(args)]

inds <- readLines(indsfile)
n_inds <- length(inds)

# robust q reader (handles different q array orientations)
read_q_mean <- function(f, Kval, n_inds){
  qraw <- h5read(f, "q")
  d <- dim(qraw)
  # find which dimension corresponds to individuals:
  ind_dim <- which(d == n_inds)
  if(length(ind_dim) == 0){
    # fallback: assume third dim is individuals (common for entropy)
    ind_dim <- which.max(d)
    if(d[ind_dim] != n_inds) warning("q dims do not match n_inds; trying to coerce")
  }
  # find K dimension (which equals Kval)
  k_dim <- which(d == Kval)
  if(length(k_dim) == 0){
    # fallback: take the smallest non-ind dimension
    other <- setdiff(1:3, ind_dim)
    k_dim <- other[which.min(d[other])]
  } else k_dim <- k_dim[1]
  iter_dim <- setdiff(1:3, c(ind_dim, k_dim))
  # permute to (iter, ind, k)
  perm <- c(iter_dim, ind_dim, k_dim)
  qperm <- aperm(qraw, perm)
  # average over iterations (1st dimension)
  qmean <- apply(qperm, c(2,3), mean) # rows = individuals, cols = K
  # ensure matrix dims n_inds x K
  qmean <- as.matrix(qmean)
  if(nrow(qmean) != n_inds) {
    if(nrow(qmean) < n_inds) qmean <- rbind(qmean, matrix(NA, n_inds-nrow(qmean), ncol(qmean)))
    else qmean <- qmean[1:n_inds,,drop=FALSE]
  }
  colnames(qmean) <- paste0("K", 1:ncol(qmean))
  qmean
}

# read all chains and average across chains
chain_qs <- lapply(chains, function(f){
  message("Reading q from ", f)
  read_q_mean(f, Kval, n_inds)
})
# average across chains (elementwise mean, handling NA)
q_avg <- Reduce(function(a,b) { 
  m <- (a + b)
  cnt <- (!is.na(a)) + (!is.na(b))
  m / pmax(cnt,1)
}, chain_qs)

# Build dataframe for plotting
df <- as.data.frame(q_avg)
df$ind <- inds
df$pop <- sub("_.*","", df$ind)               # population prefix (adjust later if needed)
df_long <- melt(df, id.vars=c("ind","pop"), variable.name="Cluster", value.name="Proportion")
# order individuals by pop then by q of cluster 1 (stable)
df_long <- df_long %>% group_by(ind) %>% mutate(meanQ1 = Proportion[Cluster==paste0("K",1)]) %>% ungroup()
ind_order <- df_long %>% distinct(ind, pop, meanQ1) %>% arrange(pop, desc(meanQ1)) %>% pull(ind)
df_long$ind <- factor(df_long$ind, levels = ind_order)

# plot
p <- ggplot(df_long, aes(x=ind, y=Proportion, fill=Cluster)) +
  geom_bar(stat="identity", width=1) +
  facet_grid(~ pop, scales="free_x", space="free_x", switch="x") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.x = element_text(angle=90, vjust=0.5)) +
  ylab("Ancestry proportion") +
  xlab("") +
  ggtitle(paste0("Admixture plot K=", Kval))

ggsave(paste0(out_prefix, "_K", Kval, "_admixture.png"), p, width=12, height=4)
# write simpleQ (transpose to individuals x clusters to match common format)
write.table(format(q_avg, digits=5), file=paste0(out_prefix, "_K", Kval, "_simpleQ.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)

message("Wrote plot and Q file: ", paste0(out_prefix, "_K", Kval, "_admixture.png"))
message("Wrote simple Q matrix: ", paste0(out_prefix, "_K", Kval, "_simpleQ.txt"))
