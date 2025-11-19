#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rhdf5)
  library(coda)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript convergence_diag_fixed.R outprefix K chain1.hdf5 chain2.hdf5 ...")
}

out <- args[1]
K <- as.integer(args[2])
chains <- args[3:length(args)]

read_chain_ll <- function(file) {
  message("Reading likelihood trace from ", file)
  if (!H5Lexists(H5Fopen(file), "likdata")) {
    stop("File ", file, " does not contain likdata dataset.")
  }
  # load only first column ? log-posterior / log-likelihood
  lik <- h5read(file, "likdata")[,1]
  return(mcmc(lik))
}

mcmc_list <- lapply(chains, read_chain_ll)
mcmc_list <- mcmc.list(mcmc_list)

# Compute Gelman–Rubin diagnostic
gr <- gelman.diag(mcmc_list, autoburnin=FALSE)

# Save results
sink(paste0(out, "_convergence.txt"))
print(gr)
sink()

cat("Saved convergence diagnostics to ", out, "_convergence.txt\n")
