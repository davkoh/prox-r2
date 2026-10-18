# Simulate all replicates of one scenario.

snakemake@source("dgp.R")

cfg <- snakemake@config
w <- snakemake@wildcards

reps <- lapply(seq_len(cfg$n_reps), function(i) {
  simulate_replicate(n = as.numeric(w$n_obs), p = as.numeric(w$p),
                     rho = as.numeric(w$rho), r2 = as.numeric(w$r2),
                     dgp_type = w$dgp, seed = cfg$seed + i)
})

saveRDS(reps, snakemake@output[[1]])
