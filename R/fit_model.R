# Fit one model to one batch of replicates and save posterior summaries.

library(cmdstanr)
library(posterior)

cfg <- snakemake@config
model <- snakemake@wildcards$model
batch <- as.integer(snakemake@wildcards$batch)

reps <- readRDS(snakemake@input$data)
mod <- cmdstan_model(exe_file = snakemake@input$exe)

# Each model reads only the data that its data block declares.
stan_data <- function(r) {
  p <- ncol(r$X)
  c(list(N = length(r$y), y = r$y, p = p, K = p + p * (p - 1) / 2,
         X = cbind(1, r$X, r$X_int),
         Ntest = length(r$y_test), ytest = r$y_test,
         Xtest = cbind(1, r$X_test, r$X_int_test),
         scale_sigma = sd(r$y)),
    cfg$hyper)
}

n_per_batch <- cfg$n_reps / cfg$n_batches
rep_ids <- (batch - 1) * n_per_batch + seq_len(n_per_batch)

summaries <- lapply(rep_ids, function(i) {
  fit <- mod$sample(
    data = stan_data(reps[[i]]),
    chains = cfg$mcmc$chains,
    parallel_chains = snakemake@threads,
    iter_warmup = cfg$mcmc$iter_warmup,
    iter_sampling = cfg$mcmc$iter_sampling,
    max_treedepth = cfg$mcmc$max_treedepth,
    adapt_delta = cfg$mcmc$adapt_delta,
    seed = cfg$seed,
    refresh = 0
  )
  s <- summarise_draws(
    fit$draws(c("beta", "sigma", "log_lik_test", "y_tilde_test")),
    default_summary_measures(), default_convergence_measures()
  )
  s$model <- model
  s$rep_id <- i
  s
})

readr::write_csv(dplyr::bind_rows(summaries), snakemake@output[[1]])
