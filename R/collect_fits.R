# Combine the fits of one scenario and add the true coefficients.

library(dplyr)

w <- snakemake@wildcards
p <- as.integer(w$p)
reps <- readRDS(snakemake@input$data)
fits <- bind_rows(lapply(snakemake@input$csv, read.csv))

# A replicate is converged if the mean R-hat over all summarised quantities
# is below 1.1 for every model.
converged <- fits %>%
  group_by(model, rep_id) %>%
  summarise(ok = mean(rhat, na.rm = TRUE) < 1.1, .groups = "drop") %>%
  group_by(rep_id) %>%
  summarise(all_converged = all(ok), .groups = "drop")

truth <- bind_rows(lapply(seq_along(reps), function(i) {
  b <- unname(c(reps[[i]]$beta_main, reps[[i]]$beta_int))
  tibble(rep_id = i, variable = sprintf("beta[%d]", seq_along(b)), true_beta = b)
}))

out <- fits %>%
  filter(grepl("^(beta|log_lik_test)\\[\\d+\\]$", variable)) %>%
  left_join(truth, by = c("rep_id", "variable")) %>%
  left_join(converged, by = "rep_id") %>%
  mutate(
    beta_idx = ifelse(startsWith(variable, "beta"), as.integer(gsub("\\D", "", variable)), NA),
    coefficient_type = ifelse(beta_idx <= p, "main", "interaction"),
    zero_status = ifelse(true_beta != 0, "true non-zero", "true zero"),
    dgp_type = w$dgp, n_obs = as.numeric(w$n_obs), p = p,
    rho = as.numeric(w$rho), dgp_r2 = as.numeric(w$r2)
  )

saveRDS(out, snakemake@output[[1]])
