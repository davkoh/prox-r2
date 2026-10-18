# Median and standard error over replicates of the RMSE of beta and of the
# mean log predictive density (MLPD) on the test data. Only replicates where
# all models converged are used.

library(dplyr)

fits <- bind_rows(lapply(snakemake@input, readRDS)) %>% filter(all_converged)
scenario <- c("n_obs", "rho", "dgp_type", "dgp_r2", "p", "model")

rmse <- fits %>%
  filter(!is.na(beta_idx)) %>%
  group_by(across(all_of(c(scenario, "rep_id", "coefficient_type", "zero_status")))) %>%
  summarise(rmse = sqrt(mean((mean - true_beta)^2, na.rm = TRUE)), .groups = "drop") %>%
  group_by(across(all_of(c(scenario, "coefficient_type", "zero_status")))) %>%
  summarise(median = median(rmse, na.rm = TRUE),
            se = sd(rmse, na.rm = TRUE) / sqrt(n()), .groups = "drop")

elpd <- fits %>%
  filter(is.na(beta_idx)) %>%
  group_by(across(all_of(c(scenario, "rep_id")))) %>%
  summarise(mlpd = mean(mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(across(all_of(scenario))) %>%
  summarise(median = median(mlpd, na.rm = TRUE),
            se = sd(mlpd, na.rm = TRUE) / sqrt(n()), .groups = "drop")

saveRDS(rmse, snakemake@output$rmse)
saveRDS(elpd, snakemake@output$elpd)
