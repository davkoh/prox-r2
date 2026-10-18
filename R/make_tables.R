# LaTeX tables for beta: MCMC diagnostics and coverage of the central 90%
# posterior interval [q5, q95]. These tables use all replicates.

library(dplyr)
library(tidyr)

cfg <- snakemake@config
model_labels <- unlist(cfg$models)
models <- names(model_labels)
dgp_labels <- c(sparse_correct = "Sparse, correct", dense_correct = "Dense, correct",
                sparse_misp = "Sparse, misp.", dense_misp = "Dense, misp.")
zero_levels <- c("true non-zero", "true zero")
coef_labels <- c(main = "main", interaction = "int.")

fmt_rhat <- function(x) {
  if (is.na(x)) return("--")
  s <- sprintf("%.3f", x)
  if (x > 1.01) paste0("\\textbf{", s, "}") else s
}
fmt_ess <- function(x) {
  if (is.na(x)) return("--")
  format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
}
fmt_cov <- function(x) {
  if (is.na(x)) return("--")
  s <- sprintf("%.2f", x)
  if (abs(x - 0.90) > 0.05) paste0("\\textbf{", s, "}") else s
}

# Body rows grouped by DGP. row_cells(sub, i) gives the cells after the DGP column.
table_body <- function(wide, row_cells) {
  dgps <- intersect(cfg$dgp_types, wide$dgp_type)
  unlist(lapply(seq_along(dgps), function(d) {
    sub <- wide[wide$dgp_type == dgps[d], ]
    rows <- vapply(seq_len(nrow(sub)), function(i) {
      first <- if (i == 1) sprintf("\\multirow{%d}{*}{%s}", nrow(sub), dgp_labels[[dgps[d]]]) else ""
      paste0(paste(c(first, row_cells(sub, i)), collapse = " & "), " \\\\")
    }, "")
    if (d > 1) c("\\midrule", rows) else rows
  }))
}

# Landscape table with two columns per model.
write_table <- function(body, lead_headers, sub_headers, caption, label, file) {
  n_lead <- length(lead_headers)
  first_col <- n_lead + 2 * seq_along(models) - 1
  header1 <- paste(c(rep("", n_lead), sprintf("\\multicolumn{2}{c}{%s}", model_labels)),
                   collapse = " & ")
  header2 <- paste(c(lead_headers, rep(sub_headers, length(models))), collapse = " & ")
  writeLines(c(
    "% LaTeX packages required: booktabs, multirow, pdflscape",
    "\\begin{landscape}",
    "\\begin{table}[ht]",
    "\\centering",
    "\\small",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    sprintf("\\begin{tabular}{@{}%s%s@{}}", strrep("l", n_lead), strrep("rr", length(models))),
    "\\toprule",
    paste0(header1, " \\\\"),
    paste(sprintf("\\cmidrule(lr){%d-%d}", first_col, first_col + 1), collapse = ""),
    paste0(header2, " \\\\"),
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "\\end{landscape}"
  ), file)
}

write_diagnostics <- function(diag, caption, label, file) {
  wide <- diag %>%
    mutate(model = factor(model, models)) %>%
    pivot_wider(names_from = model, values_from = c(max_rhat, min_ess), names_expand = TRUE) %>%
    arrange(dgp_r2)
  body <- table_body(wide, function(sub, i) {
    c(sprintf("%.1f", sub$dgp_r2[i]),
      unlist(lapply(models, function(m) {
        c(fmt_rhat(sub[[paste0("max_rhat_", m)]][i]), fmt_ess(sub[[paste0("min_ess_", m)]][i]))
      })))
  })
  write_table(body, c("DGP", "$R^2$"), c("$\\bar{\\hat{R}}$", "$\\overline{\\text{ESS}}$"),
              caption, label, file)
}

write_coverage <- function(cov, caption, label, file) {
  wide <- cov %>%
    mutate(model = factor(model, models), zero_status = factor(zero_status, zero_levels),
           coefficient_type = factor(coefficient_type, names(coef_labels))) %>%
    pivot_wider(names_from = c(model, zero_status), values_from = coverage,
                names_sep = "__", names_expand = TRUE) %>%
    arrange(dgp_r2, coefficient_type)
  body <- table_body(wide, function(sub, i) {
    r2 <- sub$dgp_r2
    r2_cell <- if (i == 1 || r2[i] != r2[i - 1]) {
      sprintf("\\multirow{%d}{*}{%.1f}", sum(r2 == r2[i]), r2[i])
    } else ""
    c(r2_cell, coef_labels[[as.character(sub$coefficient_type[i])]],
      unlist(lapply(models, function(m) {
        vapply(zero_levels, function(z) fmt_cov(sub[[paste(m, z, sep = "__")]][i]), "")
      })))
  })
  write_table(body, c("DGP", "$R^2$", "Coef."), c("non-0", "zero"), caption, label, file)
}

betas <- bind_rows(lapply(snakemake@input, readRDS)) %>% filter(!is.na(beta_idx))
scenario <- c("n_obs", "rho", "dgp_type", "dgp_r2", "p", "model")

# Mean over replicates of the maximum R-hat and the minimum bulk ESS over beta
diag <- betas %>%
  group_by(across(all_of(c(scenario, "rep_id")))) %>%
  summarise(max_rhat = max(rhat, na.rm = TRUE), min_ess = min(ess_bulk, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(across(all_of(scenario))) %>%
  summarise(max_rhat = mean(max_rhat, na.rm = TRUE), min_ess = mean(min_ess, na.rm = TRUE),
            .groups = "drop")

cov <- betas %>%
  group_by(across(all_of(c(scenario, "coefficient_type", "zero_status")))) %>%
  summarise(coverage = mean(q5 <= true_beta & true_beta <= q95, na.rm = TRUE), .groups = "drop")

diag_caption <- paste0(
  "For each DGP condition and model: mean (over simulation repetitions) of the ",
  "per-repetition maximum $\\hat{R}$ and minimum bulk ESS across $\\beta$ coefficients")
settings <- distinct(betas, n_obs, rho)

for (i in seq_len(nrow(settings))) {
  s <- settings[i, ]
  suffix <- sprintf("nobs_%s_rho_%s", s$n_obs, s$rho)
  out <- function(name) file.path("output/tables", sprintf("%s_%s.tex", name, suffix))
  diag_s <- semi_join(diag, s, by = c("n_obs", "rho"))
  cov_s <- semi_join(cov, s, by = c("n_obs", "rho"))

  write_diagnostics(
    diag_s %>%
      group_by(dgp_type, dgp_r2, model) %>%
      summarise(max_rhat = max(max_rhat, na.rm = TRUE), min_ess = min(min_ess, na.rm = TRUE),
                .groups = "drop"),
    caption = paste0(
      "MCMC diagnostics for $\\boldsymbol{\\beta}$. ", diag_caption,
      ", aggregated as worst-case over all values of $p$. ",
      "Bold entries indicate $\\bar{\\hat{R}} > 1.01$."),
    label = "tab:mcmc_diagnostics_beta",
    file = out("diagnostics_beta"))

  for (p_val in sort(unique(diag_s$p))) {
    write_diagnostics(
      diag_s %>% filter(p == p_val) %>% select(dgp_type, dgp_r2, model, max_rhat, min_ess),
      caption = paste0(
        "MCMC diagnostics for $\\boldsymbol{\\beta}$ ($p = ", p_val, "$). ", diag_caption,
        ". Bold entries indicate $\\bar{\\hat{R}} > 1.01$."),
      label = sprintf("tab:mcmc_diagnostics_beta_p%d", p_val),
      file = out(sprintf("diagnostics_beta_p%d", p_val)))

    write_coverage(
      cov_s %>% filter(p == p_val) %>%
        select(dgp_type, dgp_r2, model, coefficient_type, zero_status, coverage),
      caption = paste0(
        "Empirical coverage of the central 90\\% posterior interval ",
        "$[q_{5}, q_{95}]$ for $\\boldsymbol{\\beta}$ ($p = ", p_val, "$). ",
        "Coverage is pooled over coefficients and simulation repetitions, ",
        "split by coefficient type (main vs.\\ interaction) and zero status ",
        "(true non-zero vs.\\ true zero). ",
        "Bold entries deviate from the nominal 0.90 by more than 0.05."),
      label = sprintf("tab:coverage_beta_p%d", p_val),
      file = out(sprintf("coverage_beta_p%d", p_val)))
  }
}
