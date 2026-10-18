# Figures of RMSE and ELPD relative to the benchmark model, one block per DGP
# in a 2 x 2 grid.

library(dplyr)
library(ggplot2)
library(patchwork)

cfg <- snakemake@config
bench <- cfg$benchmark_model
model_labels <- unlist(cfg$models)
models <- setdiff(names(model_labels), bench)
model_colors <- c(prox_r2 = "#F8766D", normal_gamma = "#7CAE00", r2d2 = "#00BFC4", rhs = "#C77CFF")
zero_linetypes <- c("true non-zero" = "solid", "true zero" = "22")
zero_shapes <- c("true non-zero" = 16, "true zero" = 17)
dgp_labels <- c(
  sparse_correct = "Sparse, correctly specified (H1)",
  dense_correct = "Dense, correctly specified (H1)",
  sparse_misp = "Sparse, misspecified (H2)",
  dense_misp = "Dense, misspecified (H2)"
)
label_r2 <- function(x) paste0("DGP R2 = ", x)
fmt <- function(digits) function(x) ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))

# Breaks on a grid of 0.1 or coarser, so that one-decimal labels are distinct.
breaks_one_decimal <- function(limits) {
  span <- diff(limits)
  if (!is.finite(span) || span <= 0) return(round(limits, 1))
  steps <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
  step <- steps[which.min(abs(span / steps - 4))]
  seq(ceiling(limits[1] / step) * step, floor(limits[2] / step) * step, by = step)
}

# Axis titles, axis labels and strips show only on the outer edges of the grid.
block_theme <- function(left, right, top, bottom) {
  axis_line <- element_line(color = "grey30", linewidth = 0.3)
  theme_minimal(base_family = "sans", base_size = 15.4) +
    theme(
      plot.title = element_text(size = 16.5, face = "bold", hjust = 0.5),
      axis.title.y = if (left) element_text(size = 14.3) else element_blank(),
      axis.text.y = element_text(size = 12.7),
      axis.text.x = if (bottom) element_text(size = 12.7) else element_blank(),
      axis.line.y = axis_line,
      axis.ticks.x = if (bottom) axis_line else element_blank(),
      axis.ticks.y = axis_line,
      strip.text.x = if (top) element_text(size = 14.3) else element_blank(),
      strip.text.y = if (right) element_text(size = 14.3) else element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.13, "in"),
      legend.position = "none",
      plot.margin = margin(3, 3, 3, 3)
    )
}

# Dotted lines show the estimate +/- SE.
se_lines <- list(
  geom_line(aes(y = rel + rel_se), linetype = "dotted", linewidth = 0.55),
  geom_line(aes(y = rel - rel_se), linetype = "dotted", linewidth = 0.55)
)

rmse_block <- function(df) {
  ggplot(df, aes(p, rel, color = model, group = interaction(model, zero_status))) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey45", linewidth = 0.3) +
    se_lines +
    geom_line(aes(linetype = zero_status), linewidth = 0.55) +
    geom_point(aes(shape = zero_status), size = 1.15, stroke = 0.2) +
    facet_grid(coefficient_type ~ dgp_r2, scales = "free_y",
               labeller = labeller(dgp_r2 = label_r2)) +
    scale_y_continuous(breaks = breaks_one_decimal, labels = fmt(1)) +
    scale_linetype_manual(values = zero_linetypes) +
    scale_shape_manual(values = zero_shapes) +
    coord_cartesian(clip = "off")
}

# The y range is at least [-4, 4] and at most [-50, max].
elpd_block <- function(df) {
  ylim <- c(max(-50, min(-4, min(df$rel, na.rm = TRUE))), max(4, max(df$rel, na.rm = TRUE)))
  ggplot(df, aes(p, rel, color = model, group = model)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey45", linewidth = 0.3) +
    geom_hline(yintercept = c(-4, 4), linetype = "dashed", color = "grey45", linewidth = 0.3) +
    se_lines +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.15, stroke = 0.2) +
    facet_wrap(~dgp_r2, nrow = 1, labeller = labeller(dgp_r2 = label_r2)) +
    scale_y_continuous(labels = fmt(0)) +
    coord_cartesian(ylim = ylim, clip = "on")
}

legend_row <- function(show_zero_status) {
  m <- tibble(model = factor(models, models), x = seq_along(models), label = model_labels[models])
  g <- ggplot() +
    geom_segment(data = m, aes(x = x - 0.28, xend = x + 0.28, y = 1, yend = 1, color = model),
                 linewidth = 0.75) +
    geom_point(data = m, aes(x = x, y = 1, color = model), size = 1.4) +
    geom_text(data = m, aes(x = x, y = 0.72, label = label), size = 5.03, family = "sans") +
    scale_color_manual(values = model_colors)
  x_max <- length(models) + 0.8
  if (show_zero_status) {
    z <- tibble(status = names(zero_shapes), x = length(models) + c(1.0, 2.1))
    g <- g +
      geom_segment(data = z, aes(x = x - 0.28, xend = x + 0.28, y = 1, yend = 1, linetype = status),
                   color = "grey15", linewidth = 0.75) +
      geom_point(data = z, aes(x = x, y = 1, shape = status), color = "grey15", size = 1.4) +
      geom_text(data = z, aes(x = x, y = 0.72, label = status), size = 5.03, family = "sans") +
      scale_linetype_manual(values = zero_linetypes) +
      scale_shape_manual(values = zero_shapes)
    x_max <- length(models) + 2.7
  }
  g + coord_cartesian(xlim = c(0.2, x_max), ylim = c(0.55, 1.45), clip = "off") +
    theme_void(base_family = "sans") +
    theme(legend.position = "none", plot.margin = margin(0, 4, 0, 4))
}

plot_grid <- function(rel, block_fun, y_title, show_zero_status) {
  dgps <- intersect(cfg$dgp_types, unique(rel$dgp_type))
  blocks <- lapply(seq_along(dgps), function(i) {
    df <- filter(rel, dgp_type == dgps[i])
    block_fun(df) +
      scale_x_continuous(breaks = sort(unique(df$p))) +
      scale_color_manual(values = model_colors) +
      labs(title = dgp_labels[[dgps[i]]], x = NULL, y = y_title) +
      block_theme(left = i %% 2 == 1, right = i %% 2 == 0,
                  top = i <= 2, bottom = i > length(dgps) - 2)
  })
  wrap_plots(wrap_plots(blocks, ncol = 2), legend_row(show_zero_status),
             ncol = 1, heights = c(1, 0.12))
}

save_figure <- function(plot, file, width, height, tikz_width, tikz_height) {
  ggsave(paste0(file, ".png"), plot, width = width, height = height, dpi = 400, bg = "white")
  ggsave(paste0(file, ".pdf"), plot, width = width, height = height,
         device = grDevices::pdf, bg = "white")
  # The TikZ version is optional. It needs tikzDevice and a LaTeX installation.
  if (requireNamespace("tikzDevice", quietly = TRUE)) {
    tryCatch({
      tikzDevice::tikz(paste0(file, ".tex"), width = tikz_width, height = tikz_height)
      print(plot)
    }, error = function(e) warning("TikZ output failed: ", conditionMessage(e)),
    finally = if (grDevices::dev.cur() > 1) grDevices::dev.off())
  }
}

# Add the benchmark median as column `bench` and keep the other models.
relative <- function(summary, keys) {
  bench_median <- summary %>% filter(model == bench) %>% select(all_of(keys), bench = median)
  summary %>%
    filter(model %in% models, dgp_r2 %in% cfg$plot_dgp_r2) %>%
    left_join(bench_median, by = keys) %>%
    mutate(model = factor(model, models), dgp_r2 = factor(dgp_r2, sort(unique(dgp_r2))))
}

rmse <- readRDS(snakemake@input$rmse)
elpd <- readRDS(snakemake@input$elpd)
bench_label <- model_labels[[bench]]
settings <- distinct(rmse, n_obs, rho)

for (i in seq_len(nrow(settings))) {
  s <- settings[i, ]
  suffix <- sprintf("rel_%s_nobs_%s_rho_%s", bench, s$n_obs, s$rho)

  rmse_rel <- rmse %>%
    semi_join(s, by = c("n_obs", "rho")) %>%
    relative(c("dgp_type", "p", "dgp_r2", "coefficient_type", "zero_status")) %>%
    mutate(rel = median / bench, rel_se = se / bench,
           coefficient_type = factor(coefficient_type, c("interaction", "main"))) %>%
    filter(is.finite(rel), is.finite(rel_se))
  save_figure(plot_grid(rmse_rel, rmse_block, sprintf("RMSE / RMSE(%s)", bench_label), TRUE),
              file.path("output/figures", paste0("rmse_", suffix)), 10, 7.6, 11, 7)

  elpd_rel <- elpd %>%
    semi_join(s, by = c("n_obs", "rho")) %>%
    relative(c("dgp_type", "p", "dgp_r2")) %>%
    mutate(rel = n_obs * (median - bench), rel_se = n_obs * se) %>%
    filter(is.finite(rel))
  save_figure(plot_grid(elpd_rel, elpd_block, sprintf("ELPD - ELPD(%s)", bench_label), FALSE),
              file.path("output/figures", paste0("elpd_", suffix)), 10, 6, 11, 5.5)
}
