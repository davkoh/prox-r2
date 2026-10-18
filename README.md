# Prox-R2

This repository contains the simulation study of the paper "Proximal R2D2 Priors for Structured Variance Allocation" by J. E. Aguilar, D. Kohns, E. Yanchenko, S. Sugasawa, K. Irie, A. Vehtari and P.-C. Bürkner.

The study compares four priors for linear regression with  pairwise interactions. It measures the RMSE of the coefficients and the expected log predictive density (ELPD) on test data.

## Models

| Key | Prior | Stan file |
|---|---|---|
| `prox_r2` | Prox-R2, the proposed prior | `stan/prox_r2.stan` |
| `r2d2` | R2-D2 (Aguilar and Bürkner, 2023) | `stan/r2d2.stan` |
| `normal_gamma` | Normal-gamma-gamma with heredity (Griffin and Brown, 2017) | `stan/normal_gamma.stan` |
| `rhs` | Regularized horseshoe (Piironen and Vehtari, 2017) | `stan/rhs.stan` |

## Contents

| Path | Contents |
|---|---|
| `Snakefile` | Workflow |
| `config/config.yaml` | Settings of the study in the paper |
| `config/test.yaml` | Small settings for a test run |
| `stan/` | Stan models |
| `R/` | R scripts for the workflow steps |
| `slurm/` | Snakemake profile and job script for SLURM clusters |
| `environment.yml` | Conda environment for Snakemake |

## Requirements

- R with the packages cmdstanr, posterior, mvtnorm, truncnorm, dplyr, tidyr, readr, ggplot2 and patchwork. The package tikzDevice is optional.
- CmdStan. To install it, use `cmdstanr::install_cmdstan()`.
- Snakemake 8. To install it, use `conda env create -f environment.yml`.

The results in the paper used R 4.3.0, CmdStan 2.35.0, cmdstanr 0.8.0, posterior 1.6.1 and Snakemake 8.18.0. On other hardware or software versions, the results can differ within the Monte Carlo error.

## Run the study

Run all commands from the root folder of the repository.

1. Do a test run with small settings:
   ```
   snakemake --cores 4 --configfile config/test.yaml
   ```
2. Remove the output of the test run: `rm -r output`.
3. Run the full study:
   ```
   snakemake --cores 4
   ```

The full study has 48 scenarios, 4 models and 100 replicates for each scenario. Use a cluster for the full study.

To run the study on a SLURM cluster:

1. Set the partition in `slurm/config.yaml`.
2. Make R, CmdStan and Snakemake available in `slurm/submit.sh`.
3. Submit the workflow: `sbatch slurm/submit.sh`.

## Simulation design

- The predictors have an AR(1) correlation `rho`. The design matrix has the `p` main effects and all `p(p-1)/2` pairwise products.
- The main effects are sparse (0.5 for the first three and the last three, 0 for the others) or dense (0.05 for all).
- The interactions satisfy weak heredity (`*_correct`) or violate it (`*_misp`).
- The residual standard deviation gives the target DGP-R2.
- Each replicate has `n_obs` training observations and `n_obs` test observations.

`config/config.yaml` contains all settings, the prior hyperparameters and the MCMC settings.

## Workflow

| Rule | Script | Output |
|---|---|---|
| `compile` | `R/compile_model.R` | `output/bin/<model>` |
| `generate_data` | `R/generate_data.R`, `R/dgp.R` | `output/data/<scenario>.RDS` |
| `fit` | `R/fit_model.R` | `output/fits/<scenario>/<model>_batch_<b>.csv` |
| `collect` | `R/collect_fits.R` | `output/fits/<scenario>.RDS` |
| `summarise` | `R/summarise.R` | `output/summaries/rmse.RDS`, `output/summaries/elpd.RDS` |
| `figures` | `R/plot_figures.R` | `output/figures/` |
| `tables` | `R/make_tables.R` | `output/tables/` |

A permutation of any dgp has the name `<dgp>_nobs_<n>_p_<p>_rho_<rho>_dgpr2_<R2>`.

The fit step keeps the posterior summaries of `beta`, `sigma`, `log_lik_test` and `y_tilde_test`. The figures use only the replicates where all models converge (mean R-hat less than 1.1). The tables use all replicates.

## Output

| File | Contents |
|---|---|
| `output/figures/rmse_rel_<benchmark>_nobs_<n>_rho_<rho>.*` | Median RMSE of the coefficients, relative to the benchmark model |
| `output/figures/elpd_rel_<benchmark>_nobs_<n>_rho_<rho>.*` | Difference in ELPD to the benchmark model |
| `output/tables/diagnostics_beta*.tex` | R-hat and bulk ESS of the coefficients |
| `output/tables/coverage_beta_p<p>*.tex` | Coverage of the central 90% posterior intervals of the coefficients |

The figures are PNG, PDF and TikZ (`.tex`) files. The TikZ files need tikzDevice and a LaTeX installation.

## Citation

If you use this code, cite the paper (link to arXiv will arrive soon).
