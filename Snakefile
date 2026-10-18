configfile: "config/config.yaml"

MODELS = list(config["models"])
BENCH = config["benchmark_model"]
SCENARIO = "{dgp}_nobs_{n_obs}_p_{p}_rho_{rho}_dgpr2_{r2}"
SCENARIOS = expand(SCENARIO, dgp=config["dgp_types"], n_obs=config["n_obs"],
                   p=config["p"], rho=config["rho"], r2=config["dgp_r2"])
SETTINGS = expand("nobs_{n_obs}_rho_{rho}", n_obs=config["n_obs"], rho=config["rho"])
FITS = expand("output/fits/{scenario}.RDS", scenario=SCENARIOS)

wildcard_constraints:
    dgp=r"[a-z_]+", n_obs=r"\d+", p=r"\d+", rho=r"[\d.]+", r2=r"[\d.]+",
    model=r"\w+", batch=r"\d+"


rule all:
    input:
        expand("output/figures/{metric}_rel_{bench}_{setting}.{ext}",
               metric=["rmse", "elpd"], bench=BENCH, setting=SETTINGS, ext=["png", "pdf"]),
        expand("output/tables/diagnostics_beta_{setting}.tex", setting=SETTINGS),
        expand("output/tables/{table}_beta_p{p}_{setting}.tex",
               table=["diagnostics", "coverage"], p=config["p"], setting=SETTINGS)


rule compile:
    input: "stan/{model}.stan"
    output: "output/bin/{model}"
    script: "R/compile_model.R"


rule generate_data:
    output: f"output/data/{SCENARIO}.RDS"
    script: "R/generate_data.R"


rule fit:
    input:
        data=f"output/data/{SCENARIO}.RDS",
        exe="output/bin/{model}"
    output: f"output/fits/{SCENARIO}/{{model}}_batch_{{batch}}.csv"
    threads: config["mcmc"]["chains"]
    script: "R/fit_model.R"


rule collect:
    input:
        data=f"output/data/{SCENARIO}.RDS",
        csv=expand(f"output/fits/{SCENARIO}/{{model}}_batch_{{batch}}.csv", model=MODELS,
                   batch=range(1, config["n_batches"] + 1), allow_missing=True)
    output: f"output/fits/{SCENARIO}.RDS"
    script: "R/collect_fits.R"


rule summarise:
    input: FITS
    output:
        rmse="output/summaries/rmse.RDS",
        elpd="output/summaries/elpd.RDS"
    script: "R/summarise.R"


rule figures:
    input:
        rmse="output/summaries/rmse.RDS",
        elpd="output/summaries/elpd.RDS"
    output:
        expand("output/figures/{metric}_rel_{bench}_{setting}.{ext}",
               metric=["rmse", "elpd"], bench=BENCH, setting=SETTINGS, ext=["png", "pdf"])
    script: "R/plot_figures.R"


rule tables:
    input: FITS
    output:
        expand("output/tables/diagnostics_beta_{setting}.tex", setting=SETTINGS),
        expand("output/tables/{table}_beta_p{p}_{setting}.tex",
               table=["diagnostics", "coverage"], p=config["p"], setting=SETTINGS)
    script: "R/make_tables.R"
