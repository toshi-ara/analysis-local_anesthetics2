#!/usr/bin/env Rscript

########################################
## simulation the correlation between Mean and logSigma
########################################

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
source("common_functions.R")

dir.create("../results/03", recursive = TRUE, showWarnings = FALSE)


########################################
# Set parameters of each individual
########################################

r <- seq(-1, 1, by = 0.25)
Time <- seq(0, 180, by = 5)
Time_pred <- seq(0, 120, by = 1)
N <- 100

seed_param <- 1:8
seed_sim <- 1:10



########################################
# Get predict score from probability curve
########################################

condition <- expand_grid(
    seed_param = seed_param,
    r = r
)

result_sim_prob <- condition |>
    mutate(## set parameters
           param_drug_sim = pmap(
               list(r = r, seed = seed_param),
               set_parameters_one,
               N = N,
               parameters = param_drug[2, ],
               d_sigma = d_sigma,
               drug_name = "Lid")) |>
    unnest(cols = param_drug_sim) |>
    mutate(predProb = pmap(
               list(Mean, Sigma, adr),
               function(time, Mean, Sigma, adr) {
                   data.frame(
                       time = time,
                       prob = predProb(time, Mean, Sigma, adr)
               )},
               time = Time_pred)) |>
    dplyr::select(seed_param, r, Drug, ID, predProb) |>
    unnest(cols = predProb)



########################################
# Get duration
########################################

condition_duration <- expand_grid(
    seed_param = seed_param,
    seed_sim = seed_sim,
    r = r
)

param_duration <- condition_duration |>
    mutate(## set parameters
           param_drug_sim = pmap(
               list(r = r, seed = seed_param),
               set_parameters_one,
               N = N,
               parameters = param_drug[2, ],
               d_sigma = d_sigma,
               drug_name = "Lid"))


## violin plot + points
result_sim_duration <- param_duration |>
    mutate(## simulation
           res = pmap(
               list(x = param_drug_sim, seed = seed_sim),
               doSimulation,
               time = Time, n_stim = 6),
           duration = map(res, ~ groupSurvData(., k = 6))
    ) |>
    dplyr::select(-param_drug_sim, -res) |>
    unnest(cols = duration)


## summary (mean, SD)
result_sim_duration_summary <- result_sim_duration |>
    group_by(seed_param, seed_sim, r) |>
    summarise(mean = mean(time),
              SD = sd(time),
              median = median(time),
              Q1 = quantile(time, prob = 0.25),
              Q3 = quantile(time, prob = 0.75),
              IQR = Q3 - Q1) |>
    ungroup()

write_csv(result_sim_duration_summary,
          file = "../results/03/sim_duration_summary.csv")



########################################
## Survival analysis of simulation data
##   correlation of Mean and logSigma in various r
########################################

## function
get_surv_df <- function(.x) {
    survival:::survmean(.x, rmean = "none")$matrix |>
        data.frame() |>
        tibble::rownames_to_column() |>
        dplyr::select(!!sym("rowname"),
                      !!sym("n.start"),
                      !!sym("events"), !!sym("median"),
                      starts_with("X0.95")) |>
        separate(col = !!sym("rowname"),
                 into = c(NA, "r"), sep = "=") |>
        mutate(r = as.numeric(r)) |>
        rename(n = !!sym("n.start"),
               lwr = !!sym("X0.95LCL"), upr = !!sym("X0.95UCL"))
}

result_sim_survival <- result_sim_duration |>
    group_nest(seed_param, seed_sim) |>
    mutate(surv = map(data, ~
                      survfit(Surv(time, status) ~ r, data = .)))

result_sim_survival_tbl <- result_sim_survival |>
    mutate(surv_tbl = map(surv, ~ get_surv_df(.))) |>
    dplyr::select(seed_param, seed_sim, surv_tbl) |>
    unnest(cols = surv_tbl)

write_csv(result_sim_survival_tbl,
          file = "../results/03/sim_duration_survival.csv")



########################################
## save results
########################################

save(
    result_sim_prob,
    result_sim_duration,
    result_sim_duration_summary,
    result_sim_survival,
    result_sim_survival_tbl,
    file = "rda/sim_correlation_Mean_logSigma.RData"
)

