#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(survminer)
library(ggplot2)
library(GGally)
library(ggpubr)
library(patchwork)
source("common_functions.R")

dir.create("../results/04", recursive = TRUE, showWarnings = FALSE)


########################################
# Set parameters of each individual
########################################

Time <- seq(0, 180, by = 5)
N <- 200


seed_param <- 1:8
condition <- expand_grid(
    seed_param = seed_param,
    r_mean = seq(0, 1, by = 0.2),
    r_logSigma = seq(0, 1, by = 0.2)
)
# # A tibble: 288 × 3
#    seed_param r_mean r_logSigma
#         <int>  <dbl>      <dbl>
#  1          1    0          0  
#  2          1    0          0.2
#  3          1    0          0.4
#  4          1    0          0.6
#  5          1    0          0.8
#  6          1    0          1  
#  7          1    0.2        0  
#  8          1    0.2        0.2
#  9          1    0.2        0.4
# 10          1    0.2        0.6
# # ℹ 278 more rows

seed_sim <- 1234
set.seed(seed_sim)
result_sim_Pro_Lid <- condition |>
    mutate(param_drug = pmap(list(r_mean = r_mean,
                                  r_logSigma = r_logSigma,
                                  seed = seed_param),
                             set_parameters_r_Pro_Lid,
                             N = N,
                             parameters = param_drug[1:2, ],
                             d_sigma = d_sigma,
                             drug_name = c("Pro", "Lid"))) |>
    mutate(res = map(param_drug, ~ doSimulation(., Time, n_stim = 6)),
           duration = map(res, ~ groupSurvData(., k = 6)),
           duration_w = map(duration, ~
               pivot_wider(., id_cols = ID,
                           names_from = Drug, values_from = time)),
           r = map_dbl(duration_w, ~ cor(.$Pro, .$Lid)))

result_sim_duration_Pro_Lid <- result_sim_Pro_Lid |>
    dplyr::select(seed_param, r_mean, r_logSigma, duration_w) |>
    unnest(cols = duration_w)

write_csv(result_sim_duration_Pro_Lid,
          file = "../results/04/sim_duration_Pro_Lid.csv")



########################################
## save results
########################################

save(
    N, seed_param, seed_sim,
    result_sim_Pro_Lid,
    result_sim_duration_Pro_Lid,
    file = "rda/sim_correlation_Pro_Lid.RData"
)

