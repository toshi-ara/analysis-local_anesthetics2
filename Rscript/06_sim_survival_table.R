#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
source("common_functions.R", local = TRUE)

dir.create("../results/06", recursive = TRUE, showWarnings = FALSE)


########################################
## survival analysys in various r between mean and logSigma
########################################

dir.create("../results/Fig", recursive = TRUE, showWarnings = FALSE)
dir.create("../results/Table", recursive = TRUE, showWarnings = FALSE)
load("rda/animal_duration_correlation.RData")
load("rda/sim_correlation_drugs.RData")


lab_condition <- c("Raw data", paste("Condition", 1:4))

## function
get_surv_df <- function(.x) {
    survival:::survmean(.x, rmean = "none")$matrix |>
        data.frame() |>
        tibble::rownames_to_column() |>
        dplyr::select(!!sym("rowname"),
                      !!sym("n.start"),
                      !!sym("events"), !!sym("median"),
                      starts_with("X0.95")) |>
        separate_wider_delim(col = !!sym("rowname"),
                             delim = ",",
                             names = c("Drug", "Condition")) |>
        separate_wider_delim(col = !!sym("Drug"),
                             delim = "=",
                             names = c(NA, "Drug")) |>
        separate_wider_delim(col = !!sym("Condition"),
                             delim = "=",
                             names = c(NA, "Condition")) |>
        rename(n = !!sym("n.start"),
               lwr = !!sym("X0.95LCL"), upr = !!sym("X0.95UCL")) |>
        mutate(Drug = factor(Drug, levels = drug_name),
               Condition = Condition |>
                   stringr::str_trim() |>
                   factor(levels = lab_condition, labels = lab_condition))
}


survdata_animal <- duration |>
    mutate(condition = 0) |>
    dplyr::select(condition, everything())

result_survival <- result_sim_duration |>
    dplyr::select(seed_param, seed_sim, condition, duration) |>
    group_nest(seed_param, seed_sim) |>
    ## set data for survival analysis
    mutate(survdata = map(data, ~
                          unnest(., cols = duration)),
           survdata = map(survdata, ~
                          bind_rows(survdata_animal, .)),
           survdata = map(survdata, ~
                          mutate(., condition = factor(condition,
                                                       labels = lab_condition)))
    ) |>
    ## preform survival analysis
    mutate(res_surv = map(survdata, ~
                     survfit(Surv(time, status) ~ Drug + condition,
                             data = .)),
           res = map(res_surv, ~ get_surv_df(.))) |>
    dplyr::select(-data)



#################
## save as CSV
#################
result_survival_tbl <- result_survival |>
    dplyr::select(seed_param, seed_sim, res) |>
    unnest(cols = res)

write_csv(result_survival_tbl,
          file = "../results/06/animal_sim_survival.csv")



########################################
## save results
########################################

save(
    lab_condition,
    result_survival, result_survival_tbl,
    file = "rda/animal_sim_survival.RData"
)

