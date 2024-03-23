#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(purrr)
library(xtable)


########################################
## Coefficients of parameters among drugs for simulation
########################################

dir.create("Table", recursive = TRUE, showWarnings = FALSE)
load("rda/animal_predProb.RData")
load("rda/animal_duration_correlation.RData")
load("rda/sim_correlation_drugs.RData")

## function
judge_Pro_munus_Lid_param <- function(x) {
    res <- ifelse(x > 0, "Pro > Lid", "Pro < Lid")
    factor(res, levels = c("Pro > Lid", "Pro < Lid"))
}

judge_Lid_munus_Pro <- function(x) {
    res <- ifelse(is.nan(x), "both Inf",
           ifelse(x > 0, "Pro < Lid",
                  ifelse(x == 0, "Pro = Lid",
                         "Pro > Lid")))
    factor(res, levels = c("Pro < Lid", "Pro = Lid", "Pro > Lid", "both Inf"))
}


####################
## paramter
####################

## animal experiment
animal_param_diff <- dat_param1 |>
    pivot_wider(id_cols = ID,
                names_from = "drug",
                values_from = "Mean") |>
    dplyr::select(ID, Pro, Lid) |>
    mutate(judge = judge_Pro_munus_Lid_param(Pro - Lid),
           total = n()) |>
    group_by(judge) |>
    summarise(n = n(),
              rate = n / total[1],
              n_rate = sprintf("%d (%.1f\\%%)", n, rate * 100)) |>
    dplyr::select(judge, n_rate)

## simulation
sim_param_diff <- result_sim_duration |>
    mutate(param = map(param_drug_sim, ~
                       pivot_wider(., id_cols = ID,
                                   names_from = "Drug",
                                   values_from = "Mean") |>
                       mutate(judge = judge_Pro_munus_Lid_param(Pro - Lid))
                   )
    ) |>
    unnest(cols = param) |>
    dplyr::select(seed_param, seed_sim, condition, Pro, Lid, judge) |>
    ## summarise
    group_nest(seed_param, seed_sim, condition) |>
    mutate(res = map(data, ~ mutate(., total = n())),
           judge = map(res, ~
                       group_by(., judge) |>
                           summarise(n = n(),
                                     rate = n / total[1],
                                     n_rate = sprintf("%d (%.1f\\%%)",
                                                      n, rate * 100)) |>
                       dplyr::select(judge, n_rate)
                   )) |>
    dplyr::select(-data, -res)

## merge
res_param_diff_tbl <- sim_param_diff |>
    group_nest(seed_param, seed_sim) |>
    mutate(res = map(data, ~
                     ## join animal / simulation data
                     left_join(
                         animal_param_diff,
                         pivot_wider(unnest(., cols = judge),
                                     names_from = "condition",
                                     values_from = "n_rate",
                                     names_prefix = "Condition"),
                         by = "judge"))) |>
    dplyr::select(-data) |>
    unnest(cols = res) |>
    ## change row / col names
    rename(raw = n_rate) |>
    mutate(judge = factor(judge,
                          labels = c("Pro $>$ Lid", "Pro $<$ Lid")))

res_param_diff_latex <- res_param_diff_tbl |>
    rename(
        Comparison = judge,
        `Raw data` = raw,
        `Condition 1` = Condition1,
        `Condition 2` = Condition2,
        `Condition 3` = Condition3,
        `Condition 4` = Condition4
    ) |>
    group_nest(seed_param, seed_sim, .key = "table")


####################
## duration
####################

## animal experiment
animal_duration_diff <- duration_w |>
    dplyr::select(ID, Pro, Lid) |>
    mutate(diff = Lid - Pro,
           judge = judge_Lid_munus_Pro(diff),
           total = n()) |>
    group_by(judge) |>
    summarise(n = n(),
              rate = n / total[1],
              n_rate = sprintf("%d (%.1f\\%%)", n, rate * 100)) |>
    dplyr::select(judge, n_rate)

## simulation
sim_duration_diff <- result_sim_duration_w |>
    mutate(diff = Lid - Pro,
           judge = judge_Lid_munus_Pro(diff)) |>
    group_nest(seed_param, seed_sim, condition) |>
    mutate(res = map(data, ~ mutate(., total = n())),
           judge = map(res, ~
                       group_by(., judge) |>
                           summarise(n = n(),
                                     rate = n / total[1],
                                     n_rate = sprintf("%d (%.1f\\%%)",
                                                      n, rate * 100)) |>
                       dplyr::select(judge, n_rate)
                   )) |>
    dplyr::select(-data, -res)

## merge
res_diff_tbl <- sim_duration_diff |>
    group_nest(seed_param, seed_sim) |>
    mutate(res = map(data, ~
                     ## join animal / simulation data
                     left_join(
                         animal_duration_diff,
                         pivot_wider(unnest(., cols = judge),
                                     names_from = "condition",
                                     values_from = "n_rate",
                                     names_prefix = "Condition"),
                         by = "judge"))) |>
    dplyr::select(-data) |>
    unnest(cols = res) |>
    ## change row / col names
    rename(raw = n_rate) |>
    mutate(judge = factor(judge,
                          labels = c("Pro $<$ Lid", "Pro $=$ Lid",
                                     "Pro $>$ Lid", "both censored")))

res_diff_latex <- res_diff_tbl |>
    rename(
        Comparison = judge,
        `Raw data` = raw,
        `Condition 1` = Condition1,
        `Condition 2` = Condition2,
        `Condition 3` = Condition3,
        `Condition 4` = Condition4
    ) |>
    group_nest(seed_param, seed_sim, .key = "table")


####################
## merge parameter / duration
####################

tbl1_param <- res_param_diff_latex |>
    dplyr::filter(seed_param == 1 & seed_sim == 4)

tbl1_duration <- res_diff_latex |>
    dplyr::filter(seed_param == 1 & seed_sim == 4)


tbl_param <- tibble(
    ` ` = c("Parameter", " "),
    tbl1_param$table[[1]]
)

tbl_duration <- tibble(
    ` ` = c("Duration", rep(" ", 3)),
    tbl1_duration$table[[1]]
)

tbl1 <- bind_rows(tbl_param, tbl_duration)


## output
print(file = "../results/Table/Table6.tex",
    xtable(tbl1),
    # size = "small",
    NA.string = "---",
    include.rownames = FALSE,
    only.contents = TRUE,
    hline.after = c(0, 2, 6),
    booktabs = TRUE,
    timestamp = NULL,
    sanitize.text.function = identity
)

