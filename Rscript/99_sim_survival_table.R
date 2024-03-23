#!/usr/bin/env Rscript

library(dplyr)
library(xtable)


########################################
## survival analysys in various r between mean and logSigma
########################################

dir.create("Table", recursive = TRUE, showWarnings = FALSE)
load("rda/sim_correlation_Mean_logSigma.RData")



tbl1 <- result_sim_survival_tbl |>
    dplyr::filter(seed_param == 1 & seed_sim == 1) |>
    dplyr::select(-seed_param, -seed_sim)

colnames(tbl1) <- c(
    "$r$",
    "$n$",
    "events",
    "median",
    "lower",
    "upper"
)

print(file = "../results/Table/tbl_sim_r_survival.tex",
    xtable(tbl1, digits = c(0, 2, 0, 0, 1, 1, 1)),
    math.style.negative = TRUE,
    include.rownames = FALSE,
    only.contents = TRUE,
    booktabs = TRUE,
    timestamp = NULL,
    sanitize.text.function = identity
)

