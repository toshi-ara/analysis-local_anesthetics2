#!/usr/bin/env Rscript

library(dplyr)
library(xtable)


########################################
## Coefficients of parameters among drugs for simulation
########################################

dir.create("Table", recursive = TRUE, showWarnings = FALSE)
source("common_functions.R", local = TRUE)
load("rda/sim_correlation_drugs.RData")


## combination of two drugs
combination <- outer(drug_name[1:4], drug_name[1:4],
                     function(.x, .y) {
                         paste(.y, .x, sep = "--")
                     }) |>
    get_lower_tri_vector()


## correlatation of duration
res_cor <- lapply(1:4, function(i) {
    get_lower_tri_vector(result_sim_duration_cor$r[[i]])
})

res_cor_rho <- lapply(1:4, function(i) {
    get_lower_tri_vector(result_sim_duration_cor$rho[[i]])
})


duration_r <- data.frame(
    combination,
    res_cor[[1]],
    res_cor[[2]],
    res_cor[[3]],
    res_cor[[4]]
)

duration_rho <- data.frame(
    combination,
    res_cor_rho[[1]],
    res_cor_rho[[2]],
    res_cor_rho[[3]],
    res_cor_rho[[4]]
)

print(file = "../results/Table/Table4.tex",
    xtable(duration_r, digits = c(0, 0, 3, 3, 3, 3)),
    math.style.negative = TRUE,
    include.rownames = FALSE,
    include.colnames = FALSE,
    only.contents = TRUE,
    booktabs = TRUE,
    timestamp = NULL,
    sanitize.text.function = identity
)

print(file = "../results/Table/STable3.tex",
    xtable(duration_rho, digits = c(0, 0, 3, 3, 3, 3)),
    math.style.negative = TRUE,
    include.rownames = FALSE,
    include.colnames = FALSE,
    only.contents = TRUE,
    booktabs = TRUE,
    timestamp = NULL,
    sanitize.text.function = identity
)

