#!/usr/bin/env Rscript

library(dplyr)
library(xtable)


dir.create("Table", recursive = TRUE, showWarnings = FALSE)
source("common_functions.R", local = TRUE)
load("rda/sim_correlation_drugs.RData")


## functions
get_param_from_sigma_matrix <- function(x) {
    n <- nrow(x)

    Mean <- x[1:(n / 2), 1:(n / 2)] |>
        get_lower_tri_vector()
    logSigma <- x[(n / 2 + 1):n, (n / 2 + 1):n] |>
        get_lower_tri_vector()
    r <- x[1:(n / 2), (n / 2 + 1):n] |>
        diag()

    return(list(
        mean_logSigma = data.frame(
            mean = Mean,
            logSigma = logSigma
        ),
        cor_mean_logSigma = data.frame(
            r = r
        )
    ))
}


combination <- outer(drug_name[1:4], drug_name[1:4],
                     function(.x, .y) {
                         paste(.y, .x, sep = "--")
                     }) |>
    get_lower_tri_vector()

sim_params <- get_param_from_sigma_matrix(sigma_matrix[[4]])


sim_params_tbl <- tibble(
    drug = drug_name[1:4],
    param_drug |>
        dplyr::select(-Drug, -adr),
    r_mu_logSigma = sim_params$cor_mean_logSigma$r
)
sim_params_tbl[5:6, ] <- NA

sim_params_tbl <- sim_params_tbl |>
    bind_cols(
        tibble(
            combination = combination,
            sim_params$mean_logSigma
        )
    )


colnames(sim_params_tbl) <- c(
    "Drug",
    "$\\mu_{0}$",
    "$s_{\\mu_{0}}$",
    "$\\log \\sigma_{0}$",
    "$s_{\\log \\sigma_{0}}$",
    "$r_{\\mu - \\log \\sigma}$",
    "Combination",
    "$r_{\\mu}$",
    "$r_{\\log \\sigma}$"
)

print(file = "../results/Table/Table3.tex",
    xtable(sim_params_tbl, digits = c(0, 0, 0, 0, 1, 1, 2, 0, 2, 2)),
    include.rownames = FALSE,
    only.contents = TRUE,
    booktabs = TRUE,
    timestamp = NULL,
    sanitize.text.function = identity
)

