#!/usr/bin/env Rscript

library(dplyr)
library(xtable)


dir.create("Table", recursive = TRUE, showWarnings = FALSE)
source("common_functions.R", local = TRUE)
load("rda/animal_correlation_mean_logSigma.RData")


colnames(res_cor_mean_logSigma) <- c(
    "Drug",
    "all data ($n = 51$)",
    "without outliers ($n = 49$)"
)

print(file = "../results/Table/Table1.tex",
    xtable(res_cor_mean_logSigma, digits = c(0, 0, 3, 3)),
    math.style.negative = TRUE,
    include.rownames = FALSE,
    only.contents = TRUE,
    booktabs = TRUE,
    timestamp = NULL,
    sanitize.text.function = identity
)

