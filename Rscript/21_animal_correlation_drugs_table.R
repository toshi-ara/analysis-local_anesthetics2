#!/usr/bin/env Rscript

library(dplyr)
library(xtable)


########################################
## combine (PCP / Mean / logSigma)
########################################

dir.create("Table", recursive = TRUE, showWarnings = FALSE)
source("common_functions.R", local = TRUE)
load("rda/animal_correlation_mean_logSigma.RData")
load("rda/animal_duration_correlation.RData")



########################################
## for correlation coefficient
########################################

####################
## correlation coefficient of parameters among drugs
####################

## function
get_cor_df <- function(.x) {
    mat_mean <- .x$Mean[[1]]
    mat_logSigma <- .x$logSigma[[1]]

    ## combination of names
    name_row <- rownames(mat_mean)
    combination <- outer(name_row, name_row,
                         function(.x, .y) {
                             paste(.y, .x, sep = "--")
                         })

    ## only upper tri
    diag(mat_mean) <- NA
    mat_mean[upper.tri(mat_mean)] <- NA
    diag(mat_logSigma) <- NA
    mat_logSigma[upper.tri(mat_logSigma)] <- NA

    ## drop NA
    dat <- data.frame(
        combination = as.vector(combination),
        Mean = as.vector(mat_mean),
        logSigma = as.vector(mat_logSigma)
    )
    return(dat[complete.cases(dat), ])
}

res_cor_df <- lapply(res_cor, get_cor_df)
res_cor_df_combine <- left_join(
    res_cor_df[[1]], res_cor_df[[2]], by = "combination"
)

colnames(res_cor_df_combine) <- c("Combination",
                 "$r_{\\mu}$",
                 "$r_{\\log \\sigma}$",
                 "$r_{\\mu}$",
                 "$r_{\\log \\sigma}$"
)

print(file = "../results/Table/Table2.tex",
    xtable(res_cor_df_combine, digits = 3),
    include.rownames = FALSE,
    only.contents = TRUE,
    booktabs = TRUE,
    timestamp = NULL,
    sanitize.text.function = identity
)



####################
## correlation coefficient of duration among drugs
####################

## function
get_cor_df2 <- function(x, y) {
    ## combination of names
    name_row <- rownames(x)
    combination <- outer(name_row, name_row,
                         function(.x, .y) {
                             paste(.y, .x, sep = "--")
                         })

    ## only upper tri
    diag(x) <- NA
    x[upper.tri(x)] <- NA
    diag(y) <- NA
    y[upper.tri(y)] <- NA

    ## drop NA
    dat <- data.frame(
        combination = as.vector(combination),
        x = as.vector(x),
        y = as.vector(y)
    )
    return(dat[complete.cases(dat), ])
}

## compare among drugs
mat1 <- duration_w |>
    dplyr::select(-ID, -`Lid+Adr`) |>
    cor(method = "spearman")

mat2 <- duration_w |>
    dplyr::filter(!(ID %in% outlier_id)) |>
    dplyr::select(-ID, -`Lid+Adr`) |>
    cor(method = "spearman")

res_cor_df2 <- get_cor_df2(mat1, mat2)
colnames(res_cor_df2) <- c(
    "Combination",
    "all data ($n = 51$)",
    "without outliers ($n = 49$)"
)

print(file = "../results/Table/STable1.tex",
    xtable(res_cor_df2, digits = 3),
    math.style.negative = TRUE,
    include.rownames = FALSE,
    only.contents = TRUE,
    booktabs = TRUE,
    timestamp = NULL,
    sanitize.text.function = identity
)

