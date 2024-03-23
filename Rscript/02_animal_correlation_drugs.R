#!/usr/bin/env Rscript

########################################
## analysis colleration of estimated parameters
########################################

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(ggplot2)
library(GGally)
source("common_functions.R", local = TRUE)


dir.create("../results/02", recursive = TRUE, showWarnings = FALSE)
dir.create("rda", recursive = TRUE, showWarnings = FALSE)


## read parameters estimated by model2
dat_raw <- read_csv("../data/param2_individual.csv", show_col_types = FALSE)

## dat1: for drug parameter
## dat2: for drug parameter (eliminate outliers)


########################################
## Condition 1
########################################

dat1 <- dat_raw |>
    mutate(ID = factor(ID),
           drug = factor(drug, levels = c("Pro", "Lid", "Mep", "Bup")),
           logSigma = log(Sigma)) |>
    group_by(drug) |>
    mutate(Mean_score = scale(Mean)[, 1],
           logSigma_score = scale(logSigma)[, 1])


########################################
## Condition 2
## eliminate outlier (scale > 2)
########################################

threshold <- 2
outlier_id <- dat1 |>
    dplyr::filter(abs(Mean_score) > threshold) |>
    pull(ID) |>
    unique()

dat1 <- dat1 |>
    mutate(not_outlier = !(ID %in% outlier_id))

dat2 <- dat1 |>
    dplyr::filter(not_outlier) |>
    group_by(drug) |>
    mutate(Mean_score = scale(Mean_score)[, 1],
           logSigma_score = scale(logSigma_score)[, 1],
           not_outlier = not_outlier)

dat_param <- list(dat1, dat2)


dat_param_drug <- lapply(dat_param, function(.x) {
    .x |>
        pivot_longer(cols = c(Mean, logSigma, Mean_score, logSigma_score),
                     names_to = c("parameter", "score"),
                     names_sep = "_",
                     values_to = "value",
        ) |>
        mutate(score = !is.na(score),
               parameter = factor(parameter,
                                  levels = c("Mean", "logSigma"))
        )
})

dat_param_drug_w <- lapply(dat_param_drug, function(.x) {
    .x |>
        pivot_wider(id_cols = c(score, parameter, ID),
                    names_from = drug,
                    values_from = value) |>
        mutate(not_outlier = !(ID %in% outlier_id))
})


#########################################
### Correlation of parameters (Mean and logSigma) between Mean and logSigma
#########################################

res_cor_mean_logSigma <- bind_rows(dat_param, .id = "dat") |>
    group_by(dat, drug) |>
    summarise(r = cor(Mean_score, logSigma_score)) |>
    pivot_wider(names_from = dat, values_from = r) |>
    rename(with_outlier = `1`,
           without_outlier = `2`)
# # A tibble: 4 × 3
#   drug  with_outlier without_outlier
#   <fct>        <dbl>           <dbl>
# 1 Pro        -0.308          -0.219 
# 2 Lid        -0.415          -0.301 
# 3 Mep         0.0116          0.0140
# 4 Bup        -0.154          -0.160 

write_csv(res_cor_mean_logSigma,
          file = "../results/02/cor_mean_logSigma.csv")



########################################
## for correlation coefficient / scatter plot
########################################

res_cor <- lapply(dat_param_drug_w, function(.x) {
    .x |>
        dplyr::filter(!score) |>
        group_nest(parameter) |>
        mutate(r = map(data, ~
                       cor(.[c("Pro", "Lid", "Mep", "Bup")]))) |>
        dplyr::select(-data) |>
        pivot_wider(names_from = parameter, values_from = r)
})
# [[1]]
# # A tibble: 1 × 2
#   Mean          logSigma     
#   <list>        <list>       
# 1 <dbl [4 × 4]> <dbl [4 × 4]>

# [[2]]
# # A tibble: 1 × 2
#   Mean          logSigma     
#   <list>        <list>       
# 1 <dbl [4 × 4]> <dbl [4 × 4]>



########################################
## save results
########################################

save(
    dat_raw, dat_param,
    dat_param_drug, dat_param_drug_w,
    threshold, outlier_id,
    res_cor_mean_logSigma, res_cor,
    file = "rda/animal_correlation_mean_logSigma.RData"
)



#########################################
### Correlation of duration between Mean and logSigma
#########################################

dat_score <- read_csv("../data/local_anesteshia_data.csv",
                      show_col_types = FALSE) |>
    pivot_longer(col = c(-ID, -Drug),
                 names_to = "time",
                 names_transform = list(time = as.numeric),
                 values_to = "score",
                 values_drop_na = TRUE) |>
    mutate(Drug = factor(Drug, levels = drug_name),
           drug1 = fct_collapse(Drug,
                                Pro = "Pro", Lid = c("Lid", "Lid+Adr")),
           drug2 = ifelse(Drug == "Lid+Adr", 1L, 0L))

duration <- groupSurvData(dat_score, k = 6)

duration_w <- duration |>
    pivot_wider(id_cols = ID, names_from = Drug, values_from = time)


########################################
## save results
########################################

save(
    dat_score,
    duration, duration_w,
    file = "rda/animal_duration_correlation.RData"
)

