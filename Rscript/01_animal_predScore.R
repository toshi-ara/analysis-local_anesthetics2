#!/usr/bin/env Rscript

########################################
## analysis duration by raw data
########################################

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
source("common_functions.R", local = TRUE)


dir.create("../results/01", recursive = TRUE, showWarnings = FALSE)
dir.create("rda", recursive = TRUE, showWarnings = FALSE)


########################################
# Set parameters of each individual
########################################

Time_pred <- seq(0, 150, by = 1)


## Read data & Shaping
filename <- "../data/param1_individual.csv"

dat_param1 <- read_csv(filename, show_col_types = FALSE) |>
    mutate(ID = factor(ID),
           drug = factor(drug, levels = c("Pro", "Lid", "Mep", "Bup")))

result_predProb_animal1 <- dat_param1 |>
    mutate(predProb = pmap(
               list(Mean, Sigma),
               function(time, Mean, Sigma, adr) {
                   data.frame(
                       time = time,
                       prob = predProb(time, Mean, Sigma, adr)
               )},
               time = Time_pred, adr = 0)) |>
    dplyr::select(drug, ID, predProb) |>
    unnest(cols = predProb)

save(
    dat_param1, result_predProb_animal1,
    file = "rda/animal_predProb.RData"
)

