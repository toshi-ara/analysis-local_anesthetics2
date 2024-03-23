#!/usr/bin/env Rscript

########################################
## simulation the correlation among drugs
########################################

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
source("common_functions.R")

dir.create("../results/05", recursive = TRUE, showWarnings = FALSE)


########################################
# Set parameters of each individual
########################################

Time <- seq(0, 180, by = 5)
N <- 100


## r of Mean
# > cor(dat2_mean_score)
#           Pro       Lid       Mep       Bup
# Pro 1.0000000 0.5682454 0.4673379 0.4975916
# Lid 0.5682454 1.0000000 0.5264144 0.5270112
# Mep 0.4673379 0.5264144 1.0000000 0.4198695
# Bup 0.4975916 0.5270112 0.4198695 1.0000000

## r of logSigma
# > cor(dat2_logSigma_score)
#           Pro       Lid      Mep       Bup
# Pro 1.0000000 0.4156509 0.335883 0.2571197
# Lid 0.4156509 1.0000000 0.413687 0.4661620
# Mep 0.3358830 0.4136870 1.000000 0.5594290
# Bup 0.2571197 0.4661620 0.559429 1.0000000


## r between mean and logSigma
r1 <- -0.22  # Pro
r2 <- -0.30  # Lid
r3 <- -0.01  # Mep
r4 <- -0.16  # Bup

## mean among drugs
m12 <- 0.57  # Pro vs. Lid
m13 <- 0.47  # Pro vs. Mep
m14 <- 0.50  # Pro vs. Bup
m23 <- 0.53  # Lid vs. Mep
m24 <- 0.53  # Lid vs. Bup
m34 <- 0.42  # Mep vs. Bup

## logSigma among drugs
s12 <- 0.42  # Pro vs. Lid
s13 <- 0.34  # Pro vs. Mep
s14 <- 0.26  # Pro vs. Bup
s23 <- 0.41  # Lid vs. Mep
s24 <- 0.47  # Lid vs. Bup
s34 <- 0.56   # Mep vs. Bup


## sigma matrix

## Condition 1
sigma_matrix1 <- matrix(c(
    1, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 1
) , 8, 8)

## Condition 2
sigma_matrix2 <- matrix(c(
    1,  0,  0,  0,  r1, 0,  0,  0,
    0,  1,  0,  0,  0,  r2, 0,  0,
    0,  0,  1,  0,  0,  0,  r3, 0,
    0,  0,  0,  1,  0,  0,  0,  r4,
    r1, 0,  0,  0,  1,  0,  0,  0,
    0,  r2, 0,  0,  0,  1,  0,  0,
    0,  0,  r3, 0,  0,  0,  1,  0,
    0,  0,  0,  r4, 0,  0,  0,  1
) , 8, 8)

## Condition 3
sigma_matrix3 <- matrix(c(
    1,   m12, m13, m14, 0,   0,   0,   0,
    m12, 1,   m23, m24, 0,   0,   0,   0,
    m13, m23, 1,   m34, 0,   0,   0,   0,
    m14, m24, m34, 1,   0,   0,   0,   0,
    0,   0,   0,   0,   1,   s12, s13, s14,
    0,   0,   0,   0,   s12, 1,   s23, s24,
    0,   0,   0,   0,   s13, s23, 1,   s34,
    0,   0,   0,   0,   s14, s24, s34, 1
) , 8, 8)

## Condition 4
sigma_matrix4 <- matrix(c(
    1,   m12, m13, m14, r1,  0,   0,   0,
    m12, 1,   m23, m24, 0,   r2,  0,   0,
    m13, m23, 1,   m34, 0,   0,   r3,  0,
    m14, m24, m34, 1,   0,   0,   0,   r4,
    r1,  0,   0,   0,   1,   s12, s13, s14,
    0,   r2,  0,   0,   s12, 1,   s23, s24,
    0,   0,   r3,  0,   s13, s23, 1,   s34,
    0,   0,   0,   r4,  s14, s24, s34, 1
) , 8, 8)


sigma_matrix <- list(
    sigma_matrix1,
    sigma_matrix2,
    sigma_matrix3,
    sigma_matrix4
)


seed_param <- 1:6
seed_sim <- 1:4


condition <- expand_grid(
    seed_param = seed_param,
    seed_sim = seed_sim,
    condition = seq_len(length(sigma_matrix))
) |>
    mutate(r = sigma_matrix[condition])



########################################
# Start simulation
########################################

result_sim_duration <- condition |>
    mutate(## set parameters
           param_drug_sim = pmap(
               list(r = r, seed = seed_param),
               set_parameters_r,
               N = N,
               parameters = param_drug, d_sigma = d_sigma,
               drug_name = drug_name),
           ## simulation
           res = pmap(
               list(x = param_drug_sim, seed = seed_sim),
               doSimulation,
               time = Time, n_stim = 6),
           duration = map(res, ~ groupSurvData(., k = 6)),
           duration_w = map(duration, ~
                pivot_wider(., id_cols = ID,
                            names_from = Drug, values_from = time))
)

## save duration as csv file
result_sim_duration_w <- result_sim_duration |>
    dplyr::select(seed_param, seed_sim, condition, duration_w) |>
    unnest(cols = duration_w)

write_csv(result_sim_duration_w, file = "../results/05/sim_duration.csv")



####################
## correlation coefficients of duration among drugs
####################

result_sim_duration_cor <- result_sim_duration_w |>
    group_nest(seed_param, seed_sim, condition) |>
    mutate(df = map(data, ~
                    dplyr::select(., Pro:Bup)),
           r = map(df, ~ cor(.)),
           rho = map(df, ~ cor(., method = "spearman"))) |>
    dplyr::select(seed_param, seed_sim, condition, r, rho)



####################
## compare duration between Pro and Lid
####################

## function
compare_Pro_Lid <- function(x) {
    duration_diff <- x |>
        dplyr::select(!!sym("ID"), !!sym("Pro"), !!sym("Lid")) |>
        mutate(diff = !!sym("Lid") - !!sym("Pro"),
               judge = ifelse(is.nan(diff), "both Inf",
                              ifelse(diff > 0, "Pro < Lid",
                                     ifelse(diff == 0, "Pro = Lid",
                                            "Pro > Lid"))))
    duration_diff |>
        group_by(!!sym("judge")) |>
        summarise(n = n(),
                  rate = n / length(duration_diff$ID))
}


result_sim_compare <- result_sim_duration |>
    group_nest(seed_param, seed_sim, condition) |>
    mutate(duration_diff = map(data, ~
                               compare_Pro_Lid(.$duration_w[[1]]))) |>
    unnest(cols = duration_diff) |>
    dplyr::select(-data)

write_csv(result_sim_compare, file = "../results/05/sim_duration_compare.csv")


result_sim_compare_w <-  result_sim_compare |>
    pivot_wider(id_cols = c(seed_param, seed_sim, condition),
                names_from = judge,
                values_from = rate)

result_sim_duration_score <- result_sim_duration |>
    dplyr::select(seed_param, seed_sim, condition, duration) |>
    unnest(cols = duration) |>
    dplyr::filter(Drug != "Lid+Adr") |>
    group_by(seed_param, seed_sim, condition, Drug) |>
    mutate(score = scale(time)[, 1])



########################################
## save results
########################################

save(
    N, sigma_matrix,
    result_sim_duration,
    result_sim_duration_w,
    result_sim_duration_cor,
    result_sim_compare,
    result_sim_compare_w,
    result_sim_duration_score,
    file = "rda/sim_correlation_drugs.RData"
)

