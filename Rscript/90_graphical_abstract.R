#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(GGally)
library(ggiraph)

source("common_functions.R")


dir.create("Fig", recursive = TRUE, showWarnings = FALSE)
load("rda/sim_correlation_drugs.RData")



########################################
# Set parameters of each individual
########################################

Time <- seq(0, 180, by = 5)
N <- 5000


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


condition <- tibble(
    Condition = c(1, 4),
    r = list(
        sigma_matrix1,
        sigma_matrix4
    )
)



########################################
# Start simulation
########################################

seed_param <- 1
seed_sim <- 1


result_sim <- condition |>
    mutate(## set parameters
           param_drug_sim = pmap(
               list(r = r),
               set_parameters_r,
               N = N,
               parameters = param_drug, d_sigma = d_sigma,
               drug_name = drug_name,
               seed = seed_param),
           ## simulation
           res = pmap(
               list(x = param_drug_sim),
               doSimulation,
               time = Time, n_stim = 6, seed = seed_sim),
           duration = map(res, ~ groupSurvData(., k = 6)),
           ## convert wider
           parameter_mean_w = map(param_drug_sim, ~
                 dplyr::select(., Drug, ID, Mean) |>
                 pivot_wider(id_cols = ID,
                             names_from = "Drug", values_from = Mean) |>
                 dplyr::select(-ID, -`Lid+Adr`)),
           duration_w = map(duration, ~
                pivot_wider(., id_cols = ID,
                            names_from = Drug, values_from = time) |>
                dplyr::select(-ID, -`Lid+Adr`))
    )



########################################
## Scatter plot of duration among drugs
########################################

num_pair <- combn(1:4, 2, simplify = FALSE)
drug_pair <- combn(c("Pro", "Lid", "Mep", "Bup"), 2, simplify = FALSE)

.theme <- list(
    theme_void(),
    theme(
        # axis.text = element_text(size = 11, color = "black"),
        axis.text = element_blank(),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        legend.position = "none"
    )
)

plot_2d_density <- function(.x, drug_pair, .theme) {
    lapply(drug_pair, function(.pair) {
        ggplot(.x, aes(!!sym(.pair[1]), !!sym(.pair[2]))) +
            # geom_density_2d() +
            geom_density_2d_filled() +
            .theme
    })
}


p_parameter <- result_sim |>
    transmute(p = map(parameter_mean_w, ~
                      ggpairs(., aes(alpha = 0.5),
                              upper = list( continuous = NULL),
                              lower = list( continuous = NULL),
                              progress = FALSE) +
                              .theme
             ),
             p_density = pmap(
                 list(.x = parameter_mean_w),
                 plot_2d_density,
                 drug_pair = drug_pair,
                 .theme = .theme
             )
    )

p_duration <- result_sim |>
    transmute(p = map(duration_w, ~
                      ggpairs(., aes(alpha = 0.5),
                              upper = list( continuous = NULL),
                              lower = list( continuous = NULL),
                              # xlab = "Duration (min)",
                              progress = FALSE) +
                              .theme
             ),
             p_density = pmap(
                 list(.x = duration_w),
                 plot_2d_density,
                 drug_pair = drug_pair,
                 .theme = .theme
             )
    )


## insert to matrix
for (i in 1:2) {
    for (j in seq_len(length(num_pair))) {
        p_parameter$p[[i]][num_pair[[j]][2], num_pair[[j]][1]] <-
            p_parameter$p_density[[i]][[j]]
        p_duration$p[[i]][num_pair[[j]][2], num_pair[[j]][1]] <-
            p_duration$p_density[[i]][[j]]
    }
}


for (i in 1:2) {
    dsvg(file = sprintf("Fig/parameter%d.svg", i),
         width = 6, height = 6,
         standalone = TRUE,
         fonts = list(serif = "TeX Gyre Termes",
                      sans = "TeX Gyre Heros",
                      mono = "TeX Gyre Cursor",
                      symbol = "TeX Gyre Termes Math"))
    print(p_parameter$p[[i]])
    dev.off()

    dsvg(file = sprintf("Fig/duration%d.svg", i),
         width = 6, height = 6,
         standalone = TRUE,
         fonts = list(serif = "TeX Gyre Termes",
                      sans = "TeX Gyre Heros",
                      mono = "TeX Gyre Cursor",
                      symbol = "TeX Gyre Termes Math"))
    print(p_duration$p[[i]])
    dev.off()
}

