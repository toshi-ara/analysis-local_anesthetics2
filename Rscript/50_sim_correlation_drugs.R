#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(GGally)
library(patchwork)
source("common_functions.R")

dir.create("Fig", recursive = TRUE, showWarnings = FALSE)
load("rda/sim_correlation_drugs.RData")


########################################
## Scatter plot of duration among drugs
########################################

num_pair <- combn(1:4, 2, simplify = FALSE)
drug_pair <- combn(c("Pro", "Lid", "Mep", "Bup"), 2, simplify = FALSE)

theme_scatter <- list(
    theme_bw(),
    theme(
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white")
    )
)

plot_scatter <- function(.x, drug_pair, .theme) {
    lapply(drug_pair, function(.pair) {
        ggplot(.x, aes(!!sym(.pair[1]), !!sym(.pair[2]))) +
            geom_abline(intercept = 0, slope = 1,
                        color = "gray", linewidth = 1) +
            geom_point(alpha = 0.5) +
            .theme
    })
}


p_sim_duration <- result_sim_duration_w |>
    group_nest(seed_param, seed_sim, condition) |>
    mutate(p = map(data, ~
                   ggpairs(., aes(alpha = 0.5), columns = 2:5,
                           upper = list(
                               continuous = wrap(ggally_cor,
                                                 size = 5, color = "black",
                                                 stars = FALSE)
                           ),
                           lower = list(
                               continuous = NULL
                           ),
                           xlab = "Duration (min)",
                           progress = FALSE) +
                           theme_scatter
           ),
           p_scatters = pmap(
               list(.x = data),
               plot_scatter,
               drug_pair = drug_pair,
               .theme = theme_scatter
           )
    )


## insert to matrix
for (i in seq_len(nrow(p_sim_duration))) {
    for (j in seq_len(length(num_pair))) {
        p_sim_duration$p[[i]][num_pair[[j]][2], num_pair[[j]][1]] <-
            p_sim_duration$p_scatters[[i]][[j]]
    }
}


####################
## combine
####################

## function
combine_ggpairs <- function(x1, x2, x3, x4) {
    lapply(list(x1, x2, x3, x4),
           function(i) {
               wrap_elements(
                   full = ggmatrix_gtable(i + ggtitle(""))
               )}) |>
        wrap_plots(ncol = 2) +
        plot_annotation(
            tag_prefix = "Condition", tag_levels = "1"
        ) &
        theme(
            plot.tag = element_text(size = 21, vjust = 2)
        )
}


p_sim_duration_combine <- p_sim_duration |>
    pivot_wider(id_cols =  c(seed_param, seed_sim),
                names_from = condition,
                names_prefix = "Condition",
                values_from = p) |>
    mutate(p = pmap(
               list(x1 = Condition1, x2 = Condition2,
                    x3 = Condition3, x4 = Condition4),
                combine_ggpairs)
    )  |>
    dplyr::select(seed_param, seed_sim, p)


cairo_pdf(filename = "../results/Fig/Fig5.pdf",
          width = 12, height = 12)
print(p_sim_duration_combine$p[[1]])
dev.off()

## SFig
cairo_pdf(filename = "../results/Fig/SFig3.pdf",
          width = 12, height = 12, onefile = TRUE)
print(p_sim_duration_combine$p)
dev.off()

