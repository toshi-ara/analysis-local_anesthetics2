#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(GGally)
library(patchwork)


dir.create("Fig", recursive = TRUE, showWarnings = FALSE)
source("common_functions.R")
load("rda/sim_correlation_Mean_logSigma.RData")

label_r <- expression(paste(italic("r")[{mu - log ~ sigma}]))


########################################
## parameter (Mean / log Sigma)
########################################

theme_predProb <- list(
    theme_bw(),
    theme(
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 16)
    )
)

plot_predProb <- function(x) {
    ggplot(x, aes(!!sym("time"), !!sym("prob"))) +
        geom_line(aes(group = !!sym("ID")), alpha = 0.1) +
        labs(x = "Time (min)",
             y = "Predicted Probability",
             title = label_r) +
        scale_x_continuous(breaks = seq(0, 120, by = 40)) +
        scale_y_continuous(breaks = scales::pretty_breaks()) +
        facet_wrap(. ~ r)
}

p_predProb <- result_sim_prob |>
    group_nest(seed_param) |>
    mutate(p = map(data, ~ plot_predProb(.) + theme_predProb)) |>
    dplyr::select(-data)



########################################
## Duration: summary (mean, SD)
########################################

plot_sim_duration_summary <- function(.x) {
    p1 <- ggplot(.x, aes(!!sym("r"), !!sym("mean"))) +
        geom_line() +
        geom_point() +
        labs(x = label_r, y = "Mean of duration (min)") +
        coord_cartesian(ylim = c(50, 80))
    p2 <- ggplot(.x, aes(!!sym("r"), !!sym("SD"))) +
        geom_line() +
        geom_point() +
        labs(x = label_r, y = "SD of duration (min)") +
        coord_cartesian(ylim = c(0, 20))

    return(p1 + p2)
}

plot_sim_duration_summary_all <- function(.x) {
    p1 <- ggplot(.x, aes(!!sym("r"), !!sym("mean"))) +
        geom_line(aes(group = factor(!!sym("seed_sim"))), alpha = 0.5) +
        labs(x = label_r, y = "Mean of duration (min)") +
        coord_cartesian(ylim = c(50, 80))
    p2 <- ggplot(.x, aes(!!sym("r"), !!sym("SD"))) +
        geom_line(aes(group = factor(!!sym("seed_sim"))), alpha = 0.5) +
        labs(x = label_r, y = "SD of duration (min)") +
        coord_cartesian(ylim = c(0, 20))

    return(p1 + p2)
}

theme_duration_summary <- list(
    theme_bw(),
    theme(
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16)
    )
)

p_sim_duration_summary <- result_sim_duration_summary |>
    dplyr::filter(seed_sim == 1) |>
    group_nest(seed_param) |>
    mutate(p = map(data, ~
                   plot_sim_duration_summary(.) &
                       theme_duration_summary)) |>
    dplyr::select(-data)

p_sim_duration_summary_all <- result_sim_duration_summary |>
    group_nest(seed_param) |>
    mutate(p = map(data, ~
                   plot_sim_duration_summary_all(.) &
                       theme_duration_summary)) |>
    dplyr::select(-data)



########################################
## SFig
## Duration (violin plot + jitter)
########################################

plot_sim_duration <- function(x) {
    ggplot(x, aes(factor(!!sym("r")), !!sym("time"))) +
        geom_violin() +
        geom_jitter(size = 0.5, width = 0.1, alpha = 0.2) +
        scale_x_discrete(breaks = seq(-1, 1, by = 0.5)) +
        labs(x = label_r, y = "Duration (min)")
}

theme_duration_violin <- list(
    theme_bw(),
    theme(
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white")
    )
)

p_sim_duration <- result_sim_duration |>
    group_nest(seed_param) |>
    mutate(p = map(data, ~
                   plot_sim_duration(.) +
                       facet_wrap(~ seed_sim) +
                       theme_duration_violin)) |>
    dplyr::select(-data)



########################################
## multipanel plot
########################################

p_sim_combine <- left_join(p_predProb, p_sim_duration_summary,
                           by = "seed_param") |>
    mutate(p = map2(p.x, p.y,
                    function(x, y) {
                     (x / (y & ggtitle(""))) +
                         plot_layout(ncol = 1, heights = c(3, 1)) +
                         plot_annotation(tag_levels = "A") &
                         theme(plot.tag = element_text(size = 28))
                    })
    ) |>
    dplyr::select(seed_param, p)


## SFig
p_sim_combine_suppl <- p_predProb |>
    left_join(p_sim_duration, by = "seed_param") |>
    left_join(p_sim_duration_summary_all, by = "seed_param") |>
    mutate(p = pmap(list(p1 = p.x, p2 = p.y, p3 = p),
                    function(p1, p2, p3) {
                     (p1 / p2 / (p3 & ggtitle(""))) +
                         plot_layout(ncol = 1, heights = c(3, 3, 1)) +
                         plot_annotation(tag_levels = "A") &
                         theme(plot.tag = element_text(size = 28))
                    })
    ) |>
    dplyr::select(seed_param, p)



cairo_pdf(filename = "../results/Fig/Fig3.pdf",
          width = 8, height = 10)
print(p_sim_combine$p[[1]])
dev.off()



########################################
## Supplemental Figure
########################################

cairo_pdf(filename = "../results/Fig/SFig1.pdf",
          width = 9, height = 16, onefile = TRUE)
print(p_sim_combine_suppl$p)
dev.off()

