#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)
library(patchwork)


########################################
## parameter (Mean / log Sigma)
########################################

dir.create("Fig", recursive = TRUE, showWarnings = FALSE)
source("common_functions.R", local = TRUE)
load("rda/animal_predProb.RData")
load("rda/animal_correlation_mean_logSigma.RData")



#########################################
### Correlation bwtween Mean and logSigma
#########################################

theme_scatter_mean_logSigma <- list(
    scale_x_continuous(breaks = scales::pretty_breaks()),
    scale_y_continuous(breaks = scales::pretty_breaks()),
    theme_bw(),
    theme(
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 16),
        strip.background = element_blank()
    )
)

p_param_pair_mean_logSigma <- lapply(dat_param_drug, function(.x) {
    .x |>
        dplyr::filter(!score) |>
        pivot_wider(names_from = "parameter", values_from = "value") |>
        ggplot(aes(Mean, logSigma)) +
            geom_point() +
            labs(x = label_mean, y = label_logSigma) +
            facet_wrap(~ drug, scales = "free", nrow = 1) +
            theme_scatter_mean_logSigma
})


####################
## probability curve of score
####################

theme_predProb <- list(
    scale_x_continuous(breaks = scales::pretty_breaks()),
    scale_y_continuous(breaks = scales::pretty_breaks()),
    theme_bw(),
    theme(
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        strip.background = element_blank()
    )
)

p_pred_animal1 <- ggplot(result_predProb_animal1,
                         aes(time, prob)) +
    geom_line(aes(group = ID), alpha = 0.3) +
    labs(x = "Time (min)", y = "Predicted probability") +
    facet_wrap(~ drug, nrow = 1) +
    theme_predProb


########################################
## multipanel plot
########################################

p <- (p_param_pair_mean_logSigma[[1]] / p_pred_animal1) +
    plot_annotation(tag_levels = "A") &
    theme(
        plot.tag = element_text(size = 32)
    )


cairo_pdf(file = "../results/Fig/Fig1.pdf",
          width = 9, height = 6)
print(p)
dev.off()

