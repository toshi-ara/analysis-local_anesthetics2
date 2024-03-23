#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(GGally)
library(ggpubr)
library(patchwork)

dir.create("Fig", recursive = TRUE, showWarnings = FALSE)
source("common_functions.R")
load("rda/sim_correlation_Pro_Lid.RData")



########################################
## Scatter plot of duration between Pro and Lid
########################################

label_mean2 <- expression(paste(
                   {italic("r")[mu]}, " between Pro and Lid"
               ))
label_logSigma2 <-expression(paste(
                   {italic("r")[log] [sigma]}, " between Pro and Lid"
               ))


## function
plot_scatter_Pro_Lid <- function(x, xlabel, ylabel, .theme = NULL) {
    p <- ggplot(x, aes(!!sym("Pro"), !!sym("Lid"))) +
        geom_abline(intercept = 0, slope = 1,
                    color = "gray", linewidth = 1) +
        geom_point(size = 1, alpha = 0.3) +
        labs(x = "Duration of Pro (min)", y = "Duration of Lid (min)") +
        scale_x_continuous(breaks = seq(0, 120, by = 40)) +
        scale_y_continuous(breaks = seq(0, 120, by = 40)) +
        facet_grid(r_logSigma ~ r_mean) +
        .theme

    p <- annotate_figure(p,
        top = text_grob(xlabel, size = 21),
        right = text_grob(ylabel, size = 21, rot = -90),
        fig.lab.face = "normal"
    )
    return(p)
}

theme_scatter_Pro_Lid <- list(
    theme_bw(),
    theme(
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        aspect.ratio = 1
    )
)

p_scatter_Pro_Lid <- result_sim_duration_Pro_Lid |>
    group_nest(seed_param) |>
    mutate(p = map(data, ~
                   plot_scatter_Pro_Lid(
                       .,
                       xlabel = label_mean2,
                       ylabel = label_logSigma2,
                       theme_scatter_Pro_Lid
                   ))) |>
    dplyr::select(-data)



########################################
## correlation of coefficient
########################################

theme_cor <- list(
    theme_bw(),
    theme(
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom"
     )
)

p_cor_Pro_Lid <- result_sim_Pro_Lid |>
    group_nest(seed_param) |>
    mutate(p = map(data, ~
                   ggplot(. |> mutate(r_logSigma = factor(r_logSigma)),
                          aes(r_mean, r)) +
                       geom_line(aes(group = r_logSigma, linetype = r_logSigma)) +
                       geom_point(aes(shape = r_logSigma), size = 3) +
                       labs(x = label_mean2,
                            y = "Corralation coefficient\nof duration",
                            shape = label_logSigma2,
                            linetype = label_logSigma2
                       ) +
                       scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
                       coord_cartesian(ylim = c(0, NA)) +
                       theme_cor
                  )
    )



########################################
## multipanel plot
########################################

layout <- "
AAAAAA
AAAAAA
AAAAAA
#BBBB#
"

p_sim_combine <- left_join(p_scatter_Pro_Lid, p_cor_Pro_Lid,
                           by = "seed_param") |>
    mutate(p = map2(p.x, p.y,
                    function(x, y) {
                     (x / (y & ggtitle(""))) +
                         plot_layout(design = layout) +
                         plot_annotation(tag_levels = "A") &
                         theme(plot.tag = element_text(size = 32,
                                                       face = "plain"))
                    })
    ) |>
    dplyr::select(seed_param, p)

cairo_pdf(file = "../results/Fig/Fig4.pdf",
          width = 8, height = 12.5, onefile = TRUE)
print(p_sim_combine$p[[1]])
dev.off()

cairo_pdf(file = "../results/Fig/SFig2.pdf",
          width = 8, height = 12.5, onefile = TRUE)
print(p_sim_combine$p)
dev.off()

