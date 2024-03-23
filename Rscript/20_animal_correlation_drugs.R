#!/usr/bin/env Rscript

library(ggplot2)
library(GGally)
library(patchwork)


########################################
## combine (PCP / Mean / logSigma)
########################################

dir.create("Fig", recursive = TRUE, showWarnings = FALSE)
source("common_functions.R", local = TRUE)
load("rda/animal_correlation_mean_logSigma.RData")
load("rda/animal_duration_correlation.RData")



########################################
## Correlation of parameters (Mean and logSigma) between drugs
########################################

## Parallel Coordinates Plot
## correlation between drugs
theme_PCP <- list(
    theme_bw(),
    theme(
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 18),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 18),
        strip.background = element_blank()
    )
)


## PCP (score)
p_PCP <- lapply(dat_param_drug, function(.x) {
    .x |>
        dplyr::filter(score) |>
        mutate(parameter = factor(parameter,
                                  labels = c(label_mean, label_logSigma))) |>
        ggplot(aes(drug, value)) +
            geom_line(aes(group = ID), linewidth = 0.5, alpha = 0.3) +
            labs(y = "Standardized score") +
            facet_wrap(~ parameter,
                       labeller = label_parsed,
                       scales = "free_y") +
            theme_PCP
})


####################
## scatter plot (parameters)
####################
dat_mean <- lapply(dat_param_drug_w, function(.x) {
    .x |>
        dplyr::filter(!score & parameter == "Mean")
})

dat_logSigma <- lapply(dat_param_drug_w, function(.x) {
    .x |>
        dplyr::filter(!score & parameter == "logSigma")
})


theme_scatter_drugs <- list(
    scale_shape_manual(values = c(16, 4)),
    theme_bw(),
    theme(
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "none"
    )
)

####################
## GGally::ggpairs
####################
p_param_pair_mean <- lapply(dat_mean, function(.x) {
    ggpairs(.x, columns = 4:7,
            upper = list(
                continuous = wrap(ggally_cor,
                                  size = 5, color = "black",
                                  stars = FALSE)
            ),
            lower = NULL,
            xlab = label_mean, title = "",
            progress = FALSE) +
        theme_scatter_drugs
})

p_param_pair_logSigma <- lapply(dat_logSigma, function(.x) {
    ggpairs(.x, columns = 4:7,
            upper = list(
                continuous = wrap(ggally_cor,
                                  size = 5, color = "black",
                                  stars = FALSE)
            ),
            lower = NULL,
            xlab = label_logSigma, title = "",
            progress = FALSE) +
    theme_bw() +
    theme_scatter_drugs
})


####################
## plot for replacement (shape)
####################
num_pair <- combn(1:4, 2, simplify = FALSE)
drug_pair <- combn(c("Pro", "Lid", "Mep", "Bup"), 2, simplify = FALSE)

p_pair_mean <- lapply(dat_mean, function(.x) {
    lapply(drug_pair, function(.y) {
        .x |>
            ggplot(aes(!!sym(.y[1]), !!sym(.y[2]), shape = !not_outlier)) +
                geom_point() +
                scale_x_continuous(breaks = seq(0, 100, by = 25)) +
                scale_y_continuous(breaks = seq(0, 100, by = 25)) +
                theme_scatter_drugs
    })
})

p_pair_logSigma <- lapply(dat_logSigma, function(.x) {
    lapply(drug_pair, function(.y) {
        .x |>
            ggplot(aes(!!sym(.y[1]), !!sym(.y[2]), shape = !not_outlier)) +
                geom_point() +
                scale_x_continuous(breaks = 1:4) +
                scale_y_continuous(breaks = 1:4) +
                theme_scatter_drugs
    })
})

## insert to matrix
for (i in 1:2) {
    for (j in seq_len(length(num_pair))) {
        ## mean
        p_param_pair_mean[[i]][num_pair[[j]][2], num_pair[[j]][1]] <-
            p_pair_mean[[i]][[j]]
        ## logSigma
        p_param_pair_logSigma[[i]][num_pair[[j]][2], num_pair[[j]][1]] <-
            p_pair_logSigma[[i]][[j]]
    }
}



########################################
## Correlation of duration between drugs
########################################

theme_duration <- list(
    theme_bw(),
    theme(
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size = 16, hjust = 0.5)
    )
)

p_duration <- duration_w |>
    dplyr::select(-ID, -`Lid+Adr`) |>
    ggpairs(aes(alpha = 0.5), progress = FALSE) +
    labs(x = "Duration (min)") +
    geom_abline(intercept = 0, slope = 1,
                color = "gray", linewidth = 0.6) +
    theme_duration



########################################
## multipanel plot
########################################

layout <- "
AAAAA#
BBBCCC
BBBCCC
DDD###
DDD###
"

p_comb <- vector("list", 2)
for (i in seq_len(length(p_comb))) {
    p_comb[[i]] <- p_PCP[[i]] /
        wrap_elements(plot = ggmatrix_gtable(p_param_pair_mean[[i]])) /
        wrap_elements(plot = ggmatrix_gtable(p_param_pair_logSigma[[i]])) /
        wrap_elements(plot = ggmatrix_gtable(p_duration)) +
        plot_layout(design = layout) +
        plot_annotation(
            tag_levels = "A",
        ) &
        theme(plot.tag = element_text(size = 32))
}

cairo_pdf(file = "../results/Fig/Fig2.pdf",
          width = 10, height = 15)
print(p_comb[[1]])
dev.off()

cairo_pdf(file = "../results/Fig/Fig2_without_outlier.pdf",
          width = 10, height = 15)
print(p_comb[[2]])
dev.off()

