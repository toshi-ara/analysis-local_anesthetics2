#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(RColorBrewer)
library(xtable)
source("common_functions.R", local = TRUE)


########################################
## survival analysys in various r between mean and logSigma
########################################

dir.create("../results/Fig", recursive = TRUE, showWarnings = FALSE)
dir.create("../results/Table", recursive = TRUE, showWarnings = FALSE)
load("rda/animal_sim_survival.RData")



#################
## save as LaTeX
#################

## function
strMedianCI <- function(x1, x2, x3, str = "\\\\text{--}") {
    sprintf("%.1f [%.1f, %.1f]", x1, x2, x3) |>
        gsub(pattern = "NA", replacement = str)
}


tbl1 <- result_survival_tbl |>
    dplyr::filter(seed_param == 1 & seed_sim == 1) |>
    mutate(Drug = "",
           median_ci = strMedianCI(median, lwr, upr)) |>
    dplyr::select(Drug, Condition, n, events, median_ci)

tbl1$Drug[(0:4) * 5 + 1] <- drug_name

colnames(tbl1) <- c(
    "Drug",
    "Condition",
    "$n$",
    "Events",
    "Median [95\\% CI]"
)


print(file = "../results/Table/Table5.tex",
    xtable(tbl1, digits = c(0, 0, 0, 0, 0, 0)),
    math.style.negative = TRUE,
    NA.string = "-",
    include.rownames = FALSE,
    # include.colnames = FALSE,
    only.contents = TRUE,
    hline.after = (0:5) * 5,
    booktabs = TRUE,
    timestamp = NULL,
    sanitize.text.function = identity
)


#################
## plot
#################

theme_survival <- list(
    labs(x = "Time (min)",
         y = "Duration of drugs effects",
         color = "Condition"),
    geom_hline(yintercept = 0.5, linetype = "dotted"),
    coord_cartesian(xlim = c(0, 180)),
    scale_color_brewer(palette = "Dark2", labels = lab_condition),
    theme(
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)
    )
)


p_surv <- lapply(seq_len(nrow(result_survival)), function(i) {
    p <- ggsurvplot_facet(
        result_survival$res_surv[[i]],
        result_survival$survdata[[i]],
        size = 0.7, alpha = 0.8,
        surv.scale = "percent",
        censor.shape = "|",
        short.panel.labs = TRUE,
        panel.labs.font = list(size = 16),
        panel.labs.background = list(color = "white", fill = "white"),
        facet.by = "Drug", nrow = 2,
        ggtheme = theme_bw())

    p$data$Drug <- factor(p$data$Drug, levels = drug_name)
    p + theme_survival
})


## Save plots (PDF)
cairo_pdf(filename = "../results/Fig/SFig4.pdf",
          width = 10, height = 5, onefile = TRUE)
print(p_surv)
dev.off()

