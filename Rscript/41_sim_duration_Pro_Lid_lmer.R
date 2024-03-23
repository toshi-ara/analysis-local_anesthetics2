#!/usr/bin/env Rscript

library(xtable)
library(dplyr)
library(lmerTest)
library(broom)
library(broom.mixed)

load("rda/sim_correlation_Pro_Lid.RData")


#####################################
## use broom::tidy / broom.mixed::tidy
#####################################

## model with interaction
lmer1 <- lmer(r ~ r_mean * r_logSigma + (1 | seed_param),
              data = result_sim_Pro_Lid |>
                         mutate(seed_param = factor(seed_param)))

## model without interaction
lmer2 <- lmer(r ~ r_mean + r_logSigma + (1 | seed_param),
              data = result_sim_Pro_Lid |>
                         mutate(seed_param = factor(seed_param)))

LRtest <- anova(lmer1, lmer2)



#####################################
## LaTeX
#####################################

print(file = "../results/Table/STable2.tex",
    xtable(tidy(lmer1),
           digits = c(0, 0, 0, 0, 3, 3, 3, 1, 3),
           caption = "Analysis by Linear Mixed-Effects Models with interaction",
           label = "Stbl_lmer1"),
    include.rownames = FALSE,
    caption.placement = "top",
    size = "small",
    booktabs = TRUE,
    timestamp = NULL
)


print(file = "../results/Table/tbl_sim_correlation_lmer1.tex",
    xtable(tidy(lmer1),
           digits = c(0, 0, 0, 0, 3, 3, 3, 1, 3),
           caption = "Analysis by Linear Mixed-Effects Models with interaction",
           label = "Stbl_lmer1"),
    include.rownames = FALSE,
    caption.placement = "top",
    size = "small",
    booktabs = TRUE,
    timestamp = NULL
)

print(file = "../results/Table/tbl_sim_correlation_lmer2.tex",
    xtable(tidy(lmer2),
           digits = c(0, 0, 0, 0, 3, 3, 3, 1, 3),
           caption = "Analysis by Linear Mixed-Effects Models without interaction",
           label = "Stbl_lmer2"),
    include.rownames = FALSE,
    caption.placement = "top",
    booktabs = TRUE,
    timestamp = NULL
)

print(file = "../results/Table/tbl_sim_correlation_LRtest.tex",
    xtable(tidy(LRtest),
           digits = c(0, 0, 0, 1, 1, 1, 1, 3, 0, 3),
           caption = "Likelihood ratio test of two models",
           label = "Stbl_lmer_LRtest"),
    include.rownames = FALSE,
    caption.placement = "top",
    booktabs = TRUE,
    timestamp = NULL
)

