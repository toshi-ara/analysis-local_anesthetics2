library(dplyr)
library(purrr)
library(ggplot2)
library(mvtnorm)


###############################
## common values
###############################

## adjusted parameters
param_drug <- data.frame(
    Drug = 1:4,
    mu0 = c(68, 61, 50, 30),
    s_mu0 = c(10, 7, 7, 13),  # about half of parameters in model 1
    log_sigma0 = c(2.2, 2.4, 2.4, 2.5),
    log_s_sigma0 = c(0.4, 0.4, 0.4, 0.5),
    adr = 0.7
)
d_sigma <- 4

drug_name <- c("Pro", "Lid", "Mep", "Bup", "Lid+Adr")


## for plot
label_mean <- expression(mu)
label_logSigma <- expression(paste(log, " ", sigma))



###############################
## Get duration
###############################

## Get duration under the specified condition from logical vector
##
## input:
##   x: logical vector [TRUE/FALSE] in each time point
##   time: vector of time
##   SEQ: number of consecutive occurrences
## return:
##   time when the effect of local anesthesia finished
##   (otherwise Inf)
##
getDuration <- function(x, time, SEQ = 3) {
    n <- length(time)
    idx <- lapply(1:(n - SEQ + 1),
                  function(i) x[i:(i + SEQ - 1)]) |>
        sapply(all) |>
        match(x = TRUE) + SEQ - 1
    return(ifelse(is.na(idx), Inf, time[idx]))
}

## Get duration with status (for survival analysis)
##
## input:
##   x: data.frame(ID, Drug, time, score)
##   k: number of at least k successes
##   SEQ: number of consecutive occurrences
## return:
##   status: 0 (no-response/alive), 1 (response/dead)
##   return is tibble
##
groupSurvData <- function(x, k, SEQ = 3) {
    res <- x |>
        group_nest(!!sym("Drug"), !!sym("ID")) |>
        mutate(time = map_dbl(data,
                              ~ getDuration((.$score >= k), .$time, SEQ)),
               status = ifelse(time == Inf, 0, 1)) |>
        dplyr::select(-data) |>
        mutate(ID = as.integer(!!(sym("ID"))))
    return(res)
}



###############################
## Prediction from simulation results
###############################

## Return predict probability
##
## input:
##    time:
##    mu: mean of normal distribution
##    sigma: sigma of normal distribution
##    adr: parameter of adrenaline
## output:
##    probabilitity of response to one stimulus
##
predProb <- function(time, mu, sigma, adr) {
    conc <- 100 - (1 - adr) * time
    pnorm(conc, mu, sigma, lower.tail = FALSE)
}

## Return predict score (decimal)
##
## input:
##    time:
##    mu: mean of normal distribution
##    sigma: sigma of normal distribution
##    adr: parameter of adrenaline
## output:
##    value (decimal) of response to stimuli
##
predScore <- function(time, mu, sigma, adr, score = 6) {
    predProb(time, mu, sigma, adr) * score
}



########################################
# Simulation
########################################

## get individual parameters from drug parameters
##
## input:
##   r: matrix of parameters correlation (8 x 8)
##   N: number of individuals
##   parameters: parameter of drugs
##     (data.frame of Drug, mu0, s_mu0, log_sigma0, log_s_sigma0, adr)
##   d_sigma: sigma of offset
##   drug_name: vector of drug names
##   seed: value for set.seed()
##
## output:
##   data.frame(Drug, ID, Mean, Sigma, adr)
##
set_parameters_r <- function(r, N, parameters, d_sigma, drug_name, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    n <- nrow(r)

    ## multivariate normal distriburion for each individual
    rand <- rmvnorm(n = N, mean = rep(0, n), sigma = r)

    ## drug parameters
    param <- lapply(seq_len(nrow(parameters)), function(i) {
        data.frame(Drug = i,
                   ID = seq_len(N),
                   Mean = parameters$mu0[i] +
                          parameters$s_mu0[i] * rand[, i],
                   Sigma = exp(parameters$log_sigma0[i] +
                               parameters$log_s_sigma0[i] * rand[, i + n / 2]),
                   adr = 0)
    })

    ## Set parameter of Lid + Adr
    param_Lid <- param[[2]] |>
        mutate(Drug = 5, adr = parameters$adr[2])

    ## Combine
    param_drug_sim <- param |>
        bind_rows() |>
        bind_rows(param_Lid) |>
        mutate(Drug = factor(!!sym("Drug"), labels = drug_name))

    return(param_drug_sim)
}

## get individual parameters from two drug parameters for Pro and Lid
## input:
##   r_mean: correlation coefficient of mean between Pro and Lid
##   r_logSigma: correlation coefficient of logSigma between Pro and Lid 
##   N: number of individuals
##   parameters: parameter of drugs
##     (data.frame of Drug, mu0, s_mu0, log_sigma0, log_s_sigma0, adr)
##   d_sigma: sigma of offset
##   drug_name: vector of drug names
##   seed: value for set.seed()
##
## output:
##   data.frame(Drug, ID, Mean, Sigma, adr)
##
set_parameters_r_Pro_Lid <- function(r_mean, r_logSigma,
                                     N, parameters,
                                     d_sigma, drug_name,
                                     seed = NULL) {
    n <- 2
    ## correlation matrix
    Sigma_mean <- matrix(c(
        1,      r_mean,
        r_mean, 1     ), n, n)

    Sigma_logSigma <- matrix(c(
        1,          r_logSigma,
        r_logSigma, 1         ), n, n)

    if (!is.null(seed)) {
        set.seed(seed)
    }

    ## multivariate normal distriburion for each individual
    rand_mean <- rmvnorm(n = N, mean = rep(0, n), sigma = Sigma_mean)
    rand_logSigma <- rmvnorm(n = N, mean = rep(0, n), sigma = Sigma_logSigma)

    ## Offset value for each individual
    d <- rnorm(n = N, mean = 0, sd = d_sigma)

    ## drug parameters
    param <- lapply(seq_len(nrow(parameters)), function(i) {
        data.frame(Drug = i,
                   ID = seq_len(N),
                   Mean = parameters$mu0[i] +
                          parameters$s_mu0[i] * rand_mean[, i] + d,
                   Sigma = exp(parameters$log_sigma0[i] +
                               parameters$log_s_sigma0[i] * rand_logSigma[, i]),
                   adr = 0)
    })

    ## Combine
    param_drug_sim <- param |>
        bind_rows() |>
        mutate(Drug = factor(!!sym("Drug"), labels = drug_name))

    return(param_drug_sim)
}


## get individual parameters from one drug parameters
##
## input:
##   r: correlation bwtween Mean and logSigma
##   N: number of individuals
##   parameter: parameter of drugs (data.frame of Drug including mu0 log_sigma0)
##   d_sigma: sigma of offset
##   drug_name: [character]
##   seed: value for set.seed()
##
## output:
##   data.frame(Drug, ID, Mean, Sigma, adr)
##
set_parameters_one <- function(r, N, parameters, d_sigma, drug_name, seed = NULL) {
    n <- 2
    ## correlation matrix
    Sigma <- matrix(c(
        1, r,
        r, 1), n, n)


    if (!is.null(seed)) {
        set.seed(seed)
    }

    ## multivariate normal distriburion for each individual
    rand <- rmvnorm(n = N, mean = rep(0, n), sigma = Sigma)

    ## Offset value for each individual
    d <- rnorm(n = N, mean = 0, sd = d_sigma)

    ## drug parameters
    param <- data.frame(Drug = drug_name,
        ID = seq_len(N),
        Mean = parameters$mu0 +
               parameters$s_mu0 * rand[, 1] + d,
        Sigma = exp(parameters$log_sigma0 +
                    parameters$log_s_sigma0 * rand[, 2]),
        adr = 0
    )

    return(param)
}



## get simulation results (score values) from parameters of drugs
##
## input:
##   x: parameter of drugs (data.frame of Drug, ID, Mean, Sigma, adr)
##   time: vector of time point
##   n_stim: number of stimulation by needle
##   seed: value for set.seed()
##
## output:
##   data.frame(ID, Drug, time, score)
##
doSimulation <- function(x, time, n_stim, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    n <- nrow(x)
    res_list <- vector("list", n)

    for (i in seq_len(n)) {
        Param <- x[i, ]
        prob <- predProb(time, Param$Mean, Param$Sigma, Param$adr)
        res_list[[i]] <- data.frame(
            ID = Param$ID,
            Drug = Param$Drug,
            time = time,
            score = rbinom(n = length(time), size = n_stim, prob = prob)
        )
    }
    res <- bind_rows(res_list) |>
        mutate(ID = factor(!!sym("ID")),
               Drug = factor(!!sym("Drug")))
    return(res)
}


## plot simulation results each drug and individual
##
## input:
##   x: data.frame(time, score, breaks)
## output:
##   ggplot
##
plot_simlation <- function(x, breaks) {
    p <- ggplot(x, aes(!!sym("time"), !!sym("score"))) +
        geom_line(color = "gray30") +
        geom_point(size = 1) +
        labs(x = "Time (min)", y = "Score") +
        scale_x_continuous(breaks = breaks) +
        facet_grid(Drug ~ ID) +
        theme_bw() +
        theme(
            axis.text = element_text(size = 11, color = "black"),
            axis.title = element_text(size = 16),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 14),
            strip.text.y = element_text(size = 16)
        )
    return(p)
}


## get lower triangle of matrix as vector
##
## input:
##   x: matrix of correlation coefficients
##
## output:
##   vector
##
get_lower_tri_vector <- function(x) {
    diag(x) <- NA
    x[upper.tri(x)] <- NA
    x <- as.vector(x)
    x[!is.na(x)]
}

