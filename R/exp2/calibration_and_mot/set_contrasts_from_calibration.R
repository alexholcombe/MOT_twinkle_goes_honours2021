#!C:/Program Files/R/R-3.6.0/bin Rscript

args <- commandArgs(trailingOnly = T)
subj_id  <- as.integer(args[1])

source("D:/Documents/git/motNoise/R/utils.R")

data_dir <- "D:/Documents/git/motNoise/data/exp2/detection/range_estim/results"
data_file <- list.files(data_dir, pattern = sprintf("%d_.*.csv", subj_id), full.names = T)
df <- read.csv(data_file)
df$correct <- df$b_key_resp_2_.corr

params_fit <- mle_SDT_based_on_range_estim_results(df)

nll_fit <- nll(c(params_fit$cT,params_fit$beta),contrasts = df$contr_level,correct = df$correct)

contrasts <- c(inv_p_corr_fun(0.85,params_fit$cT, params_fit$beta),
inv_p_corr_fun(0.90,params_fit$cT, params_fit$beta),
inv_p_corr_fun(0.95,params_fit$cT, params_fit$beta))

               
contrasts
