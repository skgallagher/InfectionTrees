## September 13, 2020
## SKG
## Analyze in bootstrap results from Biowulf (on the real data)
## These are the SE from bootstrap for the different covariates

## BASE MODEL
base <- readRDS("bootstrap_sims_base_2020-08-31.RDS")

base_se <- base$se_vec

mot <- readRDS("bootstrap_sims_mot_2020-08-31.RDS")
mot_se <- mot$se_vec

## Naive
## NEEDS FILLED IN
naive_se <- rep(NA, length(mot_se))

se_vec <- c(base_se, mot_se, naive_se)
names(se_vec) <- NULL
se_vec

saveRDS(se_vec, "biowulf_bootstrap_results_8_30_20.RDS")
