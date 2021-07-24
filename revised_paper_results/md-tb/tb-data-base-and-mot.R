#!/usr/bin/env Rscript

rep_from_lib <- TRUE

if(rep_from_lib){
    library(InfectionTrees)
} else {
    devtools::load_all()
}

library(tidyr)
library(dplyr)
library(kableExtra)
library(ggplot2)
library(forcats)
library(readr)
theme_set(theme_bw() + theme(axis.title = element_text()))




## --------------------------------------------------------------------------------------------------------
clusters <- tb_clean %>%
    dplyr::mutate(smear = ifelse(spsmear == "Positive",
                                 1, 0),
                  cluster_id = group,
                  hiv_f = ifelse(hivstatus == "Negative", "neg",
                          ifelse(hivstatus == "Positive", "pos",
                                 "unk"))) %>%
    dplyr::mutate(hiv_neg_pos = ifelse(hiv_f == "neg", 1, 0),
                  hiv_unk_pos = ifelse(hiv_f == "unk", 1, 0)) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::mutate(rel_time = as.numeric(rel_time / 365)) %>%
    dplyr::mutate(cluster_size = dplyr::n()) %>%
    dplyr::ungroup() %>%
    mutate(race_f = fct_collapse(race,
                                 white = "White",
                                 black = "Black or African American",
                                 asian = "Asian")) %>%
    mutate(race_asian_white = ifelse(race_f == "asian", 1, 0),
           race_black_white = ifelse(race_f == "black", 1, 0)) %>%
    select(cluster_id, smear,
           hiv_neg_pos,
           hiv_unk_pos,
           rel_time,
           race_asian_white,
           race_black_white,
           cluster_size)


## ----results = 'hide'------------------------------------------------------------------------------------
K <- 10000
my_seed <- 6172020
set.seed(my_seed)

## MODELS
## models
covariate_list <- vector(mode = "list", length = 5)
covariate_list[[1]] <- NA
covariate_list[[2]] <- "smear"
covariate_list[[3]] <- c("smear",
                         "hiv_neg_pos", "hiv_unk_pos")
covariate_list[[4]] <- c("smear",
                         "hiv_neg_pos", "hiv_unk_pos",
                         "rel_time")
covariate_list[[5]] <- c("smear",
                         "hiv_neg_pos", "hiv_unk_pos",
                         "rel_time",
                         "race_asian_white",
                         "race_black_white")

## Set up outputs
loglike_df <- data.frame(model = 1:length(covariate_list),
                         n_params = c(1, 2, 4, 5, 7),
                         loglike = 0,
                         aic = 0)
beta_mat1 <- matrix(0, nrow = 1, ncol = 4)
rownames(beta_mat1) <- c("Intercept")
colnames(beta_mat1) <- c("Est.", "lower", "upper", "SE")
beta_list <- vector(mode = "list", length = length(covariate_list))
beta_list[[1]] <- beta_mat1
for(ii in 2:length(covariate_list)){
    mat <- matrix(0, nrow = length(covariate_list[[ii]]) + 1,
                  ncol = 4)
    rownames(mat) <- c("Intercept", covariate_list[[ii]])
    colnames(mat) <- c("Est.", "lower", "upper", "SE")
    beta_list[[ii]] <- mat
}



## ----results = 'hide'------------------------------------------------------------------------------------
t_init <- proc.time()[3]

## Sample MC trees all at once

t0 <- proc.time()[3]
## Sample all MC clusters at once for each of the 5 models
mc_trees <-  sample_mc_trees(clusters,
                             B = K,
                             multiple_outside_transmissions = FALSE,
                             covariate_names = covariate_list[[length(covariate_list)]])
print(proc.time()[3] - t0)


## Fit each of the models

for(jj in 1:length(covariate_list)){
    covariate_names <- covariate_list[[jj]]
    print("Model:")
    print(covariate_names)
    if(is.na(covariate_names[1])){
        init_params <- 0
    } else{
        init_params <- rep(0, length(covariate_names) + 1)
    }
    ## Optimize

    print("Optimizing")
    bds <- rep(-5, length(init_params))
    if(length(covariate_names) > 5){
        bds <- rep(-4, length(init_params))
    }
    lower_bds <- bds
    upper_bds <- -bds
    cov_mat <- covariate_df_to_mat(mc_trees,
                                   cov_names = covariate_names)

    t1 <- proc.time()[3]
    best_params <- optim(par = init_params,
                         fn = general_loglike,
                         mc_trees = data.table::as.data.table(mc_trees),
                         return_neg = TRUE,
                         cov_mat = cov_mat,
                         cov_names = covariate_names,
                         use_outsider_prob = FALSE,
                         multiple_outside_transmissions = FALSE,
                         method = "L-BFGS-B",
                         lower = lower_bds,
                         upper = upper_bds,
                         hessian = TRUE
                         )
    t2 <- proc.time()[3] - t1
    print(paste("Optimization time:", round( t2 / 60, 3),
                "min"))


    beta_list[[jj]][,1] <- best_params$par
    beta_list[[jj]][, 4] <- sqrt(diag(solve(best_params$hessian))) ## SE from Fisher info



    print("best params:")
    print(beta_list[[jj]])
    loglike_df$loglike[jj] <- -best_params$val

    print(paste("Total time:", round( (proc.time()[3] - t_init)  / 3600, 3),
                "hrs"))


}



## --------------------------------------------------------------------------------------------------------

loglike_df <- loglike_df %>%
    mutate(aic = -loglike + 2 * n_params,
           model = 1:5) %>%
    select(model, everything())


loglike_df %>% kable(digits = 2,
                     col.names = c("Model", "# params.", "Log like.", "AIC")) %>%
    kable_styling(bootstrap_options = c("condensed", "hover", "striped", "responsive"),
                  full_width = FALSE, position = "center")



## --------------------------------------------------------------------------------------------------------
beta_list[[4]] %>%
    kable(digits = 2) %>%
    kable_styling(bootstrap_options = c("condensed", "hover", "striped", "responsive"),
                  full_width = FALSE, position = "center")



saveRDS(list(beta_list = beta_list,
             loglike_df = loglike_df),
        "base-results.RDS")


## ----results = 'hide'------------------------------------------------------------------------------------
K <- 10000
my_seed <- 24
set.seed(my_seed)

## MODELS
## models
covariate_list <- vector(mode = "list", length = 5)
covariate_list[[1]] <- NA
covariate_list[[2]] <- "smear"
covariate_list[[3]] <- c("smear",
                         "hiv_neg_pos", "hiv_unk_pos")
covariate_list[[4]] <- c("smear",
                         "hiv_neg_pos", "hiv_unk_pos",
                         "rel_time")
covariate_list[[5]] <- c("smear",
                         "hiv_neg_pos", "hiv_unk_pos",
                         "rel_time",
                         "race_asian_white",
                         "race_black_white")

## Set up outputs
loglike_df <- data.frame(model = 1:length(covariate_list),
                         n_params = c(1, 2, 4, 5, 7),
                         loglike = 0,
                         aic = 0)
beta_mat1 <- matrix(0, nrow = 1, ncol = 4)
rownames(beta_mat1) <- c("Intercept")
colnames(beta_mat1) <- c("Est.", "lower", "upper", "SE")
beta_list <- vector(mode = "list", length = length(covariate_list))
beta_list[[1]] <- beta_mat1
for(ii in 2:length(covariate_list)){
    mat <- matrix(0, nrow = length(covariate_list[[ii]]) + 1,
                  ncol = 4)
    rownames(mat) <- c("Intercept", covariate_list[[ii]])
    colnames(mat) <- c("Est.", "lower", "upper", "SE")
    beta_list[[ii]] <- mat
}



## ----results = 'hide'------------------------------------------------------------------------------------
t_init <- proc.time()[3]

## Sample MC trees all at once

t0 <- proc.time()[3]
## Sample all MC clusters at once for each of the 5 models
mc_trees <-  sample_mc_trees(clusters,
                             B = K,
                             multiple_outside_transmissions = TRUE,
                             covariate_names = covariate_list[[length(covariate_list)]])
print(proc.time()[3] - t0)


## Fit each of the models

for(jj in 1:length(covariate_list)){
    covariate_names <- covariate_list[[jj]]
    print("Model:")
    print(covariate_names)
    if(is.na(covariate_names[1])){
        init_params <- 0
    } else{
        init_params <- rep(0, length(covariate_names) + 1)
    }
    ## Optimize

    print("Optimizing")
    bds <- rep(-5, length(init_params))
    if(length(covariate_names) > 5){
        bds <- rep(-4, length(init_params))
    }
    lower_bds <- bds
    upper_bds <- -bds
    cov_mat <- covariate_df_to_mat(mc_trees,
                                   cov_names = covariate_names)

    t1 <- proc.time()[3]
    best_params <- optim(par = init_params,
                         fn = general_loglike,
                         mc_trees = data.table::as.data.table(mc_trees),
                         return_neg = TRUE,
                         cov_mat = cov_mat,
                         cov_names = covariate_names,
                         use_outsider_prob = FALSE,
                         multiple_outside_transmissions = TRUE,
                         method = "L-BFGS-B",
                         lower = lower_bds,
                         upper = upper_bds,
                         hessian = TRUE
                         )
    t2 <- proc.time()[3] - t1
    print(paste("Optimization time:", round( t2 / 60, 3),
                "min"))


    beta_list[[jj]][,1] <- best_params$par
    beta_list[[jj]][, 4] <- sqrt(diag(solve(best_params$hessian))) ## SE from Fisher info



    print("best params:")
    print(beta_list[[jj]])
    loglike_df$loglike[jj] <- -best_params$val

    print(paste("Total time:", round( (proc.time()[3] - t_init)  / 3600, 3),
                "hrs"))


}



## --------------------------------------------------------------------------------------------------------

loglike_df <- loglike_df %>%
    mutate(aic = -loglike + 2 * n_params,
           model = 1:5) %>%
    select(model, everything())


loglike_df %>% kable(digits = 2,
                     col.names = c("Model", "# params.", "Log like.", "AIC")) %>%
    kable_styling(bootstrap_options = c("condensed", "hover", "striped", "responsive"),
                  full_width = FALSE, position = "center")



## --------------------------------------------------------------------------------------------------------
beta_list[[4]] %>%
    kable(digits = 2) %>%
    kable_styling(bootstrap_options = c("condensed", "hover", "striped", "responsive"),
                  full_width = FALSE, position = "center")


saveRDS(list(beta_list = beta_list,
             loglike_df = loglike_df),
        "mot-results.RDS")
