#!/usr/bin/env Rscript

## SKG
## June 25, 2020
## Now with the outsiders
## Modeling the data results

## Model 1: smear
## Model 2: smear, HIV (yes/not yes)
## Model 3: smear, HIV, relative time

## Load libraries and data
devtools::load_all()
library(tidyverse)
data(tb_clean)
## library(optimParallel)
#
## Set up clusters
## cl <- makeCluster(5, type = "FORK")
##setDefaultCluster(cl)

my_seed <- 6172020
set.seed(my_seed)

t0 <- proc.time()[3]

## Format the data for all the variables

## I'm filling in the Unknown with a negative smear
## Should later check to see if results are changed with positive
tb_sub <- tb_clean %>%
    filter(PCR.Cluster != "") %>%
    mutate(smear = ifelse(spsmear == "Positive",
                          1, 0),
           cluster_id = PCR.Cluster,
           hiv_f = fct_collapse(hivstatus,
                              pos = c("Positive"),
                              neg = c("Negative"),
                              unk = c("Unknown", "Not offered",
                                      "Refused")),
           age = parse_number(ageatrept)) %>%
    mutate(hiv_neg_pos = ifelse(hiv_f == "neg", 1, 0),
           hiv_unk_pos = ifelse(hiv_f == "unk", 1, 0)) %>%
    group_by(cluster_id) %>%
    mutate(rel_time = as.numeric(datecoll - min(datecoll)) / 365) %>%
    mutate(med_rel_time = as.numeric(datecoll - median(datecoll)) / 365) %>%
    mutate(max_rel_time = as.numeric(datecoll - max(datecoll)) / 365) %>%
    mutate(age_imputed = ifelse(is.na(age), median(age, na.rm = TRUE),
                                age)) %>%
    mutate(min_time = as.numeric(min(datecoll)) / 1e4) %>%
    mutate(max_time = as.numeric(max(datecoll)) / 1e4) %>%
    mutate(time = as.numeric(datecoll) / 1e4) %>%
    mutate(person_id = paste0(cluster_id, "-", row_number())) %>%
    ungroup() %>% 
    mutate(race_f = fct_collapse(race,
                                 white = "White",
                                 black = "Black or African American",
                                 asian = "Asian")) %>%
    mutate(race_asian_white = ifelse(race_f == "asian", 1, 0),
           race_black_white = ifelse(race_f == "black", 1, 0)) %>%
    mutate(hiv_pos = ifelse(hiv_neg_pos == 0 & hiv_unk_pos == 0, 1,
                            0)) %>%
    mutate(smear_pos_hiv_pos = ifelse(hiv_pos == 1 & smear == 1,
                                      1, 0)) %>%
    mutate(age_scaled = scale(age_imputed)) %>%
    mutate(scaled_rel_time = scale(rel_time)) %>%
    mutate(scaled_med_rel_time = scale(med_rel_time)) %>%
    mutate(scaled_max_rel_time = scale(max_rel_time)) %>%
    select(cluster_id, smear,
           smear_pos_hiv_pos,
           hiv_pos,
           hiv_neg_pos,
           hiv_unk_pos,
           rel_time,
           med_rel_time,
           max_rel_time,
           race_asian_white,
           race_black_white,
           age_imputed,
           age_scaled,
           scaled_rel_time,
           scaled_max_rel_time,
           scaled_med_rel_time,
           max_time, min_time, time,
           person_id) 



## Separate parameter for outsiders
use_outsider_prob <- FALSE


## Make the lists for looping
n_mods <- 6
model_list <- vector(mode = "list", length = n_mods)
##
model_list[[1]] <- c("smear")
model_list[[2]] <- c("smear", "hiv_neg_pos", "hiv_unk_pos")
model_list[[3]] <- c("smear", "hiv_neg_pos", "hiv_unk_pos",
                     "rel_time")
model_list[[4]] <- c("smear", "hiv_neg_pos", "hiv_unk_pos",
                     "rel_time",
                     "race_asian_white", "race_black_white")
model_list[[5]] <- c("hiv_neg_pos", "hiv_unk_pos")
model_list[[6]] <- c("hiv_neg_pos", "hiv_unk_pos", "rel_time")
B <- 10000
n_mods <- length(model_list)
output_list <- vector(mode = "list", length = n_mods)
loglike_df <- data.frame(Model = 1:n_mods,
                         n_params = sapply(model_list, length) + 1,
                         loglike = 0,
                         aic = 0)


print(paste("Number of simulations to sample", B))
## Sampling all vars at once
all_vars <- c(unique(do.call('c', model_list)), "rel_time")
t1 <- proc.time()[3]
sampled_data1 <-  sample_gen_out_cond(tb_sub,
                                                B = B)

## joining features
feature_df <- tb_sub %>%
    select(person_id, one_of(all_vars))

sampled_data_all <- left_join(sampled_data1, feature_df,
              by = c("feature_id" = "person_id"))

t2 <- proc.time()[3] - t1
print(paste("Sampled data time:", round( t2 / 60, 3),
            "min"))







for(ii in 1:length(model_list)){
    print(paste("Model", ii))

    covariate_names <- model_list[[ii]]
    ## Subset data
    tb_df <- tb_sub %>%
        select(contains(c("cluster_id", covariate_names)))
    

    sampled_data <-  sampled_data_all %>%
        select(contains(c("cluster_id", "gen", "n_in_gen",
                          "inf_id", "id", "n_inf", "cluster_size",
                          "orig_id",
                          covariate_names)))

    
    ## Optimize
        
    print("Optimizing")
    init_params <- rep(0, length(covariate_names) + 1)
    bds <- rep(-5, length(covariate_names) + 1)
    if(length(covariate_names) > 5){
        bds <- rep(-4, length(covariate_names) + 1)
    }
    lower_bds <- bds
    upper_bds <- -bds
    if(use_outsider_prob){
        lower_bds <- c(.00001, lower_bds)
        upper_bds <- c(1, upper_bds)
        init_params <- c(.5, init_params)
    }
    cov_mat <- covariate_df_to_mat(sampled_data,
                                   cov_names = covariate_names)
    ## Now with data.table
    t1 <- proc.time()[3]
    best_params <- optim(par = init_params,
                     fn = general_loglike,
                     observed_data = tb_df,
                     sampled_data = data.table::as.data.table(sampled_data),
                     return_neg = TRUE,
                     cov_mat = cov_mat,
                     cov_names = covariate_names,
                     use_outsider_prob = use_outsider_prob,
                     method = "L-BFGS-B",
                     lower = lower_bds,
                     upper = upper_bds
                     )
    t2 <- proc.time()[3] - t1
    print(paste("Optimization time:", round( t2 / 60, 3),
                "min"))

    loglike_df$loglike[ii] <- - best_params$val

    ## Likelihood profile
    print("Likelihood profile")

    t1 <- proc.time()[3]

    ci_mat <- matrix(0, ncol = 2,
                     nrow = length(best_params$par))

    ci_mat <- matrix(0, ncol = 2,
                 nrow = length(best_params$par))

    for(jj in 1:nrow(ci_mat)){
        max_loglike <- -best_params$value
        best_pars <- best_params$par
        lower <- best_pars[jj] - 3
        upper <- best_pars[jj] + 3
        if(use_outsider_prob){
            lower <- .000001
            upper <- .999999
        }
        lower <- uniroot(f = like_wrapper, c(lower, best_pars[jj]),
                         max_loglike = max_loglike,
                         best_pars = best_pars,
                         beta_index = jj,
                         sampled_data = data.table::as.data.table(sampled_data),
                         covariate_names = covariate_names,
                         cov_mat = cov_mat,
                         chi_stat = 1.92)
        ##
        upper <- uniroot(f = like_wrapper, c(best_pars[jj], upper),
                         max_loglike = max_loglike,
                         best_pars = best_pars,
                         beta_index = jj,
                         sampled_data = data.table::as.data.table(sampled_data),
                         covariate_names = covariate_names,
                         cov_mat = cov_mat,
                         chi_stat = 1.92)

        ci_mat[jj, 1] <- lower$root
        ci_mat[jj, 2] <- upper$root
    }
    est_pars <- cbind(best_pars, ci_mat)
    nms <- c("Intercept",
             covariate_names)
    if(use_outsider_prob){
        nms <- c("Prob. O", nms)
    }
    rownames(est_pars) <- nms
    colnames(est_pars) <- c("Mean", "Lower95", "Uppper95")

    print(est_pars, digits = 2)
    t2 <- proc.time()[3] - t1

               
   print(paste("Likelihood profile time:", round( t2 / 60, 3),
               "min"))

   ll <- list(covariate_names = covariate_names,
               est_pars = est_pars)
    output_list[[ii]] <- ll

    t0f <- proc.time()[3] - t0
    print(paste("Total time:", round(t0f / 3600, 4), "hrs"))

}


loglike_df <- loglike_df %>%
    mutate(aic = 2 * (n_params) - loglike)
loglike_df


data_out <- list(loglike_df = loglike_df,
                model_list = model_list,
                output_list = output_list,
                sampled_data = sampled_data_all,
                seed = my_seed)


current_time <- Sys.time()
current_time <- gsub(" ", "_", current_time)
current_time <- gsub(":", ".", current_time)
fn_base <- paste0("data_output_outsider_", current_time, ".RDS")

saveRDS(data_out, fn_base)
