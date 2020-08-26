#!/usr/bin/env Rscript

## SKG
## Aug. 26, 2020
## Bootstrap simulations to estimate standard error for the data
## 5 different models

covariate_list <- vector(mode = "list",
                         length = 5)




## Load libraries and data
remove.packages("InfectionTrees")
if(!require(InfectionTrees)){
    devtools::install_github("skgallagher/InfectionTrees")
}
library(InfectionTrees)
library(tidyverse)
data(tb_clean)


## simulation parameters
K <- 10000 # of MC trees to draw for each cluster
n_models <- 5
n_boot <- 100 # number of bootstrap simulations
my_seed <- 8262020
## Each inner loop takes ~1 hr

## Format original clusters to useful variables

orig_clusters <-  tb_clean %>%
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
    mutate(cluster_size = n()) %>%
    ungroup() %>% 
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
## Set up beta matrices
beta_mat1 <- matrix(0, ncol = 1,
                    nrow = n_boot)
colnames(beta_mat1) <- c("Intercept")
beta_list <- vector(mode = "list", length = 5)
beta_list[[1]] <- beta_mat1
for(ii in 2:5){
    mat <- matrix(0, ncol = length(covariate_list[[ii]]) + 1,
                  nrow = n_boot)
    colnames(mat) <- c("Intercept", covariate_list[[ii]])
    beta_list[[ii]] <- mat
}






    


t_init <- proc.time()[3]
for(nn in 1:n_boot){
    set.seed(my_seed + nn)
    print(sprintf("Bootstrap simulation %d", nn))

    ## Sample a bootstrap cluster
    b_clusters <- bootstrap_clusters(clusters = orig_clusters)

    print("sampling MC trees")
    ## sample MC trees

    t0 <- proc.time()[3]
    ## Sample all MC clusters at once for each of the 5 models
    mc_trees <-  sample_mc_trees(b_clusters,
                                 B = K,
                                 multiple_outside_transmissions = FALSE,
                                 covariate_names = covariate_list[[length(covariate_list)]])
    print(proc.time()[3] - t0)

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
    ## Now with data.table
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
                         upper = upper_bds
                         )
    t2 <- proc.time()[3] - t1
    print(paste("Optimization time:", round( t2 / 60, 3),
                "min"))
    print("best params:")
    print(best_params$par)
    beta_list[[jj]][nn,] <- best_params$par
    print(paste("Total time:", round( (proc.time()[3] - t_init)  / 3600, 3),
                "hrs"))

    }
}

se_list <- lapply(beta_list, function(x){
    apply(x, 2, sd)})
print("Bootstrap sample error matrices")
print(se_list)
    
output_list <- list(covariate_list = covariate_list,
                    se_list = se_list,
                    beta_list = beta_list)

current_time <- Sys.time()
current_time <- gsub(" ", "_", current_time)
current_time <- gsub(":", ".", current_time)
fn_base <- paste0("base_bootstrap_sims_",
                  current_time, ".RDS")

saveRDS(output_list, fn_base)


