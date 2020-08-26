#!/usr/bin/env Rscript

## SKG
## Aug. 20, 2020
## Now with the outsiders
## Simulate outsider outbreaks
## See if we can estimate the parameters



## Load libraries and data
## devtools::load_all()
if(!require(InfectionTrees)){
    devtools::install_github("skgallagher/InfectionTrees")
}
library(InfectionTrees)
library(tidyverse)

## simulation parameters
K <- 100 # number of clusters to simulate
n_sims <- 2 # total times to simulate per set of parameters
B <- 100 # number of MC draws per cluster

## EACH OUTER LOOP TAKES >30 HOURS.

params_df <- data.frame(beta0 = c(-2.5, -2.5,
                                  -2.75, -2.75,
                                  -2.5, -2.5,
                                  -2.5, -2.5,
                                  -1.5, -1.5,
                                  -1.5, -1.5,
                                  -1.5, -2.75,
                                  -1.5),
                        beta1 = c(-.5, -.5,
                                  0, -.25,
                                  0, 0,
                                  .5, .5,
                                  -.5, -1,
                                  0, 0,
                                  .5, 3,
                                  1),
                        gamma = c(.5, .7,
                                  .1, .1,
                                  .5, .7,
                                  .5, .7,
                                  .7, .5,
                                  .7, .5,
                                  .7, .1,
                                  .5))
params_df$set <- 1:nrow(params_df)

simulation_output <- vector(mode = "list", length = nrow(params_df))


for(par_index in 1:nrow(params_df)){
    my_seed <- 8152020 + par_index
    set.seed(my_seed)  # needs fixed?

    t0 <- proc.time()[3]
    gamma <- .7
    beta0 <- params_df$beta0[par_index]
    beta1 <- params_df$beta1[par_index]
    gamma <- params_df$gamma[par_index]
    print(sprintf("beta0 = %.2f, beta1 = %.2f, gamma = %.2f",
            beta0, beta1, gamma))
    inf_params <- c(beta0, beta1)
    covariates_df <- data.frame(x = c(rep(1, gamma * 10),
                                      rep(0, (1-gamma) * 10)
                                      ))
    covariate_names <- "x"
    ## Simulate an outbreak
    best_params_mat <- matrix(0, nrow = n_sims, ncol = length(inf_params))
    cluster_sizes_list <- vector(mode = "list", length = n_sims)
    t_init <- proc.time()[3]
    for(sim in 1:n_sims){

        print("simulation number")
        print(sim)


        print("simulating data set")
        df <- simulate_bp(K = K,
                          inf_params = inf_params,
                          sample_covariates_df = covariates_df,
                          covariate_names = covariate_names,
                          max_size = 70)


        outsider_obs <- df %>%
            filter(gen != 1) %>%
            mutate(cluster_size = cluster_size - 1)

        ## Number of 'observed' clusters
        length(unique(outsider_obs$cluster_id))

        cluster_size_df <- outsider_obs %>% dplyr::group_by(cluster_id) %>%
            summarize(cluster_size = n(),
                      .groups = "drop")
        cluster_sizes_list[[sim]] <- cluster_size_df$cluster_size
                                  

        print("sampling MC trees")
        ## sample MC trees

        t0 <- proc.time()[3]
        mc_trees <-  sample_mc_trees(outsider_obs,
                                     B = B,
                                     multiple_outside_transmissions = TRUE)
        print(proc.time()[3] - t0)

        
        ## Optimize
        
        print("Optimizing")
        init_params <- rep(0, length(covariate_names) + 1)
        bds <- rep(-5, length(covariate_names) + 1)
        if(length(covariate_names) > 5){
            bds <- rep(-4, length(covariate_names) + 1)
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
                             multiple_outside_transmissions = TRUE,
                             method = "L-BFGS-B",
                             lower = lower_bds,
                             upper = upper_bds,
                             hessian = TRUE
                             )
        t2 <- proc.time()[3] - t1
        print(paste("Optimization time:", round( t2 / 60, 3),
                    "min"))
        print("best params:")
        print(best_params$par)
        best_params_mat[sim,] <- best_params$par
        print(paste("Total time:", round( (proc.time()[3] - t_init)  / 3600, 3),
                    "hrs"))

    }

    means <- colMeans(best_params_mat)
    medians <- apply(best_params_mat, 2, median)
    sd <- apply(best_params_mat, 2, sd)
    
    cluster_sizes_vec <- do.call("c",
                                 cluster_sizes_list)
    cluster_med <- median(cluster_sizes_vec)
    cluster_max <- max(cluster_sizes_vec)
    cluster_90 <- quantile(cluster_sizes_vec, probs = .90)
    names(cluster_90) <- NULL




    data_out <- list(best_params_mat = best_params_mat,
                     seed = my_seed,
                     init_params = inf_params,
                     gamma = gamma,
                     set = par_index,
                     par_set = params_df[par_index,],
                     cluster_sizes = c("med" = cluster_med,
                                       "max" = cluster_max,
                                       "q90" = cluster_90))

    simulation_output[[par_index]] <- data_out


    current_time <- Sys.time()
    current_time <- gsub(" ", "_", current_time)
    current_time <- gsub(":", ".", current_time)
    fn_base <- paste0("single_output_outsider_", "set_",
                      par_index,
                      current_time, ".RDS")

    saveRDS(data_out, fn_base)

}


current_time <- Sys.time()
current_time <- gsub(" ", "_", current_time)
current_time <- gsub(":", ".", current_time)
fn_base <- paste0("all_output_outsider_", current_time, ".RDS")
    saveRDS(simulation_output, fn_base)
