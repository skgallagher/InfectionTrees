#!/usr/bin/env Rscript


## June 15, 2021
## SKG
## Rerun the simulations
## This is for the base model with a binary covariate
## We are using alpha values of .5, .7, .1
## bet0 values of -2.5, -1.5, and -2.75
## and beta1 values of -.5, 0, -.25
## There are 15 total sets (so not all combinations are used)

rep_from_lib <- TRUE

if(rep_from_lib){
    library(InfectionTrees)
} else {
    devtools::load_all()
}



library(dplyr)

set.seed(6172021)

## Simulations to try
M <- 100 # number of simulations for each parameter configuration
B <- 5000 # number of MC samples
K <- 1000 # number of clusters in a data set


par_df1 <- data.frame(
                     M = M,
                     B = B,
                     beta_0 = c(-2.5, -2.5, -2.5,
                         -1.5, -1.5, -1.5),
                     beta_1 = c(.5, 0, -.5,
                                1, 0, -1),
                     gamma = .5)

par_df2 <- data.frame(
                     M = M,
                     B = B,
                     beta_0 = c(-2.5, -2.5, -2.5,
                         -1.5, -1.5, -1.5),
                     beta_1 = c(.5, 0, -.5,
                                .5, 0, -.5),
                     gamma = .7)

par_df3 <-  data.frame(
                     M = M,
                     B = B,
                     beta_0 = c(-2.75, -2.75, -2.75),
                     beta_1 = c(2.25, 0, -.25),
                     gamma = .1)
###########################################################################

## data frame of simulation parameter configurations
par_df <- dplyr::bind_rows(par_df1, par_df2, par_df3)

out_list <- vector(mode = "list", length = nrow(par_df))

#########################################################################

## Generate a small set of pre-sampled trees to get started
binary_cluster_summaries <- data.frame(cluster_size = c(1, 1, 2, 2, 2),
                                       x_pos = c(0, 1, 0, 1, 2),
                                       x_neg = c(1, 0, 2, 1, 0),
                                       freq = 1)
mc_trees <- sample_mc_binary_cov(B = B,
                                 observed_cluster_summaries =
                                     binary_cluster_summaries) %>%
    dplyr::select(-freq)



## Run the simulations
total_time <- proc.time()[3]
for(zz in 1:nrow(par_df)){
  t_params <- proc.time()[3]
  beta_0 <- par_df$beta_0[zz]
  beta_1 <- par_df$beta_1[zz]
  gamma <- par_df$gamma[zz]
  cat(sprintf("Parameter set %d \n", zz))
  cat(sprintf("beta_0 = %.2f, beta_1 = %.2f, gamma = %.2f\n",
                beta_0, beta_1, gamma))


  ## For a given set of parameters
  best_params_mat <- matrix(0, nrow = M, ncol = 2)
  coverage_mat <- matrix(0, nrow = M, ncol = 4)
  colnames(coverage_mat) <- c("lp_beta0", "lp_beta1",
                              "fi_beta0", "fi_beta1")
  cluster_sizes <- numeric(K * M) # store all the cluster sizes from a simulated set of
  ## observed data
  max_cluster_sizes <- matrix(M)
  for(mm in 1:M){
    t_sims <- proc.time()[3]
    cat( sprintf("Simulation %d \n", mm))
    ## Simulate a branching process data set
    ## Distribution of xs
    sample_covariates_df <- data.frame(x = c(rep(0, gamma * 10),
                                             rep(1, (1-gamma) * 10)))


    data <- simulate_bp(K = K ,
                                      inf_params = c(beta_0, beta_1),
                                      sample_covariates_df = sample_covariates_df,
                                      covariate_names = "x",
                                      max_size = 55)
    ## Summarize the data set
    binary_cluster_summaries <- summarize_binary_clusters(data)
    cluster_sizes[(mm-1) *K +  (1:K)] <- rep(binary_cluster_summaries$cluster_size,
                                         times = binary_cluster_summaries$freq)
    
    max_cluster_sizes[mm] <- max(binary_cluster_summaries$cluster_size)
    cat(sprintf("max cluster size is %d \n", max_cluster_sizes[mm]))


    ## ##################################3
    ## SAMPLE (ADDITIONAL) MC TREES
    ## #######################################
    size_npos <- with(binary_cluster_summaries,
                      paste0(cluster_size, "-", x_pos))
    sampled_size_npos <- with(mc_trees,
                              paste0(cluster_size, "-", x_pos))
    missing_sample_ids <- which(!(size_npos %in% sampled_size_npos))
    if(length(missing_sample_ids) > 0){
        new_cluster_summaries <- binary_cluster_summaries[missing_sample_ids,]
        print("Sampling the following trees")
        print(size_npos[missing_sample_ids])
        new_mc_trees <- sample_mc_binary_cov(B = B,
                                             observed_cluster_summaries =
                                                 new_cluster_summaries) %>%
            dplyr::select(-freq)
        mc_trees <- bind_rows(mc_trees, new_mc_trees) %>%
            arrange(cluster_size, x_pos)
    }
    #######################################
    ## OPTIMIZE
    ######################################
    inf_params_0 <- c(0,0)

    best_params <- optim(inf_params_0, fn = bp_loglike_binary_cov,
                         obs_data_summary = binary_cluster_summaries,
                        mc_samples_summary = mc_trees,
                                    return_neg = TRUE,
                       lower = c(-5, -5),
                       upper = c(0, 5),
                       method = "L-BFGS-B",
                       hessian = TRUE)
    best_params_mat[mm,] <- best_params$par
    ## Fisher information (approximately)
    se <- sqrt(diag(solve(best_params$hessian)))
    lower <- best_params$par - 1.96 * se
    upper <- best_params$par + 1.96 * se
    coverage_mat[mm, 3:4] <- ifelse(lower  < c(beta_0, beta_1) &
                                 upper  > c(beta_0, beta_1),
                                 1, 0)

    print("best params")
    print(best_params$par)

    ## Likelihood profiles
    ## lp_ests <- binary_likelihood_profs(best_pars = best_params$par,
    ##                                    max_loglike = -best_params$value,
    ##                                    obs_data_summary =
    ##                                        binary_cluster_summaries,
    ##                                    mc_samples_summary = mc_trees,
    ##                                    alpha = .05)

#    print("95% CI from LP")
  #  print(lp_ests)
    print("95% CI from FI")
    print(cbind(lower, upper))

    ## coverage_mat[mm, 1:2] <- ifelse(lp_ests[,2] < c(beta_0, beta_1) &
    ##                              lp_ests[,3] > c(beta_0, beta_1),
    ##                              1, 0)


    print("Data simulation time (hrs)")
    print((proc.time()[3] - t_sims) / 3600, digits = 2)

    print("Cumulative time (hrs)")
    print((proc.time()[3] - total_time) / 3600, digits = 2)

  }

  print("Params simulation time (hrs)")
  print((proc.time()[3] - t_params) / 3600, digits = 2)


  cat(sprintf("beta_0 = %.2f, beta_1 = %.2f, gamma = %.2f\n ",
          beta_0, beta_1, gamma))
  print("Mean estimates")
  print(colMeans(best_params_mat))
   print("Median estimates")
 print( apply(best_params_mat, 2, median))
  print("2.5Q")
print(  apply(best_params_mat, 2, quantile, probs = .025))
   print("2.5Q")
  print(  apply(best_params_mat, 2, quantile, probs = .975))

  print("Coverage beta0, beta1")
  print(colMeans(coverage_mat))

  print("Cluster sizes:")
  cluster_size_summary <- quantile(cluster_sizes, probs = c(.5, .9, 1))
  print(cluster_size_summary)

  sims_list <- list(best_params_mat = best_params_mat,
                      set = zz,
                    max_cluster_sizes = max_cluster_sizes,
                    coverage = colMeans(coverage_mat),
                    cluster_size_summary = cluster_size_summary)

  out_list[[zz]] <- sims_list

}


current_time <- Sys.time()
current_time <- gsub(" ", "_", current_time)
current_time <- gsub(":", ".", current_time)
fn <- paste0("simulation_results_", current_time, ".RDS")


saveRDS(out_list, fn)
