## SKG
## Testing out model where transmission tree is determined by infection order
## Updated June 23, 2021 to add a bootstrap SE

## Load libraries and data
rep_from_lib <- TRUE

if(rep_from_lib){
    library(InfectionTrees)
} else {
    devtools::load_all()
}
library(tidyverse)
data(tb_clean)

my_seed <- 6172021
set.seed(my_seed)

t0 <- proc.time()[3]

## Format the data for all the variables

## I'm filling in the Unknown with a negative smear
## Should later check to see if results are changed with positive
tb_sub <- tb_clean %>%
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
tb_orig <- tb_sub




naive_loglike <- function(inf_params, dt,
                          covariate_names, cov_mat,
                          return_neg = TRUE){

    dt <- dt[, prob_inf :=
                   1 / (1 + exp(-(cov_mat %*% inf_params)))]
    ## Trying out data.table
    loglike_df <- dt[,
                     .(loglike = naive_cluster_loglike(prob_inf, rel_time)),
                            by = .(cluster_id)]
    loglike <- sum(loglike_df$loglike)
    if(return_neg) loglike <- -loglike
    return(loglike)
    


}

# cluster contains column datecoll which is the date of collection
naive_cluster_loglike <- function(prob_inf, datecoll){
    last_ind <- which.max(datecoll)
    probs <- prob_inf[-last_ind]
    loglike <- sum(log(probs)) + sum(log(1 - prob_inf))
    return(loglike)

}



######################3
## ORIGINAL
set.seed(7272021)
beta1_vec <- numeric(10)
B <- 100


beta_mat <- matrix(0, nrow = B, ncol = 5)
cov_mat <- matrix(0, nrow = B, ncol = 5)
for(jj in 1:B){

    ## boot strap orig data
    tb_sub <- bootstrap_clusters(tb_orig)
    
   

    covariate_names <- c("smear", "hiv_neg_pos", "hiv_unk_pos", "rel_time")
    cov_mat <- covariate_df_to_mat(tb_sub,
                                   cov_names = covariate_names)
    init_params <- rep(0, length(covariate_names) + 1)
    bds <- rep(-5, length(covariate_names) + 1)
    t1 <- proc.time()[3]
    best_params <- optim(par = init_params,
                         fn = naive_loglike,
                         dt = data.table::as.data.table(tb_sub),
                         return_neg = TRUE,
                         cov_mat = cov_mat,
                         covariate_names = covariate_names,
                         method = "L-BFGS-B",
                         lower = bds,
                         upper = -bds,
                         hessian = TRUE
                         )
    t2 <- proc.time()[3] - t1
    print(paste("Optimization time:", round( t2 / 60, 3),
                "min"))
    se <- sqrt(diag(solve(best_params$hessian)))
    beta_mat[jj,] <- best_params$par
    lower <- best_params$par - 1.96 * se
    upper <- best_params$par + 1.96 * se
    cov_mat[jj,] <- ifelse(best_params$par >= lower &
                      best_params$par <= upper, 1, 0)
 
    print(-best_params$value)
    print(best_params$par, digits = 2)
}

se_boot <- apply(beta_mat, 2, sd)

saveRDS(list(beta_mat = beta_mat,
           se_boot = se_boot), "naive-tab-boot.RDS")
