## SKG
## Testing out model where transmission tree is determined by infection order

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

naive_loglike_wrapper <- function(x, inf_params, dt,
                                  covariate_names,
                                  cov_mat,
                                  beta_index,
                                  max_loglike,
                                  chi_stat = 1.92){
    inf_params[beta_index] <- x
    naive_loglike(inf_params, dt,
                          covariate_names, cov_mat,
                  return_neg = FALSE) + chi_stat - max_loglike

}

######################3
## PERMUTE relative time within clusters
set.seed(7272020)
beta1_vec <- numeric(10)
for(jj in 1:10){
    tb_sub <- tb_orig %>%
        group_by(cluster_id) %>%
        mutate(rel_time = sample(rel_time)) %>%
        ungroup()


    covariate_names <- c("smear", "hiv_neg_pos", "hiv_unk_pos")
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
                         upper = -bds
                         )
    t2 <- proc.time()[3] - t1
    print(paste("Optimization time:", round( t2 / 60, 3),
                "min"))

    ## Likelihood profile
    print("Likelihood profile")


    ci_mat <- matrix(0, ncol = 2,
                     nrow = length(best_params$par))

    ci_mat <- matrix(0, ncol = 2,
                     nrow = length(best_params$par))

    for(ii in 1:nrow(ci_mat)){
        max_loglike <- -best_params$value
        best_pars <- best_params$par
        lower <- uniroot(f = naive_loglike_wrapper, c(best_pars[ii] - 3, best_pars[ii]),
                         max_loglike = max_loglike,
                         inf_params = best_pars,
                         beta_index = ii,
                         dt = data.table::as.data.table(tb_sub),
                         covariate_names = covariate_names,
                         cov_mat = cov_mat,
                         chi_stat = 1.92)
        ##
        upper <- uniroot(f = naive_loglike_wrapper, c(best_pars[ii], best_pars[ii] + 3),
                         max_loglike = max_loglike,
                         inf_params = best_pars,
                         beta_index = ii,
                         dt = data.table::as.data.table(tb_sub),
                         covariate_names = covariate_names,
                         cov_mat = cov_mat,
                         chi_stat = 1.92)

        ci_mat[ii, 1] <- lower$root
        ci_mat[ii, 2] <- upper$root
    }
    est_pars <- cbind(best_pars, ci_mat)
    rownames(est_pars) <- c("Intercept",
                            covariate_names)
    colnames(est_pars) <- c("Mean", "Lower95", "Uppper95")
    beta1_vec[jj] <- est_pars[2,1]

    print(-best_params$value)
    print(est_pars, digits = 2)
}

#######################################3

######################3
## ORIGINAL
set.seed(7272020)
beta1_vec <- numeric(10)
for(jj in 1:1){
    tb_sub <- tb_orig 


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

    ## Likelihood profile
    print("Likelihood profile")


    ci_mat <- matrix(0, ncol = 2,
                     nrow = length(best_params$par))

    ci_mat <- matrix(0, ncol = 2,
                     nrow = length(best_params$par))

    ci_mat[,1] <- best_params$par - 1.96 * se
    ci_mat[,2] <- best_params$par + 1.96 * se
   
    est_pars <- cbind(best_pars, ci_mat, se)
    rownames(est_pars) <- c("Intercept",
                            covariate_names)
    colnames(est_pars) <- c("Mean", "Lower95", "Uppper95", "SE")
    beta1_vec[jj] <- est_pars[2,1]

    print(-best_params$value)
    print(est_pars, digits = 2)
}

saveRDS(est_pars, "naive-tab.RDS")
