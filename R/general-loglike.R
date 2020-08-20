#' Calculate the estimated general loglikelihood
#'
#' @param inf_params vector of p parameters
#' \describe{
#' \item{cluster_id}{unique cluster ID}
#' \item{person_id}{order of infection in the cluster}
#' \item{gen}{generation number (>=0)}
#' \item{inf_id}{ID of the infector}
#' \item{n_inf}{number of people infected by person}
#' \item{censored}{whether the cluster end was censored or not}
#' \item{cluster_size}{size of the cluster}
#' \item{covariates}{covariates of the individuals}
#' }
#' @param mc_trees data frame of samples that correspond to the data.  See details.
#'  This is the output of \code{general_cond_tree_sims()}.
#' @param return_neg default is TRUE.  Returns the negative loglike
#' @param cov_mat optional matrix of covariates corresponding to the mc_trees
#' @param cov_names covariate vector of length p which correspond in order to the betas
#' @param multiple_outside_transmissions logical indicating whether to use multiple outside method to compute likelihood.
#' @param use_outsider_prob a separate parameter for the outsider infections.  Default is FALSE
#' @param return_clust_loglikes if TRUE, function returns the likelihood
#'  of every individual cluster
#' @param messages Should we print messages.
#' @return estimated average loglikelihood for the observed data
#' @details This is a specialized log likelihood function where we first estimate
#' the average log likelihood of trees conditioned by their total size through sampling.
#'  This is very much dependent on the values of mc_trees.
#' @details The base likelihood for a single tree is given by
#' \deqn{L(T) = \prod_{i=1}^n (1-p_i)p_i^{N_i}}
#' where \eqn{p_i} is the probability of transmission for individual \eqn{i} and \eqn{N_i} is the number of individuals infected by individual \eqn{i}.  The approximate average likelihood for a given cluster $C_m$ is then
#' \deqn{\bar{L}_K(C_m) = \frac{1}{K}\sum_{k=1}^K L(T_k)}
#' where \eqn{K} is the number of Monte Carlo transmission tree samples \eqn{T_k} for cluster \eqn{C_m}.  Finally, the log likelihood is the sum of the log likelihoods for each cluster,
#' \deqn{\ell(C_1, \dots, C_M) = \sum_{m=1}^M log(\bar{L}_K(C_m))}
#' For the multiple outside transmissions model, the above likelihood calculation is changed only for a single tree (and the transmission trees are of a different form).  The likelihood \eqn{L_O(T)} is
#' \deqn{L_O(T) = (1-p_1)p_1^{N_1-1}\prod_{i=2}(1-p_i)p_i^{N_i}}
#' because we condition on the outsider having at least one successful infection.
#' @export
general_loglike <- function(inf_params,
                            mc_trees,
                            return_neg = TRUE,
                            cov_mat = NULL,
                            cov_names = NULL,
                            multiple_outside_transmissions = FALSE,
                            use_outsider_prob = FALSE,
                            return_clust_loglikes = FALSE,
                            messages = FALSE){
    if(use_outsider_prob){
        p0 <- inf_params[1]
        inf_params <- inf_params[-1]

    }
    if(is.null(cov_mat)){ ## have to format it.  slower
        cov_mat <- covariate_df_to_mat(mc_trees, cov_names)
        mc_trees$prob_inf <- 1 / (1 + exp(-(cov_mat %*% inf_params)))
    }
    if(multiple_outside_transmissions){  ## condition on outsider (gen 1) having one success
        mc_trees$n_inf <- with(mc_trees,
                                    ifelse(gen == 1 & n_inf > 0,
                                           n_inf - 1, n_inf))
    }

    cluster_id <- n_inf <- orig_id <- prob_inf <- NULL
    my_prob_inf <-  1 / (1 + exp(-(cov_mat %*% inf_params)))
    if(use_outsider_prob){
        my_prob_inf <- ifelse(mc_trees$gen == 1,
                           p0,
                           my_prob_inf)
    }

     ## Trying out data.table
    if(!("data.table" %in% class(mc_trees))){
        if(messages){
            print("Converting 'mc_trees' to data.table format")
        }
        mc_trees <- data.table::as.data.table(mc_trees)
    }
    mc_trees <- mc_trees[, prob_inf := my_prob_inf]
    cluster_id <- NULL

    like_df <- mc_trees[,
                            .(like = general_cluster_like.dt(prob_inf, n_inf)),
                            by = .(orig_id, cluster_id)]
    if(return_clust_loglikes){
        return(like_df)
    }
    avg_loglike_df <- like_df[, .(avg_loglike = log(mean(like))),
                              by = .(orig_id)]



    loglike <- sum(avg_loglike_df$avg_loglike)

    if(return_neg){
        loglike <- -loglike
    }

    if(is.na(loglike) | is.nan(loglike) |
       is.infinite(loglike)) {
        stop("loglike is NA/NAN/infinite")
    }

    return(loglike)


}




#' Calculate the log likelihood for a set of data, averaged by n
#'
#' @param mc_trees data frame with either (the cluster_id, person_id, n_inf, cluster_size, prob_inf) or additionally with the covariate information
#' @return data frame of average loglikelihood for a cluster of a given size
general_loglike_mc_trees <- function(mc_trees){



    loglike_df <- mc_trees %>% dplyr::group_by(.data$cluster_id,
                                                   .data$cluster_size) %>%
        tidyr::nest() %>%
        dplyr::mutate(like =  purrr::map(.data$data,
                                         general_cluster_like)) %>%
        dplyr::select(-.data$data) %>%
        tidyr::unnest(cols= like) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$cluster_size) %>%
        dplyr::summarize(avg_loglike = log(mean(like))) %>%
        dplyr::select(.data$cluster_size, .data$avg_loglike)

    return(loglike_df)


}

#' Get the likelihood for a single cluster
#'
#' @param data data frame with the following columns
#' \describe{
#' \item{prob_inf}{probability of onward transmission}
#' \item{n_inf}{number infected by this person}
#' }
#' @return a single number, the likelihood
general_cluster_like <- function(data){

    like <- prod((1 - data$prob_inf) * (data$prob_inf^data$n_inf))

    ## Old version
    ## like <- data %>%
    ##     dplyr::mutate(like =
    ##                       (1-.data$prob_inf) *
    ##                       (.data$prob_inf)^.data$n_inf) %>%
    ##     dplyr::select(.data$like) %>% prod()
    return(like)


}

#' Get the likelihood for a single cluster
#'

#' @param prob_inf probability of onward transmission
#' @param n_inf number infected by this person
#' @return a single number, the likelihood
#' @details (specifically for data.table)
general_cluster_like.dt <- function(prob_inf, n_inf){

    like <- prod((1 - prob_inf) * (prob_inf^n_inf))

    return(like)


}



#' Take data frame and turn it into matrix for logistic regression
#'
#' @param mc_trees data frame with cov_names
#' @param cov_names null or string of vector
#' @return matrix of dimension n x (p+1)
#' @export
covariate_df_to_mat <- function(mc_trees, cov_names){
    if(!is.null(cov_names)){
        stopifnot(length(cov_names) ==
                  sum(cov_names %in% colnames(mc_trees)))
        sample_covariates_df <- mc_trees %>%
            dplyr::select(tidyselect::all_of(cov_names))
        cov_mat <- as.matrix(sample_covariates_df)
        cov_mat <- cbind(rep(1, nrow(cov_mat)), cov_mat)
        colnames(cov_mat) <- NULL
        return(cov_mat)
    } else{
        stop("You must provide appropriate column names for mc_trees")

    }

}
