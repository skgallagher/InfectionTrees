#' Calculate the estimated general loglikelihood
#'
#' @param inf_params vector of p parameters
#' @param observed_data data frame with the following columns
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
#' @param sampled_data data frame of samples that correspond to the data.  See details.  This is the output of \code{general_cond_tree_sims()}.
#' @param return_neg default is TRUE.  Returns the negative loglike
#' @param cov_mat optional matrix of covariates corresponding to the sampled_data
#' @param cov_names covariate vector of length p which correspond in order to the betas
#' @param use_outsider_prob a separate parameter for the outsider infections.  Default is FALSE
#' @param return_clust_loglikes if TRUE, function returns the likelihood of every individual cluster
#' @return estimated average loglikelihood for the observed data
#' @details This is a specialized log likelihood function where we first estimate the average log likelihood of trees conditioned by their total size through sampling.  This is very much dependent on the values of sampled_data.
#' @export
general_loglike <- function(inf_params,
                            observed_data,
                            sampled_data,
                            return_neg = TRUE,
                            cov_mat = NULL,
                            cov_names = NULL,
                            use_outsider_prob = FALSE,
                            return_clust_loglikes = FALSE){
    if(use_outsider_prob){
        p0 <- inf_params[1]
        inf_params <- inf_params[-1]

    }
    if(is.null(cov_mat)){ ## have to format it.  slower
        cov_mat <- covariate_df_to_mat(sampled_data, cov_names)
        sampled_data$prob_inf <- 1 / (1 + exp(-(cov_mat %*% inf_params)))
    }

    cluster_id <- n_inf <- orig_id <- prob_inf <- NULL
    ## sampled_data$prob_inf <- 1 / (1 + exp(-(cov_mat %*% inf_params)))
    my_prob_inf <-  1 / (1 + exp(-(cov_mat %*% inf_params)))
    if(use_outsider_prob){
        my_prob_inf <- ifelse(sampled_data$gen == 1,
                           p0,
                           my_prob_inf)
    }

    sampled_data <- sampled_data[, prob_inf := my_prob_inf]
    cluster_id <- NULL
    ## Trying out data.table
    like_df <- sampled_data[,
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



#' Calculate the estimated general loglikelihood
#'
#' @param inf_params vector of p parameters
#' @param obs_summarized_data data frame with the following columns
#' \describe{
#' \item{cluster_size}{size of the cluster}
#' \item{freq}{frequency that cluster size occurred in the data}
#' }
#' @param sampled_data data frame of samples that correspond to the data.  See details.  This is the output of \code{sample_general_cond_trees()}.
#' @param return_neg default is TRUE.  Returns the negative loglike
#' @param cov_mat optional matrix of covariates corresponding to the sampled_data
#' @param cov_names covariate vector of length p which correspond in order to the betas
#' @return estimated average loglikelihood for the observed data
#' @details This is a specialized log likelihood function where we first estimate the average log likelihood of trees conditioned by their total size through sampling.  This is very much dependent on the values of sampled_data.
#' @export
general_loglike_summarized <- function(inf_params,
                            obs_summarized_data,
                            sampled_data,
                            return_neg = TRUE,
                            cov_mat = NULL,
                            cov_names = NULL){


    if(is.null(cov_mat)){ ## have to format it.  slower
        cov_mat <- covariate_df_to_mat(sampled_data, cov_names)
        sampled_data$prob_inf <- 1 / (1 + exp(-(cov_mat %*% inf_params)))
    }

    sampled_data$prob_inf <- 1 / (1 + exp(-(cov_mat %*% inf_params)))
    avg_loglike_sample_df <- general_loglike_sampled_data(sampled_data)


    ## Probably can take out join and replace with cbind if we're careful..
    joined_df <- dplyr::left_join(obs_summarized_data,
                                  avg_loglike_sample_df,
                                  by = "cluster_size")
    loglike <- joined_df %>%
        dplyr::mutate(cond_loglike = .data$freq *
                          .data$avg_loglike) %>%
        dplyr::select(.data$cond_loglike) %>%
        sum()

    if(is.na(loglike)) stop("NA loglike value")




    if(return_neg){
        loglike <- -loglike
    }
    return(loglike)
}



#' Calculate the log likelihood for a set of data, averaged by n
#'
#' @param sampled_data data frame with either (the cluster_id, person_id, n_inf, cluster_size, prob_inf) or additionally with the covariate information
#' @return data frame of average loglikelihood for a cluster of a given size
general_loglike_sampled_data <- function(sampled_data){



    loglike_df <- sampled_data %>% dplyr::group_by(.data$cluster_id,
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
#' @param sampled_data data frame with cov_names
#' @param cov_names null or string of vector
#' @return matrix of dimension n x (p+1)
#' @export
covariate_df_to_mat <- function(sampled_data, cov_names){
    if(!is.null(cov_names)){
        stopifnot(length(cov_names) ==
                  sum(cov_names %in% colnames(sampled_data)))
        sample_covariates_df <- sampled_data %>%
            dplyr::select(tidyselect::all_of(cov_names))
        cov_mat <- as.matrix(sample_covariates_df)
        cov_mat <- cbind(rep(1, nrow(cov_mat)), cov_mat)
        colnames(cov_mat) <- NULL
        return(cov_mat)
    } else{
        stop("You must provide appropriate column names for sampled_data")

    }

}
