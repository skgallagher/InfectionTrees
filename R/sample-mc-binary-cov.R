## SKG
## August 17, 2020
## Samplng binary covariate trees


#' Sample binary covariate transmission trees
#'
#' @param B number of MC samples
#' @param observed_cluster_summaries data frame with the following columns
#' \describe{
#' \item{freq}{frequency of the following cluster}
#' \item{cluster_size}{size of the cluster}
#' \item{x_pos}{number of positive x (or 1)}
#' \item{x_neg}{number of negative x (or 0)}
#' }
#' @param multiple_outside_transmissions logical indicating whether to
#' sample from the multiple outside transmissions model.  Default is FALSE.
#' @return  a data frame of sampled transmission trees
#'  with the following columns (the last one is only needed for the multiple
#'  outsider model)
#'  \describe{
#'  \item{freq}{frequency of cluster with the following characteristics}
#'  \item{cluster_size}{size of the cluster}
#'  \item{x_pos}{number of x positives in the cluster}
#'  \item{x_neg}{number of x negatives in the cluster}
#'  \item{mc_freq}{frequency of MC trees with the following}
#'  \item{x_pos_trans}{number of transmissions from a positive individual in the tree}
#'  \item{x_neg_trans}{number of transmissions from a negative individual in the tree}
#'  \item{root_node_sign}{whether root node is pos or neg}
#'  }
sample_mc_binary_cov <- function(B,
                                 observed_cluster_summaries,
                                 multiple_outside_transmissions = FALSE){

    ## Get total distribution of positive vs negative
    total_x_pos <- sum(observed_cluster_summaries$x_pos *
                       observed_cluster_summaries$freq)
    total_x_neg <- sum(observed_cluster_summaries$x_neg *
                       observed_cluster_summaries$freq)
    empirical_p_pos <- total_x_pos / (total_x_neg + total_x_pos)
    
    tree_summary_list <- vector(mode = "list",
                                length = nrow(observed_cluster_summaries))
    for(ii in 1:nrow(observed_cluster_summaries)){
        x_pos <- observed_cluster_summaries$x_pos[ii]
        x_neg <- observed_cluster_summaries$x_neg[ii]
        freq <- observed_cluster_summaries$freq[ii]
        if(multiple_outside_transmissions){
            outsider_pos_freq <- stats::rbinom(n = 1, size = B, prob = empirical_p_pos)
            outsider_neg_freq <- B - outsider_pos_freq
            trees_pos <- sample_mc_binary_cov_inner(x_pos,
                                                     x_neg,
                                                     B = outsider_pos_freq,
                                                    root_node = 1)
            trees_neg <- sample_mc_binary_cov_inner(x_pos,
                                                     x_neg,
                                                     B = outsider_neg_freq,
                                                    root_node = 0)
            trees <- rbind(trees_pos,
                                  trees_neg)
            
            
        } else{
            trees <- sample_mc_binary_cov_inner(x_pos,
                                                x_neg,
                                                B = B)
        }
        tree_summary <- summarize_binary_cov_trees(trees,
                                                   multiple_outside_transmissions =
                                                       multiple_outside_transmissions)
        tree_summary$freq <- freq

        tree_summary_list[[ii]] <- tree_summary
    }
    mc_summary <- dplyr::bind_rows(tree_summary_list)
    return(mc_summary)
    

}
                                 

#' Innfer function to sample MC trees with binary covariates
#'
#' @param x_pos number of positive individuals in cluster
#' @param x_neg number of negative individuals in cluster
#' @param B number of trees to sample
#' @param root_node if NULL, use the base model.  If 1, ensure that the root node is positive.  If 0, ensure that the root node is negative.
#' @return a data frame with the following columns
#' \describe{
#' \item{freq}{frequency of the observed data}
#' \item{cluster_size}{total size of the cluster}
#' \item{x_pos}{number of individuals in the cluster with the feature of interest =1}
#' \item{x_neg}{number of individuals in the cluster with the feature of interest  =0}
#' \item{mc_freq}{number of MC frequency of clusters}
#' \item{x_pos_trans}{number of transmissions from a positive}
#' \item{x_neg_trans}{number of transmissions from a negative}
#' \item{root_node}{only for the multiple outsider.  Is 1 if the outsider is positive and 0 if it is negative}
#' }
sample_mc_binary_cov_inner <- function(x_pos,
                                 x_neg,
                                 B,
                                 root_node = NULL){
    n <- x_pos + x_neg
    type <- "binary_cov"
    if(!is.null(root_node)){
        n <- n + 1
        type <- "binary_cov_out"
    }
    tree_df <- sample_unif_trees_no_features(n, B)
    params_list <- list(x_pos = x_pos,
                        x_neg = x_neg,
                        root_node = root_node)
    
    

    tree_df <- draw_features(tree_df = tree_df,
                             feature_type = type,
                             params_list = params_list)
    return(tree_df)


}


#' Condense trees into a smaller format
#'
#' @param trees data frame with the following columns
#' \describe{
#' \item{cluster_id}{cluster ID}
#' \item{n_inf}{number of infected by this individual}
#' \item{x}{binary covariate value (0/1)}
#' \item{gen}{generation number}
#' }
#' @param multiple_outside_transmissions default is FALSE
#' @return data frame with the following columns
#' \describe{
#' \item{cluster_size}{size of cluster}
#' \item{x_pos}{number of positive x}
#' \item{x_neg}{number of negative x}
#' \item{x_pos_trans}{number of transmissions from positives}
#' \item{x_neg_trans}{number of transmissions from negatives}
#' \item{root_node}{value of the root node}
#' }
summarize_binary_cov_trees <- function(trees,
                                   multiple_outside_transmissions = FALSE){
    summary_trees <- trees %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::summarize(cluster_size = dplyr::n(),
                         x_pos = sum(.data$x == 1),
                         x_neg = sum(.data$x == 0),
                         x_pos_trans = sum((.data$x == 1) * .data$n_inf),
                         x_neg_trans = sum((.data$x == 0) * .data$n_inf),
                         root_node = .data$x[.data$gen == 1]) %>%
        dplyr::ungroup()
    if(!multiple_outside_transmissions){
        summary_trees <- summary_trees %>%
            dplyr::group_by(.data$x_pos,
                     .data$x_neg,
                     .data$x_pos_trans,
                     .data$x_neg_trans) %>%
            dplyr::summarize(mc_freq = dplyr::n()) %>%
            dplyr::ungroup()
                     
    }else{
        summary_trees <- summary_trees %>%
            dplyr::group_by(.data$x_pos,
                            .data$x_neg,
                            .data$x_pos_trans,
                            .data$x_neg_trans,
                            .data$root_node) %>%
            dplyr::summarize(mc_freq = dplyr::n()) %>%
            dplyr::ungroup()
    }
    return(summary_trees)
                                           

}
