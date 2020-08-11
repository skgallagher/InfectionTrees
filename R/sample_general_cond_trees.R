## SKG
## May 28, 2020
## A more general conditional tree sampler


#' Sample trees with same features as original but permuted through the trees
#'
#' @param observed_data data frame with the following columns
#' \describe{
#' \item{cluster_id}{unique cluster ID}
#' \item{covariates}{covariates of the individuals}
#' }
#' @param B number of trees to sample
#' @param covariate_names string names of covariates
#' @param B number of samples for each individual tree
sample_mc_trees <- function(observed_data,
                                       B = 100,
                                       covariate_names = "x"){
    if("id" %in% colnames(observed_data)){
        stop("Please rename the 'id' column")
    }
    ## this is not a pretty dplyr function
    n_clusters <- length(unique(observed_data$cluster_id))
    sampled_tree_list <- vector(mode = "list", n_clusters)
    for(ii in 1:length(unique(observed_data$cluster_id))){

        id <- unique(observed_data$cluster_id)[ii]
        observed_cluster <- observed_data %>%
            dplyr::filter(cluster_id == id) %>%
            dplyr::select(dplyr::contains(covariate_names))
        n <- nrow(observed_cluster)
        params_list <- list(covariate_df = observed_cluster,
                            weights = NULL,
                            replace = FALSE)
        sampled_clusters <- sample_general_cond_trees(n_vec = n,
                                                      B = B,
                                                      feature_type = "empirical",
                                                      params_list = params_list)
        sampled_clusters$orig_id <- id
        sampled_tree_list[[ii]] <- sampled_clusters


    }
    sampled_tree_df <- dplyr::bind_rows(sampled_tree_list)
    return(sampled_tree_df)



}


#' Sample general conditional trees
#'
#' @param n_vec vector of size of clusters
#' @param B number of sampled trees for each size
#' @param feature_type method of how we sample covariates.  Currently only "single_gauss", "empirical" are supported
#' @param params_list of additional parameters
#' @return data frame of trees
sample_general_cond_trees <- function(n_vec,
                                      B = 100,
                                      feature_type = "single_gauss",
                                      params_list = list(mean = 0,
                                                         sd = .5)){
    stopifnot(feature_type == "single_gauss" |
              feature_type == "empirical" |
              feature_type == "permute")
    if(feature_type == "single_gauss"){
        stopifnot(all(c("mean", "sd") %in% names(params_list)))
    } else if(feature_type == "empirical"){
        stopifnot(all(c("weights", "covariate_df") %in% names(params_list)))
    }

    ## Get the uniform trees
    tree_df <- sample_unif_trees_no_features(n_vec, B)

    ## Randomly draw features for the nodes on the tree
    tree_df <- draw_features(tree_df,
                             feature_type,
                             params_list)
    return(tree_df)

}


#' Sample trees
#'
#' @param n_vec vector of cluster sizes
#' @param B number of samples
#' @return data frame of sampled trees (could be very large)
sample_unif_trees_no_features <- function(n_vec, B){

    tree_list <- vector(mode = "list", length = length(n_vec))

    for(ii in 1:length(n_vec)){
        tree_list[[ii]] <- general_tree_sampler(n = n_vec[ii],
                                               B = B)
    }
    trees <- dplyr::bind_rows(tree_list)
    return(trees)

}

#' Generate trees of certain size
#'
#' @param n size of cluster
#' @param B number of times to repeat
#' @return data frame of trees
general_tree_sampler <- function(n, B){

    tree_list <- vector(mode = "list",
                        length = B)

    ## Sample generation sizes
    if(n == 1){ ## trivial case
        df <- data.frame(cluster_size = rep(1, B),
                         gen = 1,
                         n_in_gen = 1,
                         id = "1-1",
                         n_inf = 0,
                         cluster_id = paste("n1", "g1", 1:B,
                                            sep = "-"),
                         stringsAsFactors = FALSE)
        return(df)
    }

    ## Sample generation sizes
    ## only have 1 generation if there is only one person, otherwise there are at least 2
    g_vec <- 2 + stats::rbinom(n = B, size = n-2, prob = .5) # number of generations where first generation is fixed size

    g_tab <- table(g_vec)
    unique_g <- sort(unique(g_vec))
     for(ii in 1:length(unique_g)){
        ## Sample unique permutations of generation sizes given g
        g <- unique_g[ii]
        perm_mat <- sample_unique_perms(g = g, n = n,
                                        B = as.numeric(g_tab[ii]))
        tree_list[[ii]] <- sample_trees_fixed_g(perm_mat)  # generate whole tree

     }
    trees <- dplyr::bind_rows(tree_list)


    return(trees)

}

#' Sample trees with given generation sizes
#'
#' @param perm_mat g x B matrix where each column is a unique permutation of generation sizes
#' @return data frame with the following columns
#' \describe{
#' \item{cluster_id}{unique cluster ID}
#' \item{id}{person ID}
#' \item{gen}{generation number}
#' \item{n_inf}{number infected by this person}
#' }
sample_trees_fixed_g <- function(perm_mat){

    n <- sum(perm_mat[,1])
    g <- nrow(perm_mat)
    tree_list <- vector(mode = "list", length = nrow(perm_mat))
    for(ii in 1:ncol(perm_mat)){
        perm <- perm_mat[,ii]
        tree <- sample_general_tree_perm(perm)
        tree$cluster_id <- paste0("n", n, "-",
                                  "g", g, "-", ii)
        tree_list[[ii]] <- tree
    }
    tree_df <- do.call('rbind', tree_list) ## is this faster than bind_rows?
    return(tree_df)

}

#' Sample the actual tree for given permutation
#'
#' @param gen_sizes vector of generation sizes
#' @return sampled tree data framegiven the generation sizes
#' \describe{
#' \item{cluster_size}{cluster_size}
#' \item{id}{person ID}
#' \item{gen}{generation number}
#' \item{n_inf}{number infected}
#' }
sample_general_tree_perm <- function(gen_sizes){
    tree <- sample_connections(gen_sizes)
    tree$n_inf <- sapply(1:nrow(tree), function(ii){
        inds <- which(tree$inf_id == tree$id[ii])
        if(length(inds) == 0){
            return(0)
        } else{
            return(length(inds))
        }
    })
    n <- sum(gen_sizes)
    tree$cluster_size <- n

    return(tree)
}





#' Draw the features for the tree people
#'
#' @param tree_df data frame where each row is an individual
#' @param feature_type currently or supports "single_gauss" "empirical"
#' @param params_list additional parameters to pass
#' @return updated data frame with covariate features
draw_features <- function(tree_df,
                          feature_type,
                          params_list){
    n <- nrow(tree_df)

    if(feature_type == "single_gauss"){
        tree_df$x <- stats::rnorm(n, mean = params_list$mean,
                           sd = params_list$sd)

    } else if(feature_type == "empirical"){
        if(is.null(params_list$replace)){
            do_permute <- FALSE
        } else if(!params_list$replace) {
            do_permute <- TRUE
        } else{
            do_permute <- FALSE
        }
        weights <- params_list$weights
        cov_df <- params_list$covariate_df
        if(do_permute){
            n_groups <- length(unique(tree_df$cluster_id))
            sampled_inds <- as.numeric(sapply(1:n_groups,
                                   function(x) sample(1:nrow(cov_df))))

        } else if(is.null(weights)){ ## every row with equal prob
            sampled_inds <- sample(1:nrow(cov_df), size = n,
                                   replace = TRUE)
        } else{
            sampled_inds <- sample(1:nrow(cov_df), size = n,
                                   replace = TRUE,
                                   prob = weights / sum(weights))
        }
        tree_df <- cbind(tree_df, cov_df[sampled_inds,, drop = FALSE])

    } else{
        stop("No other options than 'single_gauss' and 'empirical' are supported")
    }
    rownames(tree_df) <- NULL
    return(tree_df)

}

