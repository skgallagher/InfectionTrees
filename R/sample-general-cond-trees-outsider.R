## SKG
## Outside sampling general features
## June 24, 2020

##



#' Sample trees with same features as original but permuted through the trees
#'
#' @param observed_data data.table frame with the following columns
#' \describe{
#' \item{cluster_id}{unique cluster ID}
#' \item{person_id}{unique identifier for every person}
#' \item{covariates}{covariates of the individuals}
#' }
#' @param B number of samples for each individual tree
#' @return data frame with the following columns
#' \describe{
#' \item{cluster_id}{unique cluster_id for the sampled cluster}
#' \item{orig_id}{original ID of the given cluster}
#' \item{variable}{columns from featureless tree sampling}
#' \item{person_id}{original identifier}
#' }
#' @details assumes an observed outside infection
#' @export
sample_gen_out_cond <- function(observed_data,
                                B){

    person_ids <- observed_data$person_id

    unique_clusters <- unique(observed_data$cluster_id)
    sampled_clusters <- lapply(unique_clusters, function(id){
        cluster_inds <- which(observed_data$cluster_id == id)
        cluster_person_ids <- observed_data$person_id[cluster_inds]
        sample_single_gen_out_cond(cluster_person_ids = cluster_person_ids,
                                   B = B,
                                   person_ids = person_ids,
                                   orig_id = id)

        })

    sampled_data_df <- do.call('rbind', sampled_clusters)
    return(sampled_data_df)



}

#' Make the sampled clusters for a single cluster
#'
#' @param cluster_person_ids ids of the people WITHIN the cluster
#' @param B number of times to sample a cluster
#' @param person_ids person_ids of all people
#' @param orig_id original ID
#' @return data frame of transmissions trees for the original cluster
sample_single_gen_out_cond <- function(cluster_person_ids,
                                       B,
                                       person_ids,
                                       orig_id){


    cluster_size <- length(cluster_person_ids) + 1
    trees <- sample_unif_trees_no_features(n_vec = cluster_size,
                                           B = B)

    trees$orig_id <- orig_id
    outsider_inds <- (0:(B - 1)) * cluster_size + 1

    outsider_sampled_inds <- sample(1:length(person_ids),
                                    size = B, replace = TRUE)
    permuted_sampled_inds <- as.numeric(sapply(1:B,
                                               function(x){
                                                   sample(1:length(cluster_person_ids))
                                               }))




    feature_id <- numeric(cluster_size * B)
    feature_id[outsider_inds] <- person_ids[outsider_sampled_inds]
    feature_id[-outsider_inds] <- cluster_person_ids[permuted_sampled_inds]
    trees$feature_id <- feature_id
    return(trees)



}


