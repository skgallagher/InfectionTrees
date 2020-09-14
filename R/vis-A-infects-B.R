## SKG
## 9/14/20
## Utility functions to visualize the MC trees
## Specifically make a tidy data frame of prob A -> B


## TODO
## NEED to join people by prob of transmission,
## Identify people by prob of transmission
## factor by sort of prob_trans

## need to make MATRIX based on original_person_id and not id



#' Looks at most likely transmissions from A to B for a cluster
#'
#' @param my_id ID of cluster to analyze
#' @param par_ests parameter estimates (must be named and
#' have corresponding columns in \code{mc_trees})
#' @param mc_trees output from \code{\link{sample_mc_trees}}
#' @param output either \code{"tidy"} or \code{"matrix"}
#' @param multiple_outside_transmissions logical indicating which model to calculate likelihood
#' @return either a \eqn{n \times n}  matrix where \eqn{n} is the number of individuals
#'  in the cluster and entry
#' \eqn{i,j} is the probability that individual A infects individual B
#' or a data frame with the following columns
#' \describe{
#' \item{infector}{ID of the infector}
#' }
mc_trees_to_A_infects_B <- function(my_id,
                                    par_ests,
                                    mc_trees,
                                    output = "tidy",
                                    multiple_outside_transmissions = FALSE){


  ## Filter to original cluster id
  sub_data <- mc_trees %>%
    dplyr::filter(.data$orig_id == my_id)

  covariate_names <- names(par_ests)[-1]
  ## get prob of onward infection
  sub_data$prob_trans <- calculate_transmission_prob(data = sub_data,
                                                inf_params = par_ests,
                                                covariate_names = covariate_names)
  sub_data$feature_id <- factor(sort(sub_data$prob_trans))

  ## calculate like per cluster
  like_df <- general_loglike(inf_params = par_ests,
                             mc_trees = sub_data,
                             cov_names = covariate_names,
                             return_clust_loglikes = TRUE,
                             multiple_outside_transmissions =
                               multiple_outside_transmissions
                            )
  like_df$like2 <- like_df$like / sum(like_df$like)

  if(multiple_outside_transmissions){
    sub_data <- sub_data %>%
      dplyr::filter(.data$gen != 1)
  }

  ## Join likelihood to sub_data
  joined_data <- dplyr::left_join(sub_data,
                                  like_df,
                                  by = c("cluster_id",
                                         "orig_id"))

  ## Compute adjacency matrix
  A_infects_B_mat <- compute_A_infects_B(joined_data)

  if(output == "matrix"){
    return(A_infects_B_mat)
  }

  ## Convert to long data frame

  nms <- round(as.numeric(colnames(A_infects_B_mat)), 4)
  rownames(A_infects_B_mat) <- NULL
  colnames(A_infects_B_mat) <- paste0("infected_id-", nms)

  A_infects_B_wide <- as.data.frame(A_infects_B_mat)
  A_infects_B_wide$infector_id <- nms
  A_infects_B <- A_infects_B_wide %>%
    tidyr::pivot_longer(-.data$infector_id,
                        names_prefix = "infected_id-",
                        names_to = "infected_id",
                        values_to = "prob") %>%
      dplyr::mutate(infected_id =
                        as.numeric(.data$infected_id))




  return(A_infects_B)
}


#' Get the probability of transmission given a set of covariates
#'
#' @param data data frame with variables corresponding to those in cov_names
#' @param inf_params values of parameters, possibly named
#' @param covariate_names covariate names
#' @return probability of infection, the result of a logit function
calculate_transmission_prob <- function(data,
                            inf_params,
                            covariate_names = names(inf_params)[-1]){
  if(length(covariate_names) == 0){
    covariate_names <- NA
  }

  cov_mat <- covariate_df_to_mat(data, covariate_names)
  prob_trans <- as.numeric(1 / (1 + exp(- cov_mat %*% inf_params)))
  return(prob_trans)


}


#' Compute probability from samples that person A infects person B
#'
#' @param df data frame with the following columns
#' \describe{
#' \item{cluster_id}{unique ID of cluster}
#' \item{prob_trans}{individual probability of transmission}
#' \item{like}{likelihood of cluster}
#' \item{id}{id of person within cluster}
#' \item{inf_id}{infector ID of person}
#' \item{feature_id}{a unique person/or set of features across the MC clusters}
#' }
#' @return a n x n matrix where n is the number of individuals in the cluster.
#' The rows correspond to the INFECTOR id and the columns to the INFECTED id.
#' So entry i,j is the probability i infects j
compute_A_infects_B <- function(df){

  ## This is the AVERAGE likelihood in sample
  ## Initialize matrix
  nms <- factor(sort(unique(df$feature_id)))
  nms2 <- sort(unique(df$feature_id))
  mat <- matrix(0, ncol = length(nms),
                nrow = length(nms))
  rownames(mat) <- nms2
  colnames(mat) <- nms2
  ## Loop through to see who infected whom
  for(ii in 1:nrow(df)){
    my_cluster_id <- df$cluster_id[ii]
    sub_df <- df %>%
      dplyr::filter(.data$cluster_id == my_cluster_id)
    index_id <- as.numeric(factor(df$feature_id[ii],
                       levels = levels(nms)))
    infector_index <- which(sub_df$id == df$inf_id[ii])
    infector_id <- as.numeric(factor(df$feature_id[infector_index],
                          levels = levels(nms)))
    like <- sub_df$like2[1]

    if(length(infector_index) > 0){
      infector_prob <- sub_df$prob_trans[infector_index]

      mat[infector_id, index_id] <-  mat[infector_id, index_id] +
        like * infector_prob
    }




  }
  return(mat)


}
