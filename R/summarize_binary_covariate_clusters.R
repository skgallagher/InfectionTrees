## SKG
## August 17, 2020
## Summarize binary covariate clusters


#' Summarize binary covariate clusters
#'
#' @param df data frame with the following columns
#' \describe{
#' \item{cluster_id}{unique cluster ID}
#' \item{{covariate_name}}{actual name of feature to summarize over, a binary (0/1) covariate}
#' }
#' @param covariate_name name of the single binary covariate
#' @return a data frame with the following columns
#' \describe{
#' \item{freq}{frequency of the following clusters}
#' \item{cluster_size}{total size of the cluster}
#' \item{x_pos}{number of individuals in the cluster with the feature of interest =1}
#' \item{x_neg}{number of individuals in the cluster with the feature of interest = 0}
#' }
#' @details Condense data from data frames about each individuals to the
#' summary of the number of indivduals who have a particular covariate feature (1/0).
#' This assumes the trees are in order by generation.
#' @export
#' @examples
#' example_cluster <- data.frame(cluster_id = c(1, 1, 1,
#' 2, 2,
#' 3, 3, 3, 3,
#' 4),
#' x = c(0, 1, 1,
#' 0, 0,
#' 1, 0, 1, 1,
#' 0))
#' summarize_binary_clusters(example_cluster)
summarize_binary_clusters <- function(df,
                                      covariate_name = "x"){

  stopifnot(covariate_name %in% colnames(df))
  sub_df <- df[, c("cluster_id", covariate_name)]
  colnames(sub_df) <- c("cluster_id", "x")
  stopifnot(any(sub_df$x %in% c(0,1)))


  out_df <- sub_df %>%
    dplyr::group_by(.data$cluster_id) %>%
    dplyr::summarize(cluster_size = dplyr::n(),
                     x_pos = sum(.data$x == 1),
                     x_neg = sum(.data$x == 0)) %>%
    dplyr::ungroup() %>%
      dplyr::group_by(.data$cluster_size, .data$x_pos, .data$x_neg) %>%
    dplyr::summarize(freq = dplyr::n()) %>%
    dplyr::ungroup() %>%
      dplyr::select(.data$freq, .data$cluster_size,
                    .data$x_pos, .data$x_neg)


  ## Add a binary covariate class to help us out in sampling mc and likelihood
  class(out_df) <- c("binary_cov", class(out_df))


  return(out_df)


}
