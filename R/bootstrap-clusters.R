## SKG
## BOOTSTRAP
## August 26, 2020


#' Resample clusters from an observed data set
#' 
#' @param n_boot number of clusters to resample.  If NULL, will re-sample as many in original cluster
#' @param clusters data frame of observed data with the following columns
#' \describe{
#' \item{cluster_id}{cluster ID}
#' \item{other covariates}{any other covariates}
#' }
#' @return data frame with new cluster_id variable but can be traced back to original cluster
#' through the column \code{original_id}
#' @export
#' @examples 
#'   data <- data.frame(cluster_id = c("a", "a", "a",
#' "b", "b",
#' "c", "c",
#' "d"),
#' x = c(1, 0, 1,
#'      1, 1,
#'      0, 1,
#'      0))
#' boot_clusters <- bootstrap_clusters(clusters = data)
#' boot_clusters
bootstrap_clusters <- function(clusters,
                               n_boot = NULL){
  if(is.null(n_boot)){
    n_boot <- length(unique(clusters$cluster_id))
  }
  
  boot_list <- vector(mode = "list", length = n_boot)
  cluster_ids <- unique(clusters$cluster_id)
  cluster_samples <- sample(cluster_ids, size = n_boot,
                            replace = TRUE)
  boot_list <- lapply(1:n_boot, function(ii){
      id <- cluster_samples[ii]
      clusters %>% 
        dplyr::filter(.data$cluster_id == id) %>%
        dplyr::mutate(original_id = id,
                      cluster_id = ii)
    
  })
  boot_clusters <- dplyr::bind_rows(boot_list) %>%
    dplyr::select(.data$cluster_id, dplyr::everything())
  return(boot_clusters)
}