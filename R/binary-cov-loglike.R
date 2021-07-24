## binary covariate log like
## SKG
## Aug 17, 2020
## Update June 15, 2021
## This function is missing the data  part and is calculating the log likelihood incorrectly
## because I'm only looking at MC trees.

#' Calculate the likelihood for binary covariate MC trees
#'
#' @param inf_params vector of length 2, the beta0 and beta1 coefficients
#' @param obs_data_summary
#' \describe{
#' \item{cluster_size}{number in cluster}
#' \item{x_pos}{number of positive smears in cluster}
#' \item{x_neg}{number of negative smears in the cluster}
#' \item{freq}{number of times this frequency occurs}
#' }
#' @param mc_samples_summary a data frame of sampled transmission trees
#' with the following columns,
#' \describe{
#' \item{x_pos_trans}{number of total transmissions from any smear positive person}
#' \item{x_pos}{number of positive smears}
#' \item{x_neg}{number of negative smears}
#' \item{cluster_size}{cluster size}
#' \item{mc_freq}{number of samples}
#' }
#' @param return_neg logical indicating whether to return the negative loglike.  Default is FALSE
#' @param messages logical indicating whether we should print messages.  Default is FALSE.
#' @return the log average likelihood over the observed clusters
#' @details WARNING: obs_data_summary and mc_samples_summary must not both have a 'freq' column.
#' @export
#' @examples
#' library(dplyr)
#' inf_params <- c(-1, 1)
#' mc_samples_summary <- data.frame(freq = c(10, 5, 6, 1, 1),
#'                                 mc_freq = c(10, 10, 10, 3, 7),
#'                                 cluster_size = c(1, 1, 2, 2, 2),
#'                                 x_pos = c(1, 0, 2, 1, 1),
#'                                 x_neg = c(0, 1, 0, 1, 1),
#'                                 x_pos_trans = c(0, 0, 1, 1, 0),
#'                                 x_neg_trans = c(0, 0, 0, 0, 1))
#' obs_data_summary <- mc_samples_summary %>%
#'  group_by(cluster_size, x_pos, x_neg) %>%
#'  summarize(freq = sum(freq), .groups = "drop")
#' bp_loglike_binary_cov(inf_params = inf_params,
#' obs_data_summary = obs_data_summary,
#' mc_samples_summary = mc_samples_summary %>%
#'   dplyr::select(-freq),
#' return_neg = FALSE)
#'
#'
#'
bp_loglike_binary_cov <- function(inf_params,
                                  obs_data_summary,
                                  mc_samples_summary,
                                  return_neg = FALSE,
                                  messages = FALSE){

   # browser()

    df_join <- dplyr::left_join(obs_data_summary,
                                mc_samples_summary,
                                by = c("cluster_size", "x_pos", "x_neg"))

    p_neg <- 1 / (1 + exp(-inf_params[1]))
    p_pos <- 1 / (1 + exp(-sum(inf_params)))
    prob_inf <- NULL

  if(!("data.table" %in% class(df_join))){
    if(messages){
      print("Converting 'sampled_data' to data.table format")
    }
    df_join <- data.table::as.data.table(df_join)
  }

    my_prob_inf <- (1-p_pos)^df_join$x_pos *
        (1-p_neg)^df_join$x_neg *
        p_pos^df_join$x_pos_trans *
        p_neg^df_join$x_neg_trans


#   browser()

    x_pos <- x_neg <-  NULL
    mc_freq <- freq <- cluster_size <- NULL
  df_join <- df_join[, prob_inf := my_prob_inf]
  like_df <- df_join[,
                                .(mc_avg_like = sum(mc_freq * prob_inf) /
                                      sum(mc_freq),
                                  new_freq = sum(mc_freq)),
                                by = .(freq, cluster_size, x_pos)]

  loglike <- sum(like_df$freq * log(like_df$mc_avg_like))
  if(return_neg){
    loglike <- -loglike
  }

  if(is.na(loglike) | is.nan(loglike) |
     is.infinite(loglike)) {
    stop("loglike is NA/NAN/infinite")
  }

  return(loglike)


}
