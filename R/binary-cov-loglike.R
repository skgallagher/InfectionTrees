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
#' \item{n}{number in cluster}
#' \item{n_pos}{number of positive smears in cluster}
#' \item{freq}{number of times this frequency occurs}
#' }
#' @param mc_samples_summary a data frame of sampled transmission trees
#' with the following columns,
#' \describe{
#' \item{x_pos_trans}{number of transmission from a smear positive person}
#' \item{x_pos}{number of positive smears}
#' \item{cluster_size}{cluster size}
#' \item{mc_freq}{number of samples from}
#' }
#' @param return_neg logical indicating whether to return the negative loglike.  Default is FALSE
#' @param messages logical indicating whether we should print messages.  Default is FALSE.
#' @return the estimated log average likelihood over the observed clusters
#' @export
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


#    browser()

    x_pos <- x_neg <-  NULL
    mc_freq <- freq <- NULL
  df_join <- df_join[, prob_inf := my_prob_inf]
  like_df <- df_join[,
                                .(mc_avg_like = sum(mc_freq * prob_inf)/
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
