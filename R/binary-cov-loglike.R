## binary covariate log like
## SKG
## Aug 17, 2020


#' Calculate the likelihood for binary covariate MC trees
#' 
#' @param inf_params vector of length 2, the beta0 and beta1 coefficients
#' @param mc_samples_summary a data frame of sampled transmission trees
#' with the following columns, where the last one is only needed for the multiple
#' outsider model.
#' 
#' @param return_neg logical indicating whether to return the negative loglike.  Default is FALSE
#' @param multiple_outsider_transmissions logical indicating whether to use the outsider model likelihood or not
#' @return the estimated log average likelihood over the observed clusters
bp_loglike_binary_cov <- function(inf_params,
                                  mc_samples_summary,
                                  return_neg = FALSE,
                                  multiple_outsider_transmissions = FALSE){

  total_p_pos <- sum(mc_samples_summary$freq *
                       mc_samples_summary$x_pos)
  total_p_neg <- sum(mc_samples_summary$freq *
                       mc_samples_summary$x_neg)
  emp_p_pos <- total_p_pos / (total_p_neg + total_p_pos)

  p_neg <- 1 / (1 + exp(-inf_params[1]))
    p_pos <- 1 / (1 + exp(-sum(inf_params)))
    prob_inf <- NULL

  if(!("data.table" %in% class(mc_samples_summary))){
    print("Converting 'sampled_data' to data.table format")
    mc_samples_summary <- data.table::as.data.table(mc_samples_summary)
  }
  if(!multiple_outsider_transmissions){
    my_prob_inf <- (1-p_pos)^mc_samples_summary$x_pos *
      (1-p_neg)^mc_samples_summary$x_neg *
      p_pos^mc_samples_summary$x_pos_trans *
      p_neg^mc_samples_summary$x_neg_trans
  } else{ ## must condition on sign of the imputed outsider
    root_p <- ifelse(mc_samples_summary$root_node == 1, p_pos, p_neg)
    my_prob_inf <- (1-p_pos)^(mc_samples_summary$x_pos + 1) *
      (1-p_neg)^mc_samples_summary$x_neg *
      p_pos^mc_samples_summary$x_pos_trans *
      p_neg^mc_samples_summary$x_neg_trans / root_p

  }

    x_pos <- x_neg <-  NULL
    mc_freq <- freq <- NULL
  mc_samples_summary <- mc_samples_summary[, prob_inf := my_prob_inf]
  like_df <- mc_samples_summary[,
                                .(mc_avg_like = sum(mc_freq * prob_inf)/
                                      sum(mc_freq),
                                  new_freq = sum(freq)),
                                by = .(x_pos, x_neg)]

  loglike <- sum(like_df$new_freq * log(like_df$mc_avg_like))
  if(return_neg){
    loglike <- -loglike
  }

  if(is.na(loglike) | is.nan(loglike) |
     is.infinite(loglike)) {
    stop("loglike is NA/NAN/infinite")
  }

  return(loglike)


}
