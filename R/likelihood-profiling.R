## SKG
## Likelihood profile wrappers for loglike functions
## Aug. 18, 2020



#' Return likelihood profiles
#'
#' @param best_pars vector of MLE
#' @param max_loglike loglike value from MLE
#' @param mc_trees trees, output from \code{sample_mc_trees)}
#' @param covariate_names names of covariates
#' @param alpha \eqn{alpha}-value to make a (1-\eqn{alpha})% CI
#' @param multiple_outside_transmissions logical indicating to use base model
#'  or multiple outside transmissions model.
#' @return matrix of dimension (length of best_pars) x 3 where the columns are the best parameter estimate,
#'  the lower CI (\code{lower}) and upper CI (\code{upper})
#'  @export
likelihood_profs <- function(best_pars,
                                    max_loglike,
                                    mc_trees,
                                   covariate_names,
                                    alpha = .05,
                                    multiple_outside_transmissions = FALSE){

  ci_mat <- matrix(0, ncol = 2,
                   nrow = length(best_pars))

  for(jj in 1:nrow(ci_mat)){
    max_loglike <- max_loglike
    lower <- stats::uniroot(f = loglike_wrapper,
                            c(best_pars[jj] - 4, best_pars[jj]),
                            max_loglike = max_loglike,
                            best_pars = best_pars,
                            beta_index = jj,
                            mc_trees = mc_trees,
                            covariate_names = covariate_names,
                            alpha = alpha,
                            multiple_outside_transmissions = multiple_outside_transmissions)
    ##
    upper <- stats::uniroot(f = loglike_wrapper,
                            c(best_pars[jj], best_pars[jj] + 4),
                            max_loglike = max_loglike,
                            best_pars = best_pars,
                            beta_index = jj,
                            mc_trees = mc_trees,
                            covariate_names = covariate_names,
                            alpha = alpha,
                            multiple_outside_transmissions = multiple_outside_transmissions)

    ci_mat[jj, 1] <- lower$root
    ci_mat[jj, 2] <- upper$root
  }
  est_pars <- cbind(best_pars, ci_mat)
  colnames(est_pars) <- c("est.", "lower", "upper")


  return(est_pars)


}



#' Wrapper to optimize over a single variable for likelihood profiling
#'
#' @param x parameter to optimize over
#' @param max_loglike maximum loglike value from output
#' @param best_pars vector of best parameters
#' @param mc_trees output of \code{sample_mc_trees}
#' @param covariate_names covariate names to look over
#' @param alpha \eqn{\alpha}-value for computing confidence intervals.
#' @param beta_index which parameter to optimize over.
#' @param multiple_outside_transmissions whether to use base model or not
#' @return zeros of likelihood profiling for given \eqn{alpha}
#' @details internal function only.  This is for the general case
loglike_wrapper <- function(x,
                            max_loglike,
                            best_pars,
                            mc_trees,
                            covariate_names,
                            alpha = .05,
                            beta_index = 1,
                            multiple_outside_transmissions = FALSE){
  inf_params <- best_pars
  inf_params[beta_index] <- x
  loglike <- general_loglike(inf_params,
                  mc_trees = mc_trees,
                  cov_names = covariate_names,
                  return_neg = FALSE,
                  multiple_outside_transmissions = multiple_outside_transmissions)
  chi_stat <- stats::qchisq((1 - alpha), df = 1) / 2
  return(loglike -
    max_loglike + chi_stat)

}



#' Wrapper to optimize over a single variable for likelihood profiling
#'
#' @param x parameter to optimize over
#' @param max_loglike maximum loglike value from output
#' @param best_pars vector of best parameters
#' @param obs_data_summary actual data summary with cluster_size, x_pos, and x_neg and frequency of each
#' @param mc_samples_summary output of \code{sample_mc_binary_cov}
#' @param alpha \eqn{\alpha}-value for computing confidence intervals.
#' @param beta_index which parameter to optimize over.
#' @param multiple_outside_transmissions whether to use base model or not
#' @return zeros of likelihood profiling for given \eqn{alpha}
#' @details internal function only
binary_loglike_wrapper <- function(x, max_loglike,
                                   best_pars,
                                   obs_data_summary,
                                   mc_samples_summary,
                                   alpha = .05,
                                   beta_index = 1,
                                   multiple_outside_transmissions = FALSE){

  inf_params <- best_pars
  inf_params[beta_index] <- x

  chi_stat <- stats::qchisq((1 - alpha), df = 1) / 2

    loglike <- bp_loglike_binary_cov(inf_params = inf_params,
                                     obs_data_summary = obs_data_summary,
                                   mc_samples_summary = mc_samples_summary,
                                   return_neg = FALSE)

  return(loglike - max_loglike + chi_stat)

}


#' Return likelihood profiles for binary covariates special case
#'
#' @param best_pars vector of MLE
#' @param max_loglike loglike value from MLE
#' @param obs_data_summary actual data summary with cluster_size, x_pos, and x_neg and frequency of each
#' @param mc_samples_summary tree summaries, output from \code{sample_mc_binary_cov}
#' @param alpha \eqn{alpha}-value to make a (1-\eqn{alpha})% CI
#' @param multiple_outside_transmissions logical indicating to use base model
#'  or multiple outside transmissions model.
#' @return matrix of dimension (length of best_pars) x 3 where the columns are the best parameter estimate,
#'  the lower CI (\code{lower}) and upper CI (\code{upper})
#'  @export
binary_likelihood_profs <- function(best_pars,
                                    max_loglike,
                                    obs_data_summary,
                                    mc_samples_summary,
                                    alpha = .05,
                                    multiple_outside_transmissions = FALSE){

  ci_mat <- matrix(0, ncol = 2,
                   nrow = length(best_pars))

   for(jj in 1:nrow(ci_mat)){
        max_loglike <- max_loglike
        lower <- stats::uniroot(f = binary_loglike_wrapper,
                         c(best_pars[jj] - 4, best_pars[jj]),
                         max_loglike = max_loglike,
                         best_pars = best_pars,
                         beta_index = jj,
                         obs_data_summary = obs_data_summary,
                         mc_samples_summary = mc_samples_summary,
                         alpha = alpha)
        ##
        upper <- stats::uniroot(f = binary_loglike_wrapper,
                                c(best_pars[jj], best_pars[jj] + 4),
                         max_loglike = max_loglike,
                         best_pars = best_pars,
                         beta_index = jj,
                         obs_data_summary = obs_data_summary,
                         mc_samples_summary = mc_samples_summary,
                         alpha = alpha)

        ci_mat[jj, 1] <- lower$root
        ci_mat[jj, 2] <- upper$root
    }
    est_pars <- cbind(best_pars, ci_mat)
    colnames(est_pars) <- c("est.", "lower", "upper")


  return(est_pars)


}
