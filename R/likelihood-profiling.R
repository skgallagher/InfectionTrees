## SKG
## Likelihood profile wrappers for loglike functions
## Aug. 18, 2020


#' Wrapper to optimize over a single variable for likelihood profiling
#'
#' @param x parameter to optimize over
#' @param max_loglike maximum loglike value from output
#' @param best_pars vector of best parameters
#' @param mc_samples_summary output of \code{sample_mc_binary_cov}
#' @param alpha \eqn{\alpha}-value for computing confidence intervals.
#' @param beta_index which parameter to optimize over.
#' @param multiple_outside_transmissions whether to use base model or not
#' @return zeros of likelihood profiling for given \eqn{alpha}
#' @details internal function only
binary_loglike_wrapper <- function(x, max_loglike,
                                   best_pars,
                                   mc_samples_summary,
                                   alpha = .05,
                                   beta_index = 1,
                                   multiple_outside_transmissions = FALSE){

  inf_params <- best_pars
  inf_params[beta_index] <- x

  chi_stat <- stats::qchisq((1 - alpha), df = 1) / 2

  loglike <- bp_loglike_binary_cov(inf_params = inf_params,
                                   mc_samples_summary = mc_samples_summary,
                                   return_neg = FALSE,
                                   multiple_outside_transmissions =
                                     multiple_outside_transmissions)

  return(loglike - max_loglike + chi_stat)

}


#' Return likelihood profiles for binary covariates special case
#'
#' @param best_pars vector of MLE
#' @param max_loglike loglike value from MLE
#' @param mc_samples_summary tree summaries, output from \code{sample_mc_binary_cov}
#' @param alpha \eqn{alpha}-value to make a (1-\eqn{alpha})% CI
#' @param multiple_outside_transmissions logical indicating to use base model
#'  or multiple outside transmissions model.
#' @return matrix of dimension (length of best_pars) x 3 where the columns are the best parameter estimate,
#'  the lower CI (\code{lower}) and upper CI (\code{upper})
binary_likelihood_profs <- function(best_pars,
                                    max_loglike,
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
                         mc_samples_summary = mc_samples_summary,
                         alpha = .05,
                         multiple_outside_transmissions = multiple_outside_transmissions)
        ##
        upper <- stats::uniroot(f = binary_loglike_wrapper,
                                c(best_pars[jj], best_pars[jj] + 4),
                         max_loglike = max_loglike,
                         best_pars = best_pars,
                         beta_index = jj,
                         mc_samples_summary = mc_samples_summary,
                         alpha = .05,
                         multiple_outside_transmissions = multiple_outside_transmissions)

        ci_mat[jj, 1] <- lower$root
        ci_mat[jj, 2] <- upper$root
    }
    est_pars <- cbind(best_pars, ci_mat)
    colnames(est_pars) <- c("est.", "lower", "upper")


  return(est_pars)


}
