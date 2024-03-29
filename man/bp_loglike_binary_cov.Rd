% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binary-cov-loglike.R
\name{bp_loglike_binary_cov}
\alias{bp_loglike_binary_cov}
\title{Calculate the likelihood for binary covariate MC trees}
\usage{
bp_loglike_binary_cov(
  inf_params,
  obs_data_summary,
  mc_samples_summary,
  return_neg = FALSE,
  messages = FALSE
)
}
\arguments{
\item{inf_params}{vector of length 2, the beta0 and beta1 coefficients}

\item{obs_data_summary}{\describe{
\item{cluster_size}{number in cluster}
\item{x_pos}{number of positive smears in cluster}
\item{x_neg}{number of negative smears in the cluster}
\item{freq}{number of times this frequency occurs}
}}

\item{mc_samples_summary}{a data frame of sampled transmission trees
with the following columns,
\describe{
\item{x_pos_trans}{number of total transmissions from any smear positive person}
\item{x_pos}{number of positive smears}
\item{x_neg}{number of negative smears}
\item{cluster_size}{cluster size}
\item{mc_freq}{number of samples}
}}

\item{return_neg}{logical indicating whether to return the negative loglike.  Default is FALSE}

\item{messages}{logical indicating whether we should print messages.  Default is FALSE.}
}
\value{
the log average likelihood over the observed clusters
}
\description{
Calculate the likelihood for binary covariate MC trees
}
\details{
WARNING: obs_data_summary and mc_samples_summary must not both have a 'freq' column.
}
\examples{
library(dplyr)
inf_params <- c(-1, 1)
mc_samples_summary <- data.frame(freq = c(10, 5, 6, 1, 1),
                                mc_freq = c(10, 10, 10, 3, 7),
                                cluster_size = c(1, 1, 2, 2, 2),
                                x_pos = c(1, 0, 2, 1, 1),
                                x_neg = c(0, 1, 0, 1, 1),
                                x_pos_trans = c(0, 0, 1, 1, 0),
                                x_neg_trans = c(0, 0, 0, 0, 1))
obs_data_summary <- mc_samples_summary \%>\%
 group_by(cluster_size, x_pos, x_neg) \%>\%
 summarize(freq = sum(freq), .groups = "drop")
bp_loglike_binary_cov(inf_params = inf_params,
obs_data_summary = obs_data_summary,
mc_samples_summary = mc_samples_summary \%>\%
  dplyr::select(-freq),
return_neg = FALSE)



}
