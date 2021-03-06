% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis-A-infects-B.R
\name{mc_trees_to_A_infects_B}
\alias{mc_trees_to_A_infects_B}
\title{Looks at most likely transmissions from A to B for a cluster}
\usage{
mc_trees_to_A_infects_B(
  my_id,
  par_ests,
  mc_trees,
  output = "tidy",
  multiple_outside_transmissions = FALSE
)
}
\arguments{
\item{my_id}{ID of cluster to analyze}

\item{par_ests}{parameter estimates (must be named and
have corresponding columns in \code{mc_trees})}

\item{mc_trees}{output from \code{\link{sample_mc_trees}}}

\item{output}{either \code{"tidy"} or \code{"matrix"}}

\item{multiple_outside_transmissions}{logical indicating which model to calculate likelihood}
}
\value{
either a \eqn{n \times n}  matrix where \eqn{n} is the number of individuals
in the cluster and entry
\eqn{i,j} is the probability that individual A infects individual B
or a data frame with the following columns
\describe{
\item{infector}{ID of the infector}
}
}
\description{
Looks at most likely transmissions from A to B for a cluster
}
