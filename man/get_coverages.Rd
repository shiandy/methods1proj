% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation-lib.R
\name{get_coverages}
\alias{get_coverages}
\title{Get the coverage of each coefficient}
\usage{
get_coverages(res_df, true_betas)
}
\arguments{
\item{res_df}{data.frame returned by \code{\link{run_sim}}}

\item{true_betas}{A vector of true values of regression coefficients.
Should include the intercept as the first item, beta1 as the second,
etc.}
}
\value{
A vector of coverages, one per each coefficient. The vector
  is nameded, with each element corresponding to the coverage
  probability of one coefficient in one iteration of the simulation
  (is NA if that coefficient was never selected).
}
\description{
Get the coverage of each coefficient
}

