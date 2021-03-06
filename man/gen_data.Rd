% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation-lib.R
\name{gen_data}
\alias{gen_data}
\title{Generate both the covariates and the outcomes}
\usage{
gen_data(n, true_betas, gen_dist = gen_xs_default, error_dist = rnorm, ...)
}
\arguments{
\item{n}{Number of observations}

\item{true_betas}{A vector of true values of regression coefficients.
Should include the intercept as the first item, beta1 as the second,
etc.}

\item{gen_dist}{A function specifiying how to generate the
covariates. Should take \code{n} as its first argument, corresponding
to the number of observations, and \code{p} as the second,
corresponding to the number of dimensions per observation (including
the intercept, so each X_i has p-1 dimensions).}

\item{error_dist}{A function specifying how to generate the errors.
Should take \code{n} as the first argument, specifying the number of
observations}

\item{...}{Additional arguments to gen_dist}
}
\value{
A list \code{lst}, with element \code{lst$xs} corresponding to the
covariates and \code{lst$ys} corresponding to the responses.
}
\description{
Generates X using \code{gen_dist} and then generates
Y = beta_0 + XB + eps, #' where eps is generated according to
\code{error_dist}. NOT exported to the user.
}

