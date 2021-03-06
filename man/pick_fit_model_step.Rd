% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation-lib.R
\name{pick_fit_model_step}
\alias{pick_fit_model_step}
\title{Pick the best model using \code{step} and return its fitted
confidence intervals.}
\usage{
pick_fit_model_step(xs, ys, direction)
}
\arguments{
\item{xs}{covariates}

\item{ys}{responses}

\item{direction}{A valid argument to the \code{direction} argument of
the \code{step} command}
}
\value{
95% confidence intervals from the selected model.
}
\description{
Pick the best model using \code{step} and return its fitted
confidence intervals. NOT exported to the user.
}

