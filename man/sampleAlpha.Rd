% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_e0_alpha.R
\name{sampleAlpha}
\alias{sampleAlpha}
\title{Sample alpha conditional on partition and K using an
Metropolis-Hastings step with log-normal proposal.}
\usage{
sampleAlpha(N, Nk, K, alpha, s0_proposal, log_pAlpha)
}
\arguments{
\item{N}{A number; indicating the sample size.}

\item{Nk}{An integer vector; indicating the group sizes in the partition.}

\item{K}{A number; indicating the number of components.}

\item{alpha}{A numeric value; indicating the value for \eqn{\alpha}.}

\item{s0_proposal}{A numeric value; indicating the standard deviation of the random walk.}

\item{log_pAlpha}{A function; evaluating the log prior of \eqn{\alpha}.}
}
\value{
A named list containing:
\itemize{
\item \code{"alpha"}: a numeric, the new \eqn{\alpha} value.
\item \code{"acc"}: logical indicating acceptance.
}
}
\description{
Sample \eqn{\alpha} conditional on the current
partition and value of \eqn{K} using an Metropolis-Hastings
step with log-normal proposal.
}
