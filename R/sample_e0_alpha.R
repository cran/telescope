#' Sample alpha conditional on partition and K using an
#' Metropolis-Hastings step with log-normal proposal.
#'
#' @description Sample \eqn{\alpha} conditional on the current
#'     partition and value of \eqn{K} using an Metropolis-Hastings
#'     step with log-normal proposal.

#' @param N A number; indicating the sample size.
#' @param Nk An integer vector; indicating the group sizes in the partition. 
#' @param K A number; indicating the number of components.
#' @param alpha A numeric value; indicating the value for \eqn{\alpha}.
#' @param s0_proposal A numeric value; indicating the standard deviation of the random walk.
#' @param log_pAlpha A function; evaluating the log prior of \eqn{\alpha}.
#' @return A named list containing:
#' * `"alpha"`: a numeric, the new \eqn{\alpha} value.
#' * `"acc"`: logical indicating acceptance.

sampleAlpha <- function(N, Nk, K, alpha, s0_proposal, log_pAlpha) {
  log_post_alpha <- function(x)
      lgamma(x) - lgamma(N+x) + sum(lgamma(Nk+x/K) - lgamma(x/K)) +
          log_pAlpha(x)
  lalpha_p <- log(alpha) + rnorm(1, 0, s0_proposal)
  alpha_p <- exp(lalpha_p)
  lalpha1 <- log_post_alpha(alpha_p) - log_post_alpha(alpha) +
      log(alpha_p) - log(alpha)
  alpha1 <- exp(lalpha1)
  acc <- FALSE
  if (runif(1) <= alpha1) {
    alpha <- alpha_p
    acc <- TRUE
  }
  return(list(alpha = alpha, acc = acc))
}

#' Sample e0 conditional on partition and K using an
#' Metropolis-Hastings step with log-normal proposal.
#'
#' @description Sample \eqn{e_0} conditional on the current partition
#'     and value of \eqn{K} using an Metropolis-Hastings step with
#'     log-normal proposal.
#'
#' @param K A number; indicating the number of components.
#' @param Kp A number; indicating the number of filled components \eqn{K_+}.
#' @param N A number; indicating the sample size.
#' @param Nk An integer vector; indicating the group sizes in the partition. 
#' @param s0_proposal A numeric value; indicating the standard deviation of the random walk proposal.
#' @param e0 A numeric value; indicating the current value of \eqn{e_0}.
#' @param log_p_e0 A function; evaluating the log prior of \eqn{e_0}.
#' @return A named list containing:
#' * `"e0"`: a numeric, the new \eqn{e_0} value.
#' * `"acc"`: logical indicating acceptance.

sampleE0 <- function(K, Kp, N, Nk, s0_proposal, e0, log_p_e0) {
  log_post_e0 <- function(x)
      lgamma(K+1) - lgamma(K+1-Kp) + lgamma(K*x) -
          lgamma(N+K*x) + sum(lgamma(Nk+x) - lgamma(x)) +
          log_p_e0(x)
  le0_p <- log(e0) + rnorm(1, 0, s0_proposal)
  e0_p <- exp(le0_p)
  lalpha1 <- log_post_e0(e0_p) - log_post_e0(e0) + log(e0_p) - log(e0)
  alpha1 <- min(exp(lalpha1), 1)
  acc <- FALSE
  if (runif(1) <= alpha1) {
    e0 <- e0_p
    acc <- TRUE
  }
  return(list(e0 = e0, acc = acc))
}
