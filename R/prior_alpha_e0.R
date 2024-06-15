#' Specify prior on \eqn{\alpha}.
#'
#' @description Obtain a function to evaluate the log prior specified
#'     for \eqn{\alpha}.
#'
#' @param H A character indicating which specification should be used.
#' @return A named list containing:
#' * `"log_pAlpha"`: a function of the log prior of \eqn{\alpha}.
#' * `"param"`: a list with the parameters.
#' @details
#' The following prior specifications are supported:
#' * `"alpha_const"`: \eqn{\alpha} is fixed at 1.
#' * `"gam_05_05"`: \eqn{\alpha \sim} gamma(0.5, 0.5), i.e., shape = 0.5, rate = 0.5.
#' * `"gam_1_2"`: \eqn{\alpha \sim} gamma(1, 2), i.e., shape = 1, rate = 2.
#' * `"F_6_3"`: \eqn{\alpha \sim} F(6, 3), i.e., an F-distribution with degrees of freedom equal to 6 and 3.
#'
priorOnAlpha_spec <- function(H = c("alpha_const", "gam_05_05", "gam_1_2", "F_6_3")) {
    H <- match.arg(H)
    if (H == "alpha_const") {
        param <- list(a_alpha = 1,
                      b_alpha = 1,
                      s0_proposal = 1.5)
        param$alpha <- with(param, a_alpha / b_alpha)
        log_pAlpha <- function(x) log(x == param$alpha)
    }
    if (H == "gam_05_05") {
        param <- list(a_alpha = 0.5,
                      b_alpha = 0.5,
                      s0_proposal = 1.5)
        param$alpha <- with(param, a_alpha / b_alpha)
        log_pAlpha <- function(x)
            dgamma(x, shape = param$a_alpha, rate = param$b_alpha, log = TRUE)  
    }
    if (H == "gam_1_2") {
        param <- list(a_alpha = 1,
                      b_alpha = 2,
                      s0_proposal = 1.5)
        param$alpha <- with(param, a_alpha / b_alpha)
        log_pAlpha <- function(x)
            dgamma(x, shape = param$a_alpha, rate = param$b_alpha, log = TRUE)  
    }
    if (H == "F_6_3") {
        param <- list(a_alpha = 6,
                      b_alpha = 3,
                      alpha = 1,
                      s0_proposal = 2.5)
        log_pAlpha <- function(x)
            df(x, df1 = param$a_alpha, df2 = param$b_alpha, log = TRUE)
    }
    return(list(log_pAlpha = log_pAlpha,
                param = param))
}

#' Specify prior on e0.
#'
#' @description Obtain a function to evaluate the log prior specified
#'     for \eqn{e_0}.
#'
#' @param E A character indicating which specification should be used.
#' @param e0 A numeric scalar giving the fixed value of \eqn{e_0}.
#' @return A named list containing:
#' * `"log_p_e0"`: a function of the log prior of \eqn{e_0}.
#' * `"param"`: a list with the parameters.
#' @details
#' The following prior specifications are supported:
#' * `"G_1_20"`: \eqn{e_0 \sim} gamma(1, 20), i.e., shape = 1, rate = 20.
#' * `"e0const"`: \eqn{e_0} is fixed at `e0`.
priorOnE0_spec <- function(E = c("G_1_20", "e0const"), e0) {
  E <- match.arg(E)
  param <- list(b_alpha = 1,
                e0 = e0,
                s0_proposal = 1.5)
  if (E == "G_1_20") {
      log_p_e0 <- function(x)
          dgamma(x, shape = 1, rate = 20, log = TRUE)  
  }
  if (E == "e0const") {
      log_p_e0 <- function(x) log(x == e0)
  }
  return(list(log_p_e0 = log_p_e0, param = param))
}
