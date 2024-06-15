#' Sample K conditional on e0 (fixed or random, but not depending on K).
#'
#' @description This sampling step only relies on the current
#'     partition and is independent of the current component-specific
#'     parameters, see Frühwirth-Schnatter et al (2021).
#'
#' @param Kp_j A number; indicating the current value of \eqn{K_+}.
#' @param Kmax A number; indicating the maximum value of \eqn{K}, for which the conditional posterior is evaluated.
#' @param log_pK A function; evaluating the prior of \eqn{K}.
#' @param log_p_e0 A function; evaluating the log prior of \eqn{e_0}.
#' @param e0 A number; indicating the value of \eqn{e_0}.
#' @param N A number; indicating the number of observations.
#' @return A number indicating the new value of \eqn{K}.

sampleK_e0 <- function(Kp_j, Kmax, log_pK, log_p_e0, e0, N) {
    lquot <- function(x)
        lgamma(x + 1) - lgamma(x - Kp_j + 1) + lgamma(e0*x) - lgamma(e0*x + N)
    log_p_e0 <- function(x) 0  

    ## Calculate vector  w with w[K]=p(K|Kp_j,e0,N):
    log_w <- rep(-Inf, Kmax)
    log_w[Kp_j:Kmax] <- lquot(Kp_j:Kmax) + log_pK(Kp_j:Kmax) + log_p_e0(Kp_j:Kmax)  
    max_w <- max(log_w[Kp_j:Kmax])
    log_w[Kp_j:Kmax] <- log_w[Kp_j:Kmax]-max_w
    w <- exp(log_w)
    
    K_j <- sample(1:Kmax, 1, prob = w, replace = TRUE)
    return(K_j)
}

#' Sample K conditional on \eqn{\alpha} where \eqn{e0 = \alpha/K}.
#'
#' @description This sampling step only relies on the current
#'     partition and is independent of the current component-specific
#'     parameters, see Frühwirth-Schnatter et al (2021).
#'
#' @param Kp_j A number; indicating the current value of \eqn{K_+}.
#' @param Kmax A number; indicating the maximum value of \eqn{K} for which the conditional posterior is evaluated.
#' @param Nk_j A numeric vector; indicating the group sizes in the partition, i.e.,
#'   the current number of observations in the filled components. 
#' @param alpha A number; indicating the value of the parameter \eqn{\alpha}.
#' @param log_pK A function; evaluating the log prior of \eqn{K}.
#' @return A number indicating the new value of \eqn{K}.

sampleK_alpha <- function(Kp_j, Kmax, Nk_j, alpha, log_pK) {
    lquot <- function(x)
        lgamma(x+1) - lgamma(x-Kp_j+1) + sum(lgamma(Nk_j+(alpha/x))) -
            Kp_j * (log(x/alpha) + lgamma(1+alpha/x))  
    
    ## Calculate vector w with w[K]=p(K|Kp_j,e0,N):
    log_w <- rep(-Inf, Kmax)
    log_w[Kp_j:Kmax] <- sapply(Kp_j:Kmax, FUN = lquot) + log_pK(Kp_j:Kmax)
    max_w <- max(log_w[Kp_j:Kmax])
    log_w[Kp_j:Kmax] <- log_w[Kp_j:Kmax] - max_w
    w <- exp(log_w)
  
    K_j <- sample(1:Kmax,1,prob=w,replace=TRUE)
    return(K_j)
}
