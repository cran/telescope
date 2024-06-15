#' Specify prior on \eqn{K}.
#'
#' @description Obtain a function to evaluate the log prior
#'     specified for \eqn{K}.
#' @param P A character indicating which specification should be
#'     used. See Details for suitable values.
#' @param K A numeric or integer scalar specifying the fixed (if `P`
#'     equals `"fixedK"`) or maximum value (if `P` equals `"Unif"`) of
#'     \eqn{K}.
#' @return A named list containing:
#' * `"log_pK"`: a function of the log prior of \eqn{K}.
#' * `"param"`: a list with the parameters.
#' 
#' @details
#' The following prior specifications are supported:
#' * `"fixedK"`: K has the fixed value K (second argument). 
#' * `"Unif"`: \eqn{K \sim} Unif\eqn{[1,K]}, where the upper limit is given by K (second argument). 
#' * `"BNB_111"`: \eqn{K-1 \sim} BNB(1,1,1), i.e., \eqn{K-1} follows a beta-negative binomial distribution with parameters \eqn{(1,1,1)}. 
#' * `"BNB_121"`: \eqn{K-1 \sim} BNB(1,2,1), i.e., \eqn{K-1} follows a beta-negative binomial distribution with parameters \eqn{(1,2,1)}.
#' * `"BNB_143"`: \eqn{K-1 \sim} BNB(1,2,1), i.e., \eqn{K-1} follows a beta-negative binomial distribution with parameters \eqn{(1,4,3)}. 
#' * `"BNB_443"`: \eqn{K-1 \sim} BNB(4,4,3), i.e., \eqn{K-1} follows a beta-negative binomial distribution with parameters \eqn{(4,4,3)}.
#' * `"BNB_943"`: \eqn{K-1 \sim} BNB(9,4,3), i.e., \eqn{K-1} follows a beta-negative binomial distribution with parameters \eqn{(9,4,3)}.
#' * `"Pois_1"`: \eqn{K-1 \sim} pois(1), i.e., \eqn{K-1} follows a Poisson distribution with rate 1.
#' * `"Pois_4"`: \eqn{K-1 \sim} pois(4), i.e., \eqn{K-1} follows a Poisson distribution with rate 4.
#' * `"Pois_9"`: \eqn{K-1 \sim} pois(9), i.e., \eqn{K-1} follows a Poisson distribution with rate 9.
#' * `"Geom_05"`: \eqn{K-1 \sim} geom(0.5), i.e., \eqn{K-1} follows a geometric distribution with success probability \eqn{p=0.5} and density \eqn{f(x)=p(1-p)^x}.
#' * `"Geom_02"`: \eqn{K-1 \sim} geom(0.2), i.e., \eqn{K-1} follows a geometric distribution with success probability \eqn{p=0.2} and density \eqn{f(x)=p(1-p)^x}.
#' * `"Geom_01"`: \eqn{K-1 \sim} geom(0.1), i.e., \eqn{K-1} follows a geometric distribution with success probability \eqn{p=0.1} and density \eqn{f(x)=p(1-p)^x}.
#' * `"NB_11"`: \eqn{K-1 \sim} nbinom(1,0.5), i.e., \eqn{K-1} follows a negative-binomial distribution with \eqn{size=1} and \eqn{p=0.5}.
#' * `"NB_41"`: \eqn{K-1 \sim} nbinom(4,0.5), i.e., \eqn{K-1} follows a negative-binomial distribution with \eqn{size=4} and \eqn{p=0.5}.
#' * `"NB_91"`: \eqn{K-1 \sim} nbinom(9,0.5), i.e., \eqn{K-1} follows a negative-binomial distribution with \eqn{size=9} and \eqn{p=0.5}.
#' 
priorOnK_spec <- function(P = c("fixedK", "Unif",
                                "BNB_111", "BNB_121", "BNB_143", "BNB_443", "BNB_943",
                                "Pois_1", "Pois_4", "Pois_9",
                                "Geom_05", "Geom_02", "Geom_01",
                                "NB_11", "NB_41", "NB_91"), K) {
    P <- match.arg(P)
    if (P %in% c("fixedK", "Unif")) {
        stopifnot(is.numeric(K), length(K) == 1, K >= 1)
        K <- as.integer(K)
    }
    
    if (P == "fixedK") {
        param <- list(K_0 = K)
        log_pK <-function (x) log(x == K)
    }
    if (P == "Unif") {
        param <- list(Kmax = K)
        log_pK <-function (x) log(1/K)
    }
    if (P == "BNB_111") {
        alpha.B <- 1   
        a_pi <- 1
        b_pi <- 1
        param <- list(alpha.B = alpha.B, a_pi = a_pi, b_pi = b_pi)
        log_pK <- function (x)
            dbnbinom(x, size = alpha.B, alpha = a_pi, beta = b_pi, log = TRUE) - log(0.5)
    }
    if (P == "BNB_212") {
        alpha.B <- 2    
        a_pi <- 1
        b_pi <- 2
        param <- list(alpha.B = alpha.B, a_pi = a_pi, b_pi = b_pi)
        log_pK <- function (x)
            dbnbinom(x-1, size = alpha.B, alpha = a_pi, beta = b_pi, log = TRUE)
    }
    if (P == "BNB_143") {
        alpha.B <- 1         
        a_pi <- 4
        b_pi <- 3
        param <- list(alpha.B = alpha.B, a_pi = a_pi, b_pi = b_pi)
        log_pK <- function (x)
            dbnbinom(x-1, size = alpha.B, alpha = a_pi, beta = b_pi, log = TRUE)
    }
    if (P == "BNB_443") {
        alpha.B <- 4         
        a_pi <- 4
        b_pi <- 3
        param <- list(alpha.B = alpha.B, a_pi = a_pi, b_pi = b_pi)
        log_pK <- function (x)
            dbnbinom(x-1, size = alpha.B, alpha = a_pi, beta = b_pi, log = TRUE)
    }
    if (P == "BNB_943") {
        alpha.B <- 9         
        a_pi <- 4
        b_pi <- 3
        param <- list(alpha.B = alpha.B, a_pi = a_pi, b_pi = b_pi)
        log_pK <- function (x)
            dbnbinom(x-1, size = alpha.B, alpha = a_pi, beta = b_pi, log = TRUE)
    } 
    if (P == "Pois_1") {
        lambda <- 1          
        param <- list(lambda = lambda)
        log_pK <- function(x) dpois(x-1, lambda, log = TRUE)
    } 
    if (P == "Pois_4") {
        lambda <- 4          
        param <- list(lambda = lambda)
        log_pK <- function(x) dpois(x-1, lambda, log = TRUE)
    } 
    if (P == "Pois_9") {
        lambda <- 9          
        param <- list(lambda = lambda)
        log_pK <- function(x) dpois(x-1, lambda, log = TRUE)
    }
    if (P == "Geom_05") {
        p_geom <- 0.5   
        param <- list(p_geom = p_geom)
        log_pK <- function(x) dgeom(x-1, p_geom, log = TRUE)  
    }   
    if (P == "Geom_02") {
        p_geom <- 0.2
        param <- list(p_geom = p_geom)
        log_pK <- function(x) dgeom(x-1, p_geom, log = TRUE)  
    } 
    if (P == "Geom_01") {
        p_geom <- 0.1 
        param <- list(p_geom = p_geom)
        log_pK <- function(x) dgeom(x-1, p_geom, log = TRUE)  
    }  
    if (P == "NB_11") {
        alpha.nb <- 1        
        beta.nb <- 1 
        param <- list(alpha.nb = alpha.nb, beta.nb = beta.nb)
        log_pK <- function (x)
            dnbinom(x-1, size = alpha.nb, prob = beta.nb/(beta.nb+1), log = TRUE)
    }   
    if (P == "NB_41") {
        alpha.nb <- 4        
        beta.nb <- 1 
        param <- list(alpha.nb = alpha.nb, beta.nb = beta.nb)
        log_pK <- function (x)
            dnbinom(x-1, size = alpha.nb, prob = beta.nb/(beta.nb+1), log = TRUE)
    }  
    if (P == "NB_91") {
        alpha.nb <- 9   
        beta.nb <- 1 
        param <- list(alpha.nb = alpha.nb, beta.nb = beta.nb)
        log_pK <- function (x)
            dnbinom(x-1, size = alpha.nb, prob = beta.nb/(beta.nb+1), log = TRUE)
    }  
    return(list(log_pK = log_pK, param = param))
}
