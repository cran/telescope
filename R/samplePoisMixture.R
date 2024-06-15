#' Telescoping sampling of a Bayesian finite Poisson mixture with a
#' prior on the number of components K.
#'
#' @description
#'   * The MCMC scheme is implemented as suggested in Frühwirth-Schnatter et al (2021).
#'   * The priors on the model parameters are specified as in
#'     Frühwirth-Schnatter et al (2021) and Früwirth-Schnatter and
#'     Malsiner-Walli (2019), see the vignette for details and notation.
#'            
#' @param y A numeric matrix; containing the data.
#' @param S A numeric matrix; containing the initial cluster
#'     assignments.
#' @param mu A numeric matrix; containing the initial cluster-specific
#'     rate values.
#' @param eta A numeric vector; containing the initial cluster sizes.
#' @param a0 A numeric vector; hyperparameter of the prior on the rate \eqn{\mu}.
#' @param b0 A numeric vector; hyperparameter of the prior on the rate \eqn{\mu}.
#' @param h0 A numeric vector; hyperparameter of the prior on the rate \eqn{\mu}.
#' @param H0 A numeric vector; hyperparameter of the prior on the rate \eqn{\mu}.
#' @param M A numeric scalar; specifying the number of recorded
#'     iterations.
#' @param burnin A numeric scalar; specifying the number of burn-in
#'     iterations.
#' @param thin A numeric scalar; specifying the thinning used for the
#'     iterations.
#' @param Kmax A numeric scalar; the maximum number of components. 
#' @param G A character string; either `"MixDynamic"` or `"MixStatic"`.
#' @param priorOnK A named list; providing the prior on the number of components K, see [priorOnK_spec()].
#' @param priorOnWeights A named list; providing the prior on the mixture weights.
#' @param verbose A logical; indicating if some intermediate clustering
#'     results should be printed.
#' @return A named list containing:
#'   * `"Mu"`: sampled rate \eqn{\mu}.
#'   * `"Eta"`: sampled weights.
#'   * `"S"`: sampled assignments.
#'   * `"Nk"`: number of observations assigned to the different components, for each iteration.
#'   * `"K"`: sampled number of components.
#'   * `"Kplus"`: number of filled, i.e., non-empty components, for each iteration.
#'   * `"e0"`: sampled Dirichlet parameter of the prior on the weights (if \eqn{e_0} is random).
#'   * `"alpha"`: sampled Dirichlet parameter of the prior on the weights (if \eqn{\alpha} is random).
#'   * `"acc"`: logical vector indicating acceptance in the Metropolis-Hastings step when sampling either e0 or \eqn{\alpha}.
#'
#' @examples
#' N <- 200
#' z <- sample(1:2, N, prob = c(0.5, 0.5), replace = TRUE)
#' y <- rpois(N, c(1, 6)[z])
#' 
#' Mmax <- 200
#' thin <- 1
#' burnin <- 100
#' M <- Mmax/thin
#' 
#' Kmax <- 50  
#' Kinit <- 10
#' 
#' G <- "MixDynamic"
#' priorOnAlpha <- priorOnAlpha_spec("gam_1_2")
#' priorOnK <- priorOnK_spec("BNB_143")
#' 
#' a0 <- 0.1 
#' h0 <- 0.5 
#' b0 <- a0/mean(y) 
#' H0 <- h0/b0
#' 
#' cl_y <- kmeans(y, centers = Kinit, nstart = 100)
#' S_0 <- cl_y$cluster
#' mu_0 <- t(cl_y$centers)
#' eta_0 <- rep(1/Kinit, Kinit)
#' 
#' result <- samplePoisMixture(
#'   y, S_0, mu_0, eta_0, 
#'   a0, b0, h0, H0,
#'   M, burnin, thin, Kmax, 
#'   G, priorOnK, priorOnAlpha)
#' 
#' K <- result$K
#' Kplus <- result$Kplus
#' 
#' plot(seq_along(K), K, type = "l", ylim = c(0, max(K)), 
#'      xlab = "iteration", main = "",
#'      ylab = expression("K" ~ "/" ~ K["+"]), col = 1)
#' lines(seq_along(Kplus), Kplus, col = 2)
#' legend("topright", legend = c("K", expression(K["+"])),
#'        col = 1:2, lty = 1, box.lwd = 0)
#' 
samplePoisMixture <-
    function(y, S, mu, eta, a0, b0, h0, H0,
             M, burnin, thin, Kmax,
             G = c("MixDynamic", "MixStatic"),
             priorOnK, priorOnWeights,
             verbose = FALSE) {

    y <- as.matrix(y)
    ## initial number of componens
    K_j <- length(eta)

    ## prior on K and weights
    log_pK <- priorOnK$log_pK
    G <- match.arg(G)
    if (G == "MixDynamic") {
        log_pAlpha <- priorOnWeights$log_pAlpha  
        a_alpha <- priorOnWeights$param$a_alpha        
        b_alpha <- priorOnWeights$param$b_alpha        
        alpha <- priorOnWeights$param$alpha
        e0 <- alpha / K_j
    } else {
        e0 <- priorOnWeights$param$e0
        alpha <- e0 * K_j
        log_p_e0 <- priorOnWeights$log_p_e0
    }
    s0_proposal <- priorOnWeights$param$s0_proposal

    N <- nrow(y)  # number of observations
    r <- ncol(y)  # number of dimensions

  ## initializing current values
  eta_j <- eta
  mu_j <- mu
  b0_j <- b0
  S_j <- S

  Nk_j <- tabulate(S_j, K_j)
  Kp_j <- sum(Nk_j != 0)  # number of nonempty components
  acc <- FALSE

  ## generating matrices for storing the draws:
  result <- list(Eta = matrix(NA_real_, M, Kmax),
                 Mu = array(NA_real_, dim = c(M, r, Kmax)),
                 b0 = array(NA_real_, dim = c(M, r)),
                 Nk = matrix(NA_integer_, M, Kmax),
                 S = matrix(NA_integer_, M, N),
                 K = rep(NA_integer_, M),
                 Kp = rep(0L, M),
                 mixlik = rep(0, M),
                 mixprior = rep(0, M),
                 nonnormpost = rep(0, M),
                 nonnormpost_mode_list = vector("list", Kmax),
                 mixlik_mode_list = vector("list", Kmax),
                 e0 = rep(NA_real_, M),
                 alpha = rep(NA_real_, M),
                 acc = rep(NA, M)
                 )

  ## Initialising the result-objects
  result$Mu[1, , 1:K_j] <- mu_j
  result$Eta[1, 1:K_j] <- eta_j
  result$b0[1, 1] <- b0
  result$S[1, ] <- S_j
  result$K[1] <- K_j
  result$Kplus[1] <- Kp_j
  result$Nk[1, 1:K_j] <- Nk_j
  result$e0[1] <- e0
  result$alpha[1] <- alpha
  result$acc[1] <- FALSE
  for (k in 1:Kmax) {
    result$nonnormpost_mode_list[[k]] <- list(nonnormpost = -(10)^18)
    result$mixlik_mode_list[[k]] <- list(mixlik = -(10)^18)
  }

  ##---------------------- simulation ----------------------------------------------

  s <- 1
  m <- 2
  Mmax <- M * thin
  while (m <= Mmax || m <= burnin) {
    if (verbose && !(m%%1000)) {
        cat("\n", m, " ", Nk_j)  
    }

    if (m == burnin) {
        m <- 1
        burnin <- 0
    }

    ## first step: classify observations and determine new partition
    mat <- sapply(1:K_j, function(k) eta_j[k] * ((mu_j[k])^y[, 1]) * exp(-mu_j[k]))
    S_j <- apply(mat, 1, function(x) sample(1:K_j, 1, prob = x, replace = T))

    ## determine P
    Nk_j <- tabulate(S_j, K_j)  #length(Nk_j)=K_j
    Kp_j <- sum(Nk_j != 0)

    ## reorder the components
    perm <- c(which(Nk_j > 0), which(Nk_j == 0))
    mu_j <- mu_j[, perm, drop = FALSE]  #length(mu_j[1,])=K_j
    S_ <- rep(FALSE, N)
    for (i in 1:length(perm)) {
        S_[S_j == i] <- which(perm == i)
    }
    S_j <- S_
    Nk_j <- tabulate(S_j, Kp_j)  #length(Nk_j)=Kp_j

    ## second step: parameter simulation conditional on partition P=(N_1,...,N_K+):

    ## (2a) update parameters of filled components

    mu_j[1, 1:Kp_j] <- sapply(1:Kp_j, function(k) {
        rgamma(1, shape = a0 + Nk_j[k] * mean(y[S_j == k, ]), rate = b0 + Nk_j[k])
    })
    ## storing the moments for clustering the draws in the point process representation
    mean_yk <- matrix(sapply(1:Kp_j, function(k) colMeans(y[S_j == k, , drop = FALSE])), ncol = Kp_j)  
    ak <- as.vector(a0) + as.vector(Nk_j * mean_yk)
    bk <- as.vector(b0) + as.vector(Nk_j)

    ## (2b) sample hyperparameters conditional on P

    ## (i): sample b0
    b0 <- rgamma(1, shape = h0 + Kp_j * a0, rate = H0 + sum(mu_j[1, 1:Kp_j]))  #ATTENTION!Only for learnK

    ## third step: sample K and alpha (e0) conditional on partition

    if (G == "MixDynamic") {
        ## If e0=alpha/K (=dependent on K) (3a) Sample K
        K_j <- sampleK_alpha(Kp_j, Kmax, Nk_j, alpha, log_pK)
        
        ## (3b) Sample alpha, if alpha~p(a_alpha,b_alpha)
        value <- sampleAlpha(N, Nk_j, K_j, alpha, s0_proposal, log_pAlpha)
        alpha <- value$alpha
        e0 <- alpha / K_j
        acc <- value$acc
    } else {
        ## If e0 fixed or e0~G(a_e,b_e) (independent of K): (3a*) Sample K
        K_j <- sampleK_e0(Kp_j, Kmax, log_pK, log_p_e0, e0, N)

        ## (3b*) Sample e0, if e0~G(a_e,b_e) (independent of K)
        value <- sampleE0(K_j, Kp_j, N, Nk_j, s0_proposal, e0, log_p_e0)
        e0 <- value$e0
        alpha <- e0 * K_j
        acc <- value$acc
    }
    
    ## fourth step: add empty components conditional on K

    ## (4a) Add/remove empty components
    if (K_j > Kp_j) {
        Nk_j <- c(Nk_j[1:Kp_j], rep(0, (K_j - Kp_j)))

        mu_j <- cbind(mu_j[, 1:Kp_j, drop = FALSE], matrix(0, r, K_j - Kp_j))
        mu_j[, (Kp_j + 1):K_j] <- rgamma(K_j - Kp_j, shape = a0, rate = b0)
    } else {
        mu_j <- mu_j[, 1:K_j, drop = FALSE]
    }

    ## (4b): Sample eta_j:
    ek <- e0 + Nk_j
    eta_j <- MCMCpack::rdirichlet(1, ek)

    ## fifth step: evaluating the mixture likelihood and storing the values

    ## evaluating the mixture likelihood
    mat_neu <- sapply(1:K_j, function(k) eta_j[k] * dpois(y, lambda = mu_j[, k], log = FALSE))
    mixlik_j <- sum(log(rowSums(mat_neu)))

    ## evaluating the mixture prior
    mixprior_j <- log(MCMCpack::ddirichlet(as.vector(eta_j), rep(e0, K_j))) + sum(dgamma(mu_j, shape = a0,
      rate = b0, log = TRUE)) + dgamma(b0, shape = h0, rate = H0) + log_pK(K_j)

    if (burnin == 0) {
      result$mixlik[m] <- mixlik_j
      result$mixprior[m] <- mixprior_j
      result$nonnormpost[m] <- result$mixlik[m] + result$mixprior[m]
    }

    ## storing the nonnormalized posterior for having good starting
    ## points for clustering the draws in the point process repres.
    if ((burnin == 0) && (result$nonnormpost[m] > result$nonnormpost_mode_list[[Kp_j]]$nonnormpost)) {
        result$nonnormpost_mode_list[[Kp_j]] <- list(nonnormpost = result$nonnormpost[m],
                                                     mu = mu_j[, Nk_j != 0],
                                                     mean_muk = ak/bk,
                                                     var_muk = ak/(bk^2),
                                                     eta = eta_j)
    }

    ## storing the results
    if ((burnin == 0) && !(m%%thin)) {
        ## storing the new values
        result$Mu[m/thin, , 1:K_j] <- mu_j
        result$Eta[m/thin, 1:K_j] <- eta_j
        result$b0[m/thin, ] <- b0
        result$S[m/thin, ] <- S_j
        result$Nk[m/thin, 1:K_j] <- Nk_j
        result$K[m/thin] <- K_j
        result$Kplus[m/thin] <- Kp_j
        result$e0[m/thin] <- e0
        result$alpha[m/thin] <- alpha
        result$acc[m/thin] <- acc
    }
    m <- m + 1
  }

  return(result)
}


