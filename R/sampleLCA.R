#' Telescoping sampling of the LCA model where a prior on the number
#' of components K is specified.
#' 
#' @description
#'   * The MCMC scheme is implemented as suggested in Frühwirth-Schnatter et al (2021).
#'   * The priors on the model parameters are specified as in Frühwirth-Schnatter et al (2021),
#'     see the vignette for details and notation.
#'
#' @param y A numeric matrix; containing the data.
#' @param S A numeric matrix; containing the initial cluster
#'     assignments.
#' @param pi A numeric vector; containing the initial cluster-specific
#'     success probabilities.
#' @param eta A numeric vector; containing the initial cluster sizes.
#' @param a0 A numeric vector; containing the parameters of the prior on the
#'     cluster-specific success probabilities.
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
#'  * `"Pi"`: sampled component-specific success probabilities.
#'  * `"Eta"`: sampled weights.
#'  * `"S"`: sampled assignments.
#'  * `"Nk"`: number of observations assigned to the different components, for each iteration.
#'  * `"K"`: sampled number of components.
#'  * `"Kplus"`: number of filled, i.e., non-empty components, for each iteration.
#'  * `"e0"`: sampled Dirichlet parameter of the prior on the weights (if \eqn{e_0} is random).
#'  * `"alpha"`: sampled Dirichlet parameter of the prior on the weights (if \eqn{\alpha} is random).
#'  * `"acc"`: logical vector indicating acceptance in the Metropolis-Hastings step when sampling either \eqn{e_0} or \eqn{\alpha}.
#'
#' @examples
#' if (requireNamespace("poLCA", quietly = TRUE)) {
#'     data("carcinoma", package = "poLCA")
#'     y <- carcinoma
#'     N <- nrow(y)
#'     r <- ncol(y)
#'     
#'     M <- 200
#'     thin <- 1
#'     burnin <- 100
#'     Kmax <- 50  
#'     Kinit <- 10
#'     
#'     G <- "MixDynamic"
#'     priorOnAlpha <- priorOnAlpha_spec("gam_1_2")
#'     priorOnK <- priorOnK_spec("Pois_1")
#'     
#'     cat <- apply(y, 2, max)
#'     a0 <- rep(1, sum(cat))
#' 
#'     cl_y <- kmeans(y, centers = Kinit, iter.max = 20)
#'     S_0 <- cl_y$cluster
#'     eta_0 <- cl_y$size/N
#' 
#'     pi_0 <- do.call("cbind", lapply(1:r, function(j) {
#'         prop.table(table(S_0, y[, j]), 1)
#'     }))
#' 
#'     result <- sampleLCA(
#'         y, S_0, pi_0, eta_0, a0, 
#'         M, burnin, thin, Kmax, 
#'         G, priorOnK, priorOnAlpha)
#' 
#'     K <- result$K
#'     Kplus <- result$Kplus   
#'     
#'     plot(K, type = "l", ylim = c(0, max(K)),  
#'          xlab = "iteration", main = "",
#'          ylab = expression("K" ~ "/" ~ K["+"]), col = 1)
#'     lines(Kplus, col = 2)
#'     legend("topright", legend = c("K", expression(K["+"])),
#'            col = 1:2, lty = 1, box.lwd = 0)
#' }
#' 
sampleLCA <- function(y, S, pi, eta, a0,
                      M, burnin, thin, Kmax,
                      G = c("MixDynamic", "MixStatic"),
                      priorOnK, priorOnWeights,
                      verbose = FALSE) {

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
        e0 <- alpha / nrow(pi)
    } else {
        e0 <- priorOnWeights$param$e0
        alpha <- e0 * nrow(pi)
        log_p_e0 <- priorOnWeights$param$log_p_e0
    }
    s0_proposal <- priorOnWeights$param$s0_proposal

    N <- nrow(y)  # number of observations
    r <- ncol(y)  # number of dimensions
    cat <- apply(y, 2, max)  # number of categories
  
    ## index for taking the variable-specific probabilities from the probability matrix
    index <- c(0, cumsum(cat))
    low <- (index + 1)[-length(index)]
    up <- index[-1]

    ## indicator matrix for categories.
    y_I <- matrix(0, N, sum(cat))
    for (j in 1:r) {
        for (d in 1:cat[j]) {
            y_I[y[, j] == d, low[j]+d-1] <- 1
        }
    }
    
    ## initial values
    eta_j <- eta
    pi_j <- pi
    S_j <- S
    Nk_j <- tabulate(S_j, K_j)
    if (verbose) {
        cat("0 ", Nk_j)
    }
    Kp_j <- sum(Nk_j != 0)  # number of nonempty components
    
    ## generating matrices for storing the draws:
    result <- list(Eta = matrix(NA_real_, M, Kmax),
                   Pi = array(NA_real_, dim = c(M, Kmax, sum(cat))),
                   Nk = matrix(NA_integer_, M, Kmax),
                   S = matrix(NA_integer_, M, N),
                   K = rep(NA_integer_, M),
                   Kp = rep(0L, M),
                   mixlik = rep(0,  M),
                   mixprior = rep(0, M),
                   nonnormpost = rep(0, M),
                   nonnormpost_mode = vector("list", Kmax),
                   e0 = rep(NA_real_, M),
                   alpha = rep(NA_real_, M),
                   acc = rep(NA, M))
    
    ## Initialising the result object
    for (k in 1:Kmax) {
        result$nonnormpost_mode[[k]] <- list(nonnormpost = -(10)^18)
    }
    
    ##---------------------- simulation ----------------------------------------------
    
    m <- 1
    Mmax <- M * thin
    while (m <= Mmax || m <= burnin) {
        if (verbose && !(m%%500)) {
            cat("\n", m, " ", Nk_j)
        }
        
        if (m == burnin) {
            m <- 1
            burnin <- 0
        }

        ## first step: classify observations and determine new partition
        
        mat <- matrix(0, N, K_j)
        for (k in 1:K_j) {
            mat[, k] <- eta_j[k] * (apply(y_I, 1, function(x) (prod((pi_j[k, x != 0])))))
        }
        S_j <- apply(mat, 1, function(x) sample(1:K_j, 1, prob = x, replace = T))
        
        ## determine P
        Nk_j <- tabulate(S_j, K_j)  # length(Nk_j)=K_j
        Kp_j <- sum(Nk_j != 0)
        
        ## reorder the components
        perm <- c(which(Nk_j > 0), which(Nk_j == 0))
        pi_j <- pi_j[perm, , drop = FALSE]  # length(pi_j[,1])=K_j   
        S_ <- rep(FALSE, N)
        for (i in 1:length(perm)) {
            S_[S_j == i] <- which(perm == i)
        }
        S_j <- S_
        Nk_j <- tabulate(S_j, Kp_j)  # length(Nk_j)=Kp_j

        ## second step: parameter simulation conditional on partition P=(N_1,...,N_K+):
        
        ## (2a) update parameters of filled components sample r*K
        ## probability distributions pi_{k,j}:
        Nk_jd <- matrix(0, K_j, sum(cat))
        for (k in 1:K_j) {
            for (j in 1:r) {
                Nk_jd[k, low[j]:up[j]] <- tabulate(y[S_j == k, j], cat[j])
                a_kj <- a0[low[j]:up[j]] + Nk_jd[k, low[j]:up[j]]
                pi_j[k, low[j]:up[j]] <- MCMCpack::rdirichlet(1, a_kj)
            }
        }
    
        ## third step: sample K and alpha (or e0) conditional on partition
        if (G == "MixDynamic") {
            ## (3a) Sample K, if e0=alpha/K (=dependent on K)
            K_j <- sampleK_alpha(Kp_j, Kmax, Nk_j, alpha, log_pK)
            
            ## (3b) Sample alpha, if alpha~p(a_alpha,b_alpha)
            value <- sampleAlpha(N, Nk_j, K_j, alpha, s0_proposal, log_pAlpha)
            alpha <- value$alpha
            e0 <- alpha / K_j
            acc <- value$acc
        } else {
            ## (3a*) Sample K, if e0 fixed or e0~G(a_e,b_e) (independent of K):
            K_j <- sampleK_e0(Kp_j, Kmax, log_pK, log_p_e0, e0, N)
            
            ## (3b*) Sample e0, if e0~G(a_e,b_e) (independent of K)
            value <- sampleE0(K_j, Kp_j, N, Nk_j, s0_proposal, e0, log_p_e0)
            e0 <- value$e0
            alpha <- e0 * K_j
            acc <- value$acc
        }
        
        ## fourth step: add empty components conditional on K
        ## (4a) Add/remove empty components
        if (K_j > Kp_j) {Nk_j <- c(Nk_j[1:Kp_j], rep(0, (K_j - Kp_j)))  # length(Nk_j)=K_j
            pi_j <- rbind(pi_j[1:Kp_j, , drop = FALSE], matrix(0, K_j - Kp_j, sum(cat)))
            
            for (k in (Kp_j + 1):K_j) {
                for (j in 1:r) {
                    pi_j[k, low[j]:up[j]] <- MCMCpack::rdirichlet(1, a0[low[j]:up[j]])
                }
            }
        } else {
            pi_j <- pi_j[1:K_j, , drop = FALSE]
        }

        ## (4b): Sample eta_j:
        ek <- e0 + Nk_j
        eta_j <- MCMCpack::rdirichlet(1, ek)
        
        ## fifth step: evaluating the mixture likelihood and storing the values
        
        ## evaluating the mixture likelihood
        mat_neu <- mat
        mixlik_j <- sum(log(rowSums(mat_neu)))
        
        ## evaluating the mixture prior
        logprior_pi <- 0
        for (j in 1:r) {
            logprior_pi <- logprior_pi +
                DirichletReg::ddirichlet(pi_j[, low[j]:up[j]],
                                         alpha = a0[low[j]:up[j]],
                                         log = TRUE, sum.up = TRUE)
        }
        mixprior_j <- DirichletReg::ddirichlet((eta_j), rep(e0, K_j)) + logprior_pi + log_pK(K_j)
        
        if (burnin == 0) {
            result$mixlik[m] <- mixlik_j
            result$mixprior[m] <- mixprior_j
            result$nonnormpost[m] <- result$mixlik[m] + result$mixprior[m]
        }
        
        ## storing the nonnormalized posterior for having good starting points for clustering the draws in
        ## the point process repres.
        if ((burnin == 0) & (result$nonnormpost[m] > result$nonnormpost_mode[[Kp_j]]$nonnormpost)) {
            result$nonnormpost_mode[[Kp_j]] <- list(nonnormpost = result$nonnormpost[m],
                                                    pi = pi_j[Nk_j != 0, ],
                                                    eta = eta_j)
        }
        
        ## storing the results for given S_j, Nk,K_j
        if ((burnin == 0) & !(m%%thin)) {
            result$Pi[m/thin, 1:K_j, ] <- pi_j
            result$Eta[m/thin, 1:K_j] <- eta_j
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



