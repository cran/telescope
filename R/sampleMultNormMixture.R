#' Telescoping sampling of a Bayesian finite multivariate Gaussian
#' mixture where a prior on the number of components is specified.
#'
#' @description
#'   * The MCMC scheme is implemented as suggested in Frühwirth-Schnatter et al (2021). 
#'   * The priors on the model parameters are specified as in Frühwirth-Schnatter et al (2021),
#'     see the vignette for details and notation.
#'   * The parameterizations of the Wishart and inverse Wishart distribution are used as in 
#'     Frühwirth-Schnatter et al (2021), see also the vignette.
#'
#' @param y A numeric matrix; containing the data.
#' @param S A numeric matrix; containing the initial cluster
#'     assignments.
#' @param mu A numeric matrix; containing the initial cluster-specific
#'     mean values.
#' @param Sigma A numeric matrix; containing the initial cluster-specific
#'     variance covariance values.
#' @param eta A numeric vector; containing the initial cluster sizes.
#' @param c0 A numeric vector; hyperparameter of the prior on \eqn{\Sigma_k}.
#' @param g0 A numeric vector; hyperparameter of the prior on \eqn{C_0}.
#' @param G0 A numeric vector; hyperparameter of the prior on \eqn{C_0}.
#' @param C0 A numeric vector; initial value of the hyperparameter \eqn{C_0}.
#' @param b0 A numeric vector; hyperparameter of the prior on \eqn{\mu_k}.
#' @param B0 A numeric vector; hyperparameter of the prior on \eqn{\mu_k}.
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
#'   * `"Mu"`: sampled component means.
#'   * `"Eta"`: sampled weights.
#'   * `"S"`: sampled assignments.
#'   * `"Nk"`: number of observations assigned to the different components, for each iteration.
#'   * `"K"`: sampled number of components.
#'   * `"Kplus"`: number of filled, i.e., non-empty components, for each iteration.
#'   * `"e0"`: sampled Dirichlet parameter of the prior on the weights (if \eqn{e_0} is random).
#'   * `"alpha"`: sampled Dirichlet parameter of the prior on the weights (if \eqn{\alpha} is random).
#'   * `"acc"`: logical vector indicating acceptance in the Metropolis-Hastings step when sampling either \eqn{e_0} or \eqn{\alpha}.
#' @examples
#' y <- iris[, 1:4]
#' z <- iris$Species
#' r <- ncol(y)
#' 
#' M <- 50
#' thin <- 1
#' burnin <- 0
#' Kmax <- 40  
#' Kinit <- 10
#' 
#' G <- "MixStatic"      
#' priorOnE0 <- priorOnE0_spec("G_1_20", 1)
#' priorOnK <- priorOnK_spec("BNB_143")
#' 
#' R <- apply(y, 2, function(x) diff(range(x)))
#' b0 <- apply(y, 2, median)
#' B_0 <- rep(1, r)  
#' B0 <- diag((R^2) * B_0)
#' c0 <- 2.5 + (r-1)/2
#' g0 <- 0.5 + (r-1)/2
#' G0 <- 100 * g0/c0 * diag((1/R^2), nrow = r)
#' C0 <- g0 * chol2inv(chol(G0))
#' 
#' cl_y <- kmeans(y, centers = Kinit, nstart = 100)
#' S_0 <- cl_y$cluster
#' mu_0 <- t(cl_y$centers)
#' 
#' eta_0 <- rep(1/Kinit, Kinit)
#' Sigma_0 <- array(0, dim = c(r, r, Kinit))
#' Sigma_0[, , 1:Kinit] <- 0.5 * C0
#' 
#' result <- sampleMultNormMixture(
#'   y, S_0, mu_0, Sigma_0, eta_0,
#'   c0, g0, G0, C0, b0, B0,  
#'   M, burnin, thin, Kmax, G, priorOnK, priorOnE0)
#' 
#' K <- result$K
#' Kplus <- result$Kplus   
#' 
#' plot(K, type = "l", ylim = c(0, max(K)),
#'      xlab = "iteration", main = "",
#'      ylab = expression("K" ~ "/" ~ K["+"]), col = 1)
#' lines(Kplus, col = 2)
#' legend("topright", legend = c("K", expression(K["+"])),
#'        col = 1:2, lty = 1, box.lwd = 0)
#'

sampleMultNormMixture <-
    function(y, S, mu, Sigma, eta, c0, g0, G0, C0, 
             b0, B0, M, burnin, thin, Kmax,
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
        log_pAlpha <- priorOnWeights$param$log_pAlpha  
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
    invB0 <- chol2inv(chol(B0))
    eta_j <- eta
    mu_j <- mu
    Sigma_j <- Sigma
    cholSigma_j <- invSigma_j <- array(0, dim = c(r, r, K_j))
    det_invSigma_j <- rep(0, K_j)
    for (k in 1:K_j) {
        cholSigma_j[, , k] <- chol(Sigma_j[, , k])
        invSigma_j[, , k] <- chol2inv(cholSigma_j[, , k])
        det_invSigma_j[k] <- 1 / prod(diag(cholSigma_j[, , k]))^2
    }
    S_j <- S
    C0_j <- C0
    Nk_j <- tabulate(S_j, K_j)
    if (verbose) {
        cat("0 ", Nk_j)
    }
    Kp_j <- sum(Nk_j != 0)  ##number of nonempty components
    acc <- FALSE
    
    ## generating matrices for storing the draws:
    result <- list(Eta = matrix(NA_real_, M, Kmax), 
                   Mu = array(NA_real_, dim = c(M, r, Kmax)), 
                   Nk = matrix(NA_integer_, M, Kmax), 
                   S = matrix(NA_integer_, M, N),
                   K = rep(NA_integer_, M),
                   Kplus = rep(0L, M),
                   mixlik = rep(0, M), 
                   mixprior = rep(0, M), 
                   nonnormpost = rep(0, M), 
                   nonnormpost_mode = vector("list", Kmax), 
                   C0 = array(NA_real_, dim = c(M, r, r)),
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
    while (m <= Mmax | m <= burnin) {
        if (verbose && !(m%%500)) {
            cat("\n", m, " ", Nk_j)
        }
        
        if (m == burnin) {
            m <- 1
            burnin <- 0
        }
    
        ## first step: classify observations and determine new partition 
        mat <- sapply(1:K_j, function(k)
            eta_j[k] * dmvnorm(y, mu_j[, k], as.matrix(Sigma_j[, , k])))
        S_j <- apply(mat, 1, function(x) sample(1:K_j, 1, prob = x, replace = TRUE))
        
        ##determine partition
        Nk_j <- tabulate(S_j, K_j)          #length(Nk_j)=K_j
        Kp_j <- sum(Nk_j != 0)
        
        ##reorder the components
        perm <- c(which(Nk_j > 0), which(Nk_j == 0))
        mu_j <- mu_j[, perm, drop = FALSE]
        Sigma_j <- Sigma_j[,,perm, drop = FALSE] 
        S_ <- rep(NA_integer_, N)
        for (i in 1:length(perm)) {
            S_[S_j == i] <- which(perm == i)
        }
        S_j <- S_
        Nk_j <- tabulate(S_j, Kp_j)     
        
        ## second step: update parameters conditional on partition C=(N_1,...,N_K+): 
        
        ## (2a) update parameters of filled components
    
        ## (i): sample Sigma^{-1} for filled components
        Ck <- array(0, dim = c(r, r, Kp_j))
        ck <- c0 + Nk_j/2
        for (k in 1:Kp_j) {
            Ck[, , k] <- C0_j + 0.5 * crossprod(sweep(y[S_j == k, , drop = FALSE],
                                                      2, mu_j[,k], FUN = "-"))
            sig <- bayesm::rwishart(2 * ck[k], 0.5 * chol2inv(chol(Ck[, , k])))
            Sigma_j[, , k] <- sig$IW
            invSigma_j[, , k] <- sig$W
            det_invSigma_j[k] <- det(as.matrix(invSigma_j[, , k]))
        }

        ## (ii): sample mu_j for filled components
        mean_yk <- matrix(sapply(1:Kp_j, function(k) colMeans(y[S_j == k, , drop = FALSE])),
                          ncol = Kp_j)
        Bk <- array(0, dim = c(r, r, Kp_j))
        bk <- matrix(0, r, Kp_j)
        for (k in 1:Kp_j) {
            Bk[, , k] <- chol2inv(chol(invB0 + as.matrix(invSigma_j[, , k]) * Nk_j[k]))
            bk[, k] <- as.matrix(Bk[, , k]) %*%
                (invB0 %*% b0 + as.matrix(invSigma_j[, , k]) %*% mean_yk[, k] * Nk_j[k])
            mu_j[, k] <- t(chol(as.matrix(Bk[, , k]))) %*% rnorm(r) + bk[, k]
        }
    
        ## (2b) sample hyperparameters conditional on partition
        ## (i): sample C0
        gK <- g0 + Kp_j*c0
        C0_j <- bayesm::rwishart(2 * gK,
                                 0.5 * chol2inv(chol(G0 + rowSums(invSigma_j[,, 1:Kp_j, drop = FALSE], 
                                                                  dims = 2))))$W
        
        ## third step: sample K and alpha (or e0) conditional on partition
     
        if (G == "MixDynamic") {
            ## (3a) Sample K, if e0=alpha/K (=dependent on K)
            K_j <- sampleK_alpha(Kp_j, Kmax, Nk_j, alpha, log_pK)
            
            ## (3b) Sample alpha, if alpha~p(a_alpha,b_alpha) 
            value <- sampleAlpha(N,Nk_j,K_j,alpha,s0_proposal,log_pAlpha)
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
        if (K_j > Kp_j) {
        
            Nk_j <- c(Nk_j[1:Kp_j],rep(0,(K_j-Kp_j)))        #length(Nk_j)=K_j
            
            invSigma_e <- array(NA, dim = c(r, r, K_j-Kp_j))
            invSigma_j <- abind::abind(invSigma_j[,, 1:Kp_j, drop = FALSE],
                                       invSigma_e, along = 3)
            
            Sigma_e <- array(NA, dim = c(r, r, K_j-Kp_j))
            Sigma_j <- abind::abind(Sigma_j[,, 1:Kp_j, drop = FALSE],
                                    Sigma_e, along = 3)
            
            mu_j <- cbind(mu_j[, 1:Kp_j, drop = FALSE],
                          matrix(0, r, K_j - Kp_j))
            
            det_invSigma_j <- c(det_invSigma_j[1:Kp_j],
                                rep(0, K_j - Kp_j))
            
            for(k in (Kp_j+1):K_j) {          
                sig <- bayesm::rwishart(2 * c0, 0.5 * chol2inv(chol(C0_j)))
                                        #attention: rwishart(nu,v)(Rossi)=> nu=2*c0,v=0.5*C0, wishart(c0,C0) (FS)
                Sigma_j[, , k] <- sig$IW
                invSigma_j[, , k] <- sig$W
                det_invSigma_j[k] <- det(as.matrix(invSigma_j[, , k]))
                mu_j[, k] <- t(chol((B0))) %*% rnorm(r) + b0
            } 
            
        }else{
            invSigma_j <- invSigma_j[,, 1:K_j, drop = FALSE]
            Sigma_j <- Sigma_j[,, 1:K_j, drop = FALSE]
            mu_j <- mu_j[, 1:K_j, drop = FALSE]
            det_invSigma_j <- det_invSigma_j[1:K_j]
        }
        
        ## (4b): Sample eta_j:
        ek <- e0 + Nk_j
        eta_j <- MCMCpack::rdirichlet(1,ek) 

        ## fifth step: evaluating the mixture likelihood and storing the values
    
        ## evaluating the mixture likelihood:
        mat_neu <- sapply(1:K_j, function(k) eta_j[k] * dmvnorm(y, mu_j[, k], as.matrix(Sigma_j[, , k])))
        mixlik_j <- sum(log(rowSums(mat_neu)))
        
        ## evaluating the mixture prior:
        mixprior_j <- log(MCMCpack::ddirichlet(as.vector(eta_j), rep(e0, K_j))) +
            sum(dmvnorm(t(mu_j),b0,  B0, log = TRUE))
        sum(sapply(1:K_j, function(k) bayesm::lndIWishart(2 * c0,0.5 * C0_j, as.matrix(Sigma_j[, , k])))) +
            bayesm::lndIWishart(2 * g0, 0.5 * G0, C0_j) +    #ATTENTION!
            log_pK(K_j)
        
        if (burnin == 0) {
            result$mixlik[m] <- mixlik_j
            result$mixprior[m] <- mixprior_j
            result$nonnormpost[m] <- result$mixlik[m] + result$mixprior[m]
        }
        
        ## storing the nonnormalized posterior for having good starting points
        ##    when clustering the draws in the point process representation 
        if ((burnin == 0) && (result$nonnormpost[m] > result$nonnormpost_mode[[Kp_j]]$nonnormpost)) {
            result$nonnormpost_mode[[Kp_j]] <- list(nonnormpost = result$nonnormpost[m],
                                                    mu = mu_j[,Nk_j != 0],
                                                    Sigma = Sigma_j[,,Nk_j != 0],
                                                    bk = bk, 
                                                    Bk = Bk, 
                                                    eta = eta_j[Nk_j != 0])
        }
        
        
        ## storing the draws
        if ((burnin == 0)&!(m%%thin)) {
            result$Mu[m/thin, , 1:K_j] <- mu_j
            result$Eta[m/thin, 1:K_j] <- eta_j
            result$S[m/thin, ] <- S_j
            result$Nk[m/thin, 1:K_j] <- Nk_j
            result$K[m/thin] <- K_j
            result$Kplus[m/thin] <- Kp_j
            result$C0[m/thin,,] <- C0_j
            result$e0[m/thin] <- e0
            result$alpha[m/thin] <- alpha
            result$acc[m/thin] <- acc
        }
        m <- m + 1
    }
    
    return(result)
}


