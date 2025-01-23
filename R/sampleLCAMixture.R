#' Telescoping sampling of the mixture of LCA models where a prior on the
#' number of components K is specified.
#' 
#' @description
#'   * The MCMC scheme is implemented as suggested in Malsiner-Walli et al (2024).
#'   * Also the priors on the model parameters are specified as in Malsiner-Walli et al (2024),
#'     see the vignette for details and notation.
#'
#' @param y A numeric matrix; containing the data where categories are coded with numbers.
#' @param S A numeric matrix; containing the initial cluster assignments.
#' @param L A numeric scalar; specifiying the number of classes within each component.
#' @param pi A numeric matrix; containing the initial class-specific
#'     occurrence probabilities.
#' @param eta A numeric vector; containing the initial cluster sizes.
#' @param mu A numeric matrix; containing the initial central component
#'     occurrence probabilities.
#' @param phi A numeric matrix; containing the initial component- and variable-specific
#'     precisions.
#' @param a_00 A numeric scalar; specifying the prior parameter a_00.
#' @param a_mu A numeric vector; containing the prior parameter a_mu.
#' @param a_phi A numeric vector; containing the prior parameter a_phi for each variable. 
#' @param b_phi A numeric vector; containing the initial value of b_phi for each variable.
#' @param c_phi A numeric vector; containing the prior parameter c_phi for each variable.
#' @param d_phi A numeric vector; containing the prior parameter d_phi for each variable.   
#' @param M A numeric scalar; specifying the number of recorded
#'     iterations.
#' @param burnin A numeric scalar; specifying the number of burn-in
#'     iterations.
#' @param thin A numeric scalar; specifying the thinning used for the
#'     iterations.
#' @param Kmax A numeric scalar; the maximum number of components. 
#' @param s_mu A numeric scalar; specifying the standard deviation of
#'     the proposal in the Metropolis-Hastings step when sampling mu.
#' @param s_phi A numeric scalar; specifying the standard deviation of
#'     the proposal in the Metropolis-Hastings step when sampling phi.
#' @param eps A numeric scalar; a regularizing constant to bound the
#'     Dirichlet proposal away from the boundary in the
#'     Metropolis-Hastings step when sampling mu.
#' @param G A character string; either `"MixDynamic"` or `"MixStatic"`.
#' @param priorOnWeights A named list; providing the prior on the mixture weights.
#' @param d0 A numeric scalar; containing the Dirichlet prior parameter on the class weights. 
#' @param priorOnK A named list; providing the prior on the number of components K, see [priorOnK_spec()].
#' @param verbose A logical; indicating if some intermediate clustering
#'     results should be printed.
#' @return A named list containing:
#'  * `"Eta"`: sampled weights.
#'  * `"S"`: sampled assignments.
#'  * `"K"`: sampled number of components.
#'  * `"Kplus"`: number of filled, i.e., non-empty components, for each iteration.
#'  * `"Nk"`: number of observations assigned to the different components, for each iteration.
#'  * `"Nl"`: number of observations assigned to the different classes within the components, for each iteration.  
#'  * `"e0"`: sampled Dirichlet parameter of the prior on the weights (if \eqn{e_0} is random).
#'  * `"alpha"`: sampled Dirichlet parameter of the prior on the weights (if \eqn{\alpha} is random).
#'  * `"acc"`: logical vector indicating acceptance in the Metropolis-Hastings step when sampling either \eqn{e_0} or \eqn{\alpha}.
#'  * `"Mu"`: sampled central component occurrence probabilities.
#'  * `"Phi"`: sampled precisions.
#'  * `"acc_mu"`: the acceptance rate in the Metropolis-Hastings step when sampling \eqn{\mu_{k,j}}.
#'  * `"acc_phi"`: the acceptance rate in the Metropolis-Hastings step when sampling \eqn{\phi_{k,j}}.
#'  * `"nonnormpost_mode"`: parameter values corresponding to the mode of the nonnormalized posterior.
#'  * `"Pi_k"`: sampled weighted component occurrence probabilities.
#' 
#' @examples
#' data("SimData", package = "telescope")
#' y <- as.matrix(SimData[, 1:30])
#' z <- SimData[, 31]
#' N <- nrow(y)
#' r <- ncol(y)
#'     
#' M <- 5
#' thin <- 1
#' burnin <- 0
#' Kmax <- 50  
#' Kinit <- 10
#'     
#' G <- "MixDynamic"
#' priorOnAlpha <- priorOnAlpha_spec("gam_1_2")
#' priorOnK <- priorOnK_spec("Pois_1")
#' d0 <- 1  
#' 
#' cat <- apply(y, 2, max)
#' a_mu <- rep(20, sum(cat))
#' mu_0 <- matrix(rep(rep(1/cat, cat), Kinit),
#'   byrow = TRUE, nrow = Kinit)
#'
#' c_phi <- 30; d_phi <- 1; b_phi <- rep(10, r)
#' a_phi <- rep(1, r)
#' phi_0 <- matrix(cat, Kinit, r, byrow = TRUE)
#'
#' a_00 <- 0.05
#'
#' s_mu <- 2; s_phi <- 2; eps <- 0.01 
#'
#' set.seed(1234)
#' cl_y <- kmeans(y, centers = Kinit, nstart = 100, iter.max = 50)
#' S_0 <- cl_y$cluster
#' eta_0 <- cl_y$size/N
#'
#' I_0 <- rep(1L, N)
#' L <- 2
#' for (k in 1:Kinit) {
#'   cl_size <- sum(S_0 == k)
#'   I_0[S_0 == k] <- rep(1:L, length.out = cl_size)
#' }
#'
#' index <- c(0, cumsum(cat))
#' low <- (index + 1)[-length(index)]
#' up <- index[-1]
#'
#' pi_km <- array(NA_real_, dim = c(Kinit, L, sum(cat)))
#' rownames(pi_km) <- paste0("k_", 1:Kinit)
#' for (k in 1:Kinit) {
#'   for (l in 1:L) {
#'     index <- (S_0 == k) & (I_0 == l)
#'     for (j in 1:r) {
#'       pi_km[k, l, low[j]:up[j]] <- tabulate(y[index, j], cat[j]) / sum(index)
#'     }
#'   }
#' }
#' pi_0 <- pi_km 
#' 
#' result <- sampleLCAMixture(
#'     y, S_0, L,
#'     pi_0, eta_0, mu_0, phi_0,
#'     a_00, a_mu, a_phi, b_phi, c_phi, d_phi,
#'     M, burnin, thin, Kmax,
#'     s_mu, s_phi, eps,
#'     G, priorOnAlpha, d0, priorOnK)

sampleLCAMixture <- function(y, S, L, 
                             pi, eta, mu, phi,
                             a_00, a_mu, a_phi, b_phi, c_phi, d_phi,
                             M, burnin, thin, Kmax,
                             s_mu, s_phi, eps,
                             G, priorOnWeights, d0, priorOnK,
                             verbose = FALSE) {

    if (!requireNamespace("invgamma", quietly = TRUE)) {
        stop("Package invgamma is required, please install.")
    }
    
    ## initial number of componens
    Kinit <- length(eta)  

    ## prior on K and weights
    log_pK <- priorOnK$log_pK
    #G <- match.arg(G)
    if (G == "MixDynamic") {
      log_pAlpha <- priorOnWeights$log_pAlpha  
      a_alpha <- priorOnWeights$param$a_alpha        
      b_alpha <- priorOnWeights$param$b_alpha        
      alpha <- priorOnWeights$param$alpha
      e0 <- alpha / Kinit
    } else {
      e0 <- priorOnWeights$param$e0
      alpha <- e0 * Kinit
      log_p_e0 <- priorOnWeights$log_p_e0
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
    
    ## creating parameters:  
    pi_k <- matrix(NA_real_, Kinit, sum(cat)) 
  
    ## initializing parameters: 
    K_j <- Kinit   # initial number of components
    w <- matrix(1/L, Kinit, L) 
    Nk <- tabulate(S, Kinit)
    if (verbose) {
        cat("0 ", Nk)
    }
    Kp <- sum(Nk != 0)  ##number of nonempty components
    cN <- rep(c_phi, r)   
    dN <- rep(d_phi, r)
  
    ## generating matrices for storing the draws:
    result <- list(Eta = matrix(NA_real_, M, Kmax), 
                   Pi_k = array(NA_real_, dim = c(M, Kmax, sum(cat))), 
                   Nk = matrix(NA_integer_, M, Kmax), 
                   Nl = matrix(NA_integer_, M , Kmax * L), 
                   S = matrix(NA_integer_, M, N),
                   K = rep(NA_integer_, M),
                   Kplus = rep(0L, M),
                   Mu = array(NA_real_, dim = c(M, Kmax, sum(cat))),
                   Phi = array(NA_real_, dim = c(M, Kmax, r)),
                   B_phi = matrix(NA_real_,M, r),
                   acc_mu = array(FALSE, dim = c(M, Kmax, r)),
                   acc_phi = array(FALSE, dim = c(M, Kmax, r)),
                   mixlik = rep(0, M), 
                   mixprior = rep(0, M), 
                   nonnormpost = rep(0, M), 
                   nonnormpost_mode = vector("list", Kmax), 
                   e0 = rep(NA_real_, M),
                   alpha = rep(NA_real_, M),
                   acc = rep(NA, M)
                 )
 
    ## Initialising the result object
    for (k in 1:Kmax) {
        result$nonnormpost_mode[[k]] <- list(nonnormpost = -(10)^18)
    }
  
    ##---------------------- simulation ----------------------------------------------

    m <- 1
    Mmax <- M * thin
    while (m <= Mmax | m <= burnin) {
        
        if (verbose && !(m%%100)) {
            cat("\n", m, ": ", "Nk=", Nk)
        }
        
        if (m == burnin) {
            m <- 1
            burnin <- 0
        }
    
        #################### (1a): sample classification S
        ## sample S:
        matk <- array(0, dim = c(N, K_j, L))
        for (k in 1:K_j) {
            matk[, k, ] <- sapply(1:L, function(l) w[k, l] * (apply(y_I, 1, function(x) prod(pi[k,l,x!=0])) ) )
        }
        mat <- sapply(1:K_j, function(k) eta[k] * rowSums(matk[, k, , drop = FALSE]))
        S <- apply(mat, 1, function(x) sample(1:K_j, 1, prob = x))
        
        ## determine partition P
        Nk <- tabulate(S, K_j)     
        Kp <- sum(Nk != 0)
  
        ## reorder the components
        perm <- c(which(Nk > 0), which(Nk == 0))
        pi <- pi[perm,,, drop = FALSE]       
        S_<- rep(NA_integer_, N)
        for (i in 1:length(perm)) {
            S_[S == i] <- which(perm == i)
        }
        S <- S_
        Nk <- tabulate(S, Kp)
    
        ########### second step: parameter simulation conditional on partition P=(N_1,...,N_K+): 
        ## For cluster k=1,...Kp:
        I <- rep(0L, N)
        for (k in 1:Kp) {
        
            #### (2a): Sample I_i:
            ## Classification of the observations within a cluster to the subcomponents:    
            matl <- matk[S == k, k, ]
            if (is.matrix(matl)) {
                I_l <- apply(matl, 1, function(x) sample(1:L, 1, prob = x))
            } else {
                if (L==1) {
                    I_l <- rep(1, sum(S == k))
                }else{
                    I_l <- sample(1:L, 1, prob = matl)
                }
            }
            I[S== k] <- I_l
            
            ## (2b): sample subcomponent weights w_j:
            Nl <- tabulate(I_l, L)
            if (burnin == 0) {
                result$Nl[m/thin, (L * (k - 1) + 1):(L * (k - 1) + L)] <- Nl
            }
            dl <- d0 + Nl
            w[k, ] <- bayesm::rdirichlet(dl)
        
            ## (2c): sample pi_kl:  
            Nl_jd <- matrix(0, L, sum(cat))
            for (l in 1:L) {
                for (j in 1:r) {
                    Nl_jd[l, low[j]:up[j]] <- tabulate(y[(S==k),, drop = FALSE][I_l==l,j], cat[j]) 
                    a_0j <- a_00 + mu[k,low[j]:up[j]]*phi[k,j]
                    a_lj <- a_0j + Nl_jd[l,low[j]:up[j]]
                    pi[k,l,low[j]:up[j]] <- MCMCpack::rdirichlet(1, a_lj)
                }
            }
            pi_k[k, ] <- w[k, ] %*% pi[k,,]   #these are the weighted cluster probabilities, they are calculated for the ppr.  
        
            ## (2d) sample hyperparamers
        
            ## (i) sample cluster- and dimension-specific mu_kj:
            log_fullCond_mu <-function(x, a)
                sum((a-1)*log(x)) +
                    DirichletReg::ddirichlet(pi[k,,low[j]:up[j]],
                                             alpha = a_00 + as.vector(x)*phi[k, j],
                                             log = TRUE, sum.up = TRUE)
        
            for (j in 1:r) {
                mu_j <- mu[k, low[j]:up[j]]                            # current value
                mu_p <- as.vector(MCMCpack::rdirichlet(1, alpha=eps+s_mu*mu_j))  #proposal
                lalpha1 <- log_fullCond_mu(mu_p, a_mu[low[j]:up[j]]) -
                    log_fullCond_mu(mu_j, a_mu[low[j]:up[j]]) +
                    DirichletReg::ddirichlet(t(mu_j),alpha=0.00+s_mu*mu_p,log=T)-  #ATTENTION:0.01+
                    DirichletReg::ddirichlet(t(mu_p),alpha=0.00+s_mu*mu_j,log=T)   #ATTENTION:0.01+
          
                alpha1 <- min(exp(lalpha1), 1)
                alpha2 <- runif(1)
                if (alpha2 <= alpha1) {
                    mu_j <- mu_p
                    if (burnin == 0) {
                        result$acc_mu[m/thin,k,j] <- TRUE
                    }
                    mu[k,low[j]:up[j]] <- mu_j
                }
            }
        
            ## (ii):  sample cluster- and dimension-specific phi_kj:
            ## for phi ~Gamma^{-1}(a,b):
            ## x=phi[j]
            log_fullCond_phi <-function(x, a, b)
                -(a+1)*log(x) - b/x +
                    DirichletReg::ddirichlet(pi[k,,low[j]:up[j]], alpha=a_00+mu_j*x,log = TRUE,sum.up = TRUE)
        
            phi_j <- rep(NA_real_, r)
        
            for (j in 1:r) {
                mu_j <- mu[k, low[j]:up[j]]
                phi_j <- phi[k, j]
                
                lphi_p <- log(phi[k,j]) + rnorm(1, 0, s_phi)
                phi_p <- exp(lphi_p)         #proposed value
                
                ## ratio of the target distribution between the proposed value and the previous value
                ## (the proposal density is symmetric and is cancelled out)
                lalpha1 <- log_fullCond_phi(phi_p,a_phi[j],b_phi[j])-
                    log_fullCond_phi(phi_j,a_phi[j],b_phi[j])+
                    log(phi_p)-log(phi_j)
                alpha1 <- min(exp(lalpha1), 1)
                alpha2 <- runif(1)
                ## the proposed value is accepted with probability alpha2
                if (alpha2 <= alpha1) {
                    phi_j <- phi_p
                    if (burnin == 0) {
                        result$acc_phi[m/thin,k,j] <- TRUE
                    }
                }
                phi[k,j] <- phi_j
            }  
        }  # end k=1,..,Kp 
    
        ########### third step: sample b_phi
        cN <- c_phi + Kp*a_phi
        for (j in 1:r) {
            dN[j] <- d_phi + sum(1/phi[1:Kp,j])
        }
        b_phi <- rgamma(r, cN, rate = dN)
    
        ########### fourth step: sample K and alpha (or e0) conditional on partition
        if (G == "MixDynamic") {
            ## (3a) Sample K, if e0=alpha/K (=dependent on K)
            K_j <- sampleK_alpha(Kp, Kmax, Nk, alpha, log_pK)
      
            ## (3b) Sample alpha, if alpha~p(a_alpha,b_alpha)
            value <- sampleAlpha(N, Nk, K_j, alpha, s0_proposal, log_pAlpha)
            alpha <- value$alpha
            e0 <- alpha / K_j
            acc <- value$acc
        } else {
            ## (3a*) Sample K, if e0 fixed or e0~G(a_e,b_e) (independent of K):
            K_j <- sampleK_e0(Kp, Kmax, log_pK, log_p_e0, e0, N)
            
            ## (3b*) Sample e0, if e0~G(a_e,b_e) (independent of K)
            value <- sampleE0(K_j, Kp, N, Nk, s0_proposal, e0, log_p_e0)
            e0 <- value$e0
            acc <- value$acc
        }

        ########### fifth step: add empty components conditional on K
        ## (4a) Add/remove empty components    
        if (K_j > Kp) {
            Nk <- c(Nk[1:Kp], rep(0, (K_j-Kp))) 
            phi_e <- matrix(NA, K_j-Kp, r)
            phi <- rbind(phi[1:Kp,], phi_e)
            mu_e <- matrix(NA, K_j-Kp, sum(cat))
            mu <- rbind(mu[1:Kp,],mu_e)
            pi_e <- array(NA, dim = c(K_j-Kp, L, sum(cat)))
            pi <- abind::abind(pi[1:Kp,,, drop = FALSE],pi_e, along = 1)
            pi_k_e <- matrix(NA, K_j-Kp,sum(cat))
            pi_k <- rbind(pi_k[1:Kp,],pi_k_e)
            w <- rbind(w[1:Kp,],
                      MCMCpack::rdirichlet(K_j-Kp, rep(d0,L)))
            
            for(k in (Kp+1):K_j){
                for (j in 1:r) {
                    phi[k,j] <- invgamma::rinvgamma(1, shape=a_phi[j], rate = b_phi[j])
                    mu[k,low[j]:up[j]] <- as.vector(MCMCpack::rdirichlet(1, alpha = a_mu[low[j]:up[j]]))
                    a_0j=a_00+mu[k,low[j]:up[j]]*phi[k,j]
                    
                    pi[k,,low[j]:up[j]] <- MCMCpack::rdirichlet(L,a_0j)  
                }
            }
            
        } else {
            pi <- pi[1:K_j,,, drop = FALSE]
            mu <- mu[1:K_j,, drop = FALSE]
            phi <- phi[1:K_j,, drop = FALSE]
            pi_k <- pi_k[1:K_j,, drop = FALSE]
        } 
        
        ## (4b): Sample eta_j:
        ek <- e0 + Nk
        eta <- MCMCpack::rdirichlet(1,ek)
    
        ########### sixth step: evaluating the mixture likelihood and storing the values

        ## evaluating  the mixture likelihood:
        mixlik_j <- sum(log(rowSums(mat)))
    
        ## (5c) storing the mixture likelihood #nonnormalized posterior
        if (burnin == 0) {
            result$nonnorm_post[m] <- mixlik_j #+ mixprior_j
            result$mixlik[m] <- mixlik_j
            result$nonnormpost[m] <- result$mixlik[m] #+ result$mixprior[m]
        }
    
        ## storing the nonnormalized posterior for having good starting points for clustering the draws in the point process repres. 
        if ((burnin == 0) & (result$nonnormpost[m] > result$nonnormpost_mode[[Kp]]$nonnormpost)) {
            result$nonnormpost_mode[[Kp]] <- list(nonnormpost = result$nonnormpost[m],
                                                  pi_k = pi_k[Nk != 0,],eta = eta[Nk !=0])
        }
    
        ## storing the results for given S,Nk, K_j:
        if ((burnin == 0) & !(m%%thin)) {
            result$Eta[m/thin, 1:K_j] <- eta
            result$S[m/thin, ] <- S
            result$Nk[m/thin, 1:K_j] <- Nk
            
            result$K[m/thin] <- K_j
            result$Kplus[m/thin] <- Kp
            
            result$e0[m/thin] <- e0
            result$alpha[m/thin] <- alpha
            result$acc[m/thin] <- acc
            
            result$Mu[m/thin, 1:K_j, ] <- mu
            result$Phi[m/thin, 1:K_j, ] <- phi
            result$B_phi[m/thin, ] <- b_phi
            result$Pi_k[m/thin, 1:K_j, ] <- pi_k
    }
    m <- m + 1
  }
  return(result)
}




