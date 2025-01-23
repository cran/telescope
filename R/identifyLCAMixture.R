#' Solve label switching and identify mixture for a mixture of LCA models. 
#'
#' @description Clustering of the draws in the point process representation (PPR) using
#'   \eqn{k}-means clustering.
#' @param Func A numeric array of dimension \eqn{M \times d \times K}; data for clustering in the PPR.
#' @param Mu A numeric array of dimension \eqn{M \times r \times K}; draws of cluster means.
#' @param Eta A numeric array of dimension \eqn{M \times K}; draws of cluster sizes.
#' @param Phi A numeric array of dimension \eqn{M \times K \times r}; draws of precisions.
#' @param S A numeric matrix of dimension \eqn{M \times N}; draws of cluster assignments. 
#' @param centers An integer or a numeric matrix of dimension \eqn{K \times d}; used to initialize [stats::kmeans()].
#' @return A named list containing:
#'  * `"S"`: reordered assignments.
#'  * `"Mu"`: reordered Mu matrix.
#'  * `"Phi"`: reordered Phi matrix. 
#'  * `"Eta"`: reordered weights.
#'  * `"non_perm_rate"`: proportion of draws where the clustering did not
#'    result in a permutation and hence no relabeling could be
#'    performed; this is the proportion of draws discarded.
#' 
#' @details The following steps are implemented:
#' * A functional of the draws of the component-specific
#'   parameters (`Func`) is passed to the function.  The functionals
#'   of each component and iteration are stacked on top of each other in
#'   order to obtain a matrix where each row corresponds to the
#'   functional of one component.
#' * The functionals are clustered into \eqn{K_+} clusters using \eqn{k}-means
#'   clustering. For each functional a group label is obtained.
#' * The obtained labels of the functionals are used to construct
#'   a classification for each MCMC iteration.  Those classifications
#'   which are a permutation of \eqn{(1,\ldots,K_+)} are used to reorder 
#'   the Mu and Eta draws and the assignment matrix S. This results in an
#'   identified mixture model.
#'  * Note that only iterations resulting in permutations
#'   are used for parameter estimation and deriving the final
#'   partition. Those MCMC iterations where the obtained
#'   classifications of the functionals are not a permutation of
#'   \eqn{(1,\ldots,K_+)} are discarded as no unique assignment of functionals
#'   to components can be made.  If the non-permutation rate, i.e. the
#'   proportion of MCMC iterations where the obtained classifications
#'   of the functionals are not a permutation, is high, this is an
#'   indication of a poor clustering solution, as the
#'   functionals are not clearly separated.
#'
identifyLCAMixture <- function(Func, Mu, Phi, Eta, S, centers) {

    K <- length(Func[1, , 1])
    M <- length(Func[, 1, 1])
    d <- length(Func[1, 1, ])
    r <- length(Phi[1, 1, ])
    N <- length(S[1, ])
    
    ##---------------------- step 1
    ## from the functional draws, a matrix KM of size (M*K)xr is
    ## created, by putting the draws of the different clusters below
    ## each other

    KM <- do.call("rbind", lapply(1:K, function(k) matrix(Func[, k, , drop = TRUE], ncol = d)))
    colnames(KM) <- paste0("func", 1:d)

    ##---------------------- step 2
    ## applying k-means clustering with known clusters
    ## centers(=mu_func) and obtaining the classification 'class'

    cl_y <- kmeans(KM, centers = centers, nstart = 30)
    class <- cl_y$cluster
    
    ##---------------------- step 3
    ## classification matrix Rho_m is constructed: each row
    ## corresponds to one MCMC iteration and contains the labels of
    ## the group where the corresponding draw was assigned to by
    ## kmeans.
    Rho_m <- NULL
    for (l in 0:(K - 1)) {
        Rho_m <- cbind(Rho_m, class[(l * M + 1):((l + 1) * M)])
    }

    ##---------------------- step 4
    ## Identifying non-permutations, i.e. those rows where the
    ## sequence of group labels is not a permutation of 1:K
    m_rho <- NULL
    for (m in 1:M) {
        if (any(sort(Rho_m[m, ]) != 1:K))
            m_rho <- c(m_rho, m)
    }
    non_perm_rate <- length(m_rho)/M

    ##---------------------- step 5
    ## Reordering of the draws of Mu, Eta and S using the matrix
    ## Rho_m.  In this way, for the permutations, a unique labeling is
    ## achieved.
    Mu_reord <- array(0, dim = c(M, K, d))
    Phi_reord <- array(0, dim = c(M, K, r))
    Eta_reord <- matrix(0, M, K)
    S_reord <- matrix(0, M, N)
    
    for (m in 1:M) {
        Mu_reord[m,  Rho_m[m, ],] <- Mu[m, , ]
        Phi_reord[m,  Rho_m[m, ],] <- Phi[m, , ]
        Eta_reord[m, Rho_m[m, ]] <- Eta[m, ]
        S_reord[m, ] <- Rho_m[m, ][S[m, ]] 
    }
    
    ##---------------------- step 5
    ## Finally, drop draws which are not permutations:
    Mu_only_perm <- Mu_reord[setdiff(1:M, m_rho), , ]  # will discard any duplicated values from 1:M, i.e values of m_rho
    Phi_only_perm <- Phi_reord[setdiff(1:M, m_rho), , ]
    Eta_only_perm <- Eta_reord[setdiff(1:M, m_rho), ]
    S_only_perm <- S_reord[setdiff(1:M, m_rho), ]
    
    return(list(S = S_only_perm,
                Mu = Mu_only_perm,
                Phi = Phi_only_perm,
                Eta = Eta_only_perm,
                non_perm_rate = non_perm_rate))
}

