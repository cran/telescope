% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampleMultNormMixture.R
\name{sampleMultNormMixture}
\alias{sampleMultNormMixture}
\title{Telescoping sampling of a Bayesian finite multivariate Gaussian
mixture where a prior on the number of components is specified.}
\usage{
sampleMultNormMixture(
  y,
  S,
  mu,
  Sigma,
  eta,
  c0,
  g0,
  G0,
  C0,
  b0,
  B0,
  M,
  burnin,
  thin,
  Kmax,
  G = c("MixDynamic", "MixStatic"),
  priorOnK,
  priorOnWeights,
  verbose = FALSE
)
}
\arguments{
\item{y}{A numeric matrix; containing the data.}

\item{S}{A numeric matrix; containing the initial cluster
assignments.}

\item{mu}{A numeric matrix; containing the initial cluster-specific
mean values.}

\item{Sigma}{A numeric matrix; containing the initial cluster-specific
variance covariance values.}

\item{eta}{A numeric vector; containing the initial cluster sizes.}

\item{c0}{A numeric vector; hyperparameter of the prior on \eqn{\Sigma_k}.}

\item{g0}{A numeric vector; hyperparameter of the prior on \eqn{C_0}.}

\item{G0}{A numeric vector; hyperparameter of the prior on \eqn{C_0}.}

\item{C0}{A numeric vector; initial value of the hyperparameter \eqn{C_0}.}

\item{b0}{A numeric vector; hyperparameter of the prior on \eqn{\mu_k}.}

\item{B0}{A numeric vector; hyperparameter of the prior on \eqn{\mu_k}.}

\item{M}{A numeric scalar; specifying the number of recorded
iterations.}

\item{burnin}{A numeric scalar; specifying the number of burn-in
iterations.}

\item{thin}{A numeric scalar; specifying the thinning used for the
iterations.}

\item{Kmax}{A numeric scalar; the maximum number of components.}

\item{G}{A character string; either \code{"MixDynamic"} or \code{"MixStatic"}.}

\item{priorOnK}{A named list; providing the prior on the number of components K, see \code{\link[=priorOnK_spec]{priorOnK_spec()}}.}

\item{priorOnWeights}{A named list; providing the prior on the mixture weights.}

\item{verbose}{A logical; indicating if some intermediate clustering
results should be printed.}
}
\value{
A named list containing:
\itemize{
\item \code{"Mu"}: sampled component means.
\item \code{"Eta"}: sampled weights.
\item \code{"S"}: sampled assignments.
\item \code{"Nk"}: number of observations assigned to the different components, for each iteration.
\item \code{"K"}: sampled number of components.
\item \code{"Kplus"}: number of filled, i.e., non-empty components, for each iteration.
\item \code{"e0"}: sampled Dirichlet parameter of the prior on the weights (if \eqn{e_0} is random).
\item \code{"alpha"}: sampled Dirichlet parameter of the prior on the weights (if \eqn{\alpha} is random).
\item \code{"acc"}: logical vector indicating acceptance in the Metropolis-Hastings step when sampling either \eqn{e_0} or \eqn{\alpha}.
}
}
\description{
\itemize{
\item The MCMC scheme is implemented as suggested in Frühwirth-Schnatter et al (2021).
\item The priors on the model parameters are specified as in Frühwirth-Schnatter et al (2021),
see the vignette for details and notation.
\item The parameterizations of the Wishart and inverse Wishart distribution are used as in
Frühwirth-Schnatter et al (2021), see also the vignette.
}
}
\examples{
y <- iris[, 1:4]
z <- iris$Species
r <- ncol(y)

M <- 50
thin <- 1
burnin <- 0
Kmax <- 40  
Kinit <- 10

G <- "MixStatic"      
priorOnE0 <- priorOnE0_spec("G_1_20", 1)
priorOnK <- priorOnK_spec("BNB_143")

R <- apply(y, 2, function(x) diff(range(x)))
b0 <- apply(y, 2, median)
B_0 <- rep(1, r)  
B0 <- diag((R^2) * B_0)
c0 <- 2.5 + (r-1)/2
g0 <- 0.5 + (r-1)/2
G0 <- 100 * g0/c0 * diag((1/R^2), nrow = r)
C0 <- g0 * chol2inv(chol(G0))

cl_y <- kmeans(y, centers = Kinit, nstart = 100)
S_0 <- cl_y$cluster
mu_0 <- t(cl_y$centers)

eta_0 <- rep(1/Kinit, Kinit)
Sigma_0 <- array(0, dim = c(r, r, Kinit))
Sigma_0[, , 1:Kinit] <- 0.5 * C0

result <- sampleMultNormMixture(
  y, S_0, mu_0, Sigma_0, eta_0,
  c0, g0, G0, C0, b0, B0,  
  M, burnin, thin, Kmax, G, priorOnK, priorOnE0)

K <- result$K
Kplus <- result$Kplus   

plot(K, type = "l", ylim = c(0, max(K)),
     xlab = "iteration", main = "",
     ylab = expression("K" ~ "/" ~ K["+"]), col = 1)
lines(Kplus, col = 2)
legend("topright", legend = c("K", expression(K["+"])),
       col = 1:2, lty = 1, box.lwd = 0)

}
