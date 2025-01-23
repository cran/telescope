## -----------------------------------------------------------------------------
library("telescope")

## -----------------------------------------------------------------------------
data("thyroid", package = "mclust")
y <- as.matrix(thyroid[, 2:6])
z <- thyroid$Diagnosis

## -----------------------------------------------------------------------------
r <- ncol(y)
N <- nrow(y)
table(z)

## ----fig.width=7, fig.height=8------------------------------------------------
pairs(y, col = z)

## -----------------------------------------------------------------------------
Mmax <- 5000
thin <- 10
burnin <- 100

## -----------------------------------------------------------------------------
M <- Mmax/thin

## -----------------------------------------------------------------------------
Kmax <- 50  
Kinit <- 10

## -----------------------------------------------------------------------------
G <- "MixStatic"      
priorOnE0 <- priorOnE0_spec("G_1_20", 1)

## -----------------------------------------------------------------------------
priorOnK <- priorOnK_spec("BNB_143")

## -----------------------------------------------------------------------------
R <- apply(y, 2, function(x) diff(range(x)))
b0 <- apply(y, 2, median)
B_0 <- rep(1, r)  
B0 <- diag((R^2) * B_0)
c0 <- 2.5 + (r-1)/2
g0 <- 0.5 + (r-1)/2
G0 <- 100 * g0/c0 * diag((1/R^2), nrow = r)
C0 <- g0 * chol2inv(chol(G0))

## -----------------------------------------------------------------------------
set.seed(1234)
cl_y <- kmeans(y, centers = Kinit, nstart = 100, iter.max = 30)
S_0 <- cl_y$cluster
mu_0 <- t(cl_y$centers)

## -----------------------------------------------------------------------------
eta_0 <- rep(1/Kinit, Kinit)
Sigma_0 <- array(0, dim = c(r, r, Kinit))
Sigma_0[, , 1:Kinit] <- 0.5 * C0

## -----------------------------------------------------------------------------
estGibbs <- sampleMultNormMixture(
    y, S_0, mu_0, Sigma_0, eta_0,
    c0, g0, G0, C0, b0, B0,
    M, burnin, thin, Kmax,
    G, priorOnK, priorOnE0)

## -----------------------------------------------------------------------------
Mu <- estGibbs$Mu
Eta <- estGibbs$Eta
S <- estGibbs$S
Nk <- estGibbs$Nk
K <- estGibbs$K
Kplus <- estGibbs$Kplus
nonnormpost_mode <- estGibbs$nonnormpost_mode
acc <- estGibbs$acc
e0 <- estGibbs$e0
alpha <- estGibbs$alpha

## ----fig.height = 5, fig.width = 7--------------------------------------------
sum(acc)/M 
plot(1:length(e0), e0,
     type = "l", ylim = c(0, max(e0)),
     xlab = "iterations", ylab = "e0")

## ----fig.height = 5, fig.width = 7--------------------------------------------
hist(e0, freq = FALSE, breaks = 50)
quantile(e0, probs = c(0.25, 0.5, 0.75))

## ----fig.height = 5, fig.width = 7--------------------------------------------
plot(seq_along(K), K, type = "l", ylim = c(0, 30),  
     xlab = "iteration", main = "",
     ylab = expression("K" ~ "/" ~ K["+"]),
     lwd = 0.5, col = "grey")
lines(seq_along(Kplus), Kplus, col = "red3", lwd = 2, lty = 1)
legend("topright", legend=c("K", "K_+"), col=c("grey", "red3"),lwd = 2)

## ----fig.height = 5, fig.width = 7--------------------------------------------
Kplus <- rowSums(Nk != 0, na.rm = TRUE)  
p_Kplus <- tabulate(Kplus, nbins = max(Kplus))/M 
quantile(Kplus, probs = c(0.25, 0.5, 0.75))
barplot(p_Kplus/sum(p_Kplus), xlab = expression(K["+"]), names = 1:length(p_Kplus),
        col = "red3", ylab = expression("p(" ~ K["+"] ~ "|" ~ bold(y) ~ ")"))

## -----------------------------------------------------------------------------
quantile(Kplus, probs = c(0.25, 0.5, 0.75))

## -----------------------------------------------------------------------------
Kplus_hat <- which.max(p_Kplus)
Kplus_hat
M0 <- sum(Kplus == Kplus_hat)
M0          

## -----------------------------------------------------------------------------
p_K <- tabulate(K, nbins = max(K))/M
quantile(K, probs = c(0.25, 0.5, 0.75))

## ----fig.height = 5, fig.width = 7--------------------------------------------
barplot(p_K/sum(p_K), names = 1:length(p_K), xlab = "K", 
        ylab = expression("p(" ~ K ~ "|" ~ bold(y) ~ ")"))
which.max(tabulate(K, nbins = max(K)))   

## -----------------------------------------------------------------------------
index <- Kplus == Kplus_hat
Nk[is.na(Nk)] <- 0
Nk_Kplus <- (Nk[index, ] > 0)

## -----------------------------------------------------------------------------
Mu_inter <- Mu[index, , , drop = FALSE]
Mu_Kplus <- array(0, dim = c(M0, r, Kplus_hat)) 
for (j in 1:r) {
  Mu_Kplus[, j, ] <- Mu_inter[, j, ][Nk_Kplus]
}

Eta_inter <- Eta[index, ]
Eta_Kplus <- matrix(Eta_inter[Nk_Kplus], ncol = Kplus_hat)
Eta_Kplus <- sweep(Eta_Kplus, 1, rowSums(Eta_Kplus), "/")

w <- which(index)
S_Kplus <- matrix(0, M0, N)
for (i in seq_along(w)) {
    m <- w[i]
    perm_S <- rep(0, Kmax)
    perm_S[Nk[m, ] != 0] <- 1:Kplus_hat
    S_Kplus[i, ] <- perm_S[S[m, ]]
}

## -----------------------------------------------------------------------------
Func_init <- t(nonnormpost_mode[[Kplus_hat]]$mu)
identified_Kplus <- identifyMixture(
    Mu_Kplus, Mu_Kplus, Eta_Kplus, S_Kplus, Func_init)

## -----------------------------------------------------------------------------
identified_Kplus$non_perm_rate

## -----------------------------------------------------------------------------
colMeans(identified_Kplus$Mu)
colMeans(identified_Kplus$Eta)

## -----------------------------------------------------------------------------
z_sp <- apply(identified_Kplus$S, 2,
              function(x) which.max(tabulate(x, Kplus_hat)))
table(z_sp)

## -----------------------------------------------------------------------------
table(z, z_sp)
library("mclust")
classError(z_sp, z)$errorRate
adjustedRandIndex(z, z_sp)

## ----fig.width=7, fig.height=8------------------------------------------------
plotScatter(y, z_sp, label = "y", trim = 0)

