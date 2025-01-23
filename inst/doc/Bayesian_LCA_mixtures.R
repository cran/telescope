## -----------------------------------------------------------------------------
library("telescope")

## -----------------------------------------------------------------------------
data("SimData", package = "telescope")
y <- as.matrix(SimData[, 1:30])
z <- SimData[, 31]

## -----------------------------------------------------------------------------
dim(y)

## -----------------------------------------------------------------------------
table(z)

## -----------------------------------------------------------------------------
N <- nrow(y)
r <- ncol(y)
cat <- apply(y, 2, max)

## -----------------------------------------------------------------------------
Mmax <- 300
thin <- 1
burnin <- 100

## -----------------------------------------------------------------------------
M <- Mmax/thin

## -----------------------------------------------------------------------------
Kmax <- 50  
Kinit <- 10

## -----------------------------------------------------------------------------
L <- 2

## -----------------------------------------------------------------------------
G <- "MixDynamic"

## -----------------------------------------------------------------------------
priorOnAlpha <- priorOnAlpha_spec("gam_1_2")

## -----------------------------------------------------------------------------
priorOnK <- priorOnK_spec("BNB_143")

## -----------------------------------------------------------------------------
d0 <- 1  

## -----------------------------------------------------------------------------
a_mu <- rep(20, sum(cat)) 
mu_0 <- matrix(rep(rep(1/cat, cat), Kinit),
               byrow = TRUE, nrow = Kinit)

c_phi <- 30
d_phi <- 1
b_phi <- rep(10, r)

a_phi <- rep(1, r)  
phi_0 <- matrix(cat, Kinit, r, byrow = TRUE)

## regularizing constant
a_00 <- 0.05    

## proposal standard deviations for MH steps for sampling mu and phi, and regularizing constant eps to bound the Dirichlet proposal mu away from the boundary of the simplex
s_mu <- 2  
s_phi <- 2 
eps <- 0.01

## -----------------------------------------------------------------------------
set.seed(1234)
cl_y <- kmeans(y, centers = Kinit, nstart = 100, iter.max = 50)
S_0 <- cl_y$cluster
eta_0 <- cl_y$size/N

## -----------------------------------------------------------------------------
I_0 <- rep(1L, N)
if (L > 1) {
  for (k in 1:Kinit) {
    cl_size <- sum(S_0 == k)
    I_0[S_0 == k] <- rep(1:L, length.out = cl_size)
  }
} 

## -----------------------------------------------------------------------------
index <- c(0, cumsum(cat))
low <- (index + 1)[-length(index)]
up <- index[-1]

pi_km <- array(NA_real_, dim = c(Kinit, L, sum(cat)))
rownames(pi_km) <- paste0("k_", 1:Kinit)
for (k in 1:Kinit) {
    for (l in 1:L) {
        index <- (S_0 == k) & (I_0 == l)
        for (j in 1:r) {
            pi_km[k, l, low[j]:up[j]] <- tabulate(y[index, j], cat[j]) / sum(index)
        }
    }
}
pi_0 <- pi_km 

## -----------------------------------------------------------------------------
result <- sampleLCAMixture(y, S_0, L, 
                           pi_0, eta_0, mu_0, phi_0,
                           a_00, a_mu, a_phi, b_phi, c_phi, d_phi,
                           M, burnin, thin, Kmax,
                           s_mu, s_phi, eps,
                           G, priorOnAlpha, d0, priorOnK)

## -----------------------------------------------------------------------------
Eta <- result$Eta
S <- result$S
K <- result$K
Kplus <- result$Kplus

Nk <- result$Nk
Nl <- result$Nl
acc <- result$acc
e0 <- result$e0
alpha <- result$alpha

Mu <- result$Mu
Phi <- result$Phi
B_phi <- result$B_phi
acc_mu <- result$acc_mu
acc_phi <- result$acc_phi
nonnormpost_mode <- result$nonnormpost_mode
Pi_k <- result$Pi_k

## ----fig.height = 5, fig.width = 7--------------------------------------------
acc <- sum(acc)/M
acc
plot(1:length(alpha), alpha, type = "l", 
     ylim = c(0, max(alpha)), 
     xlab = "iterations", ylab = expression(alpha))

## ----fig.height = 5, fig.width = 7--------------------------------------------
hist(alpha, freq = FALSE, breaks = 50)
mean(alpha)
quantile(alpha, probs = c(0.25, 0.5, 0.75))

## ----fig.height = 5, fig.width = 7--------------------------------------------
hist(e0, freq = FALSE, breaks = 50)
mean(e0)
quantile(e0, probs = c(0.25, 0.5, 0.75))

## -----------------------------------------------------------------------------
k <- 1  # component
j <- 1  # variable
sum(acc_mu[, k, j])/M

## ----fig.height = 5, fig.width = 7--------------------------------------------
boxplot(Mu[, k, seq(2, 2*r, 2)], xlab = "mu_j")

## ----fig.height = 5, fig.width = 7--------------------------------------------
j <- 1  # variable
d <- 2  # category
k <- 1  # component
plot(Mu[, k, low[j]+(d-1)], type = "l", 
     ylab = paste0("mu_jkd, j=", j, ",d=", d," (k=", k, ")"))

## ----fig.height = 5, fig.width = 7--------------------------------------------
k <- 1  # component
j <- 1  # variable
sum(acc_phi[, k, j])/M

boxplot(Phi[, k, ], xlab = "phi_j",
        ylim = quantile(Phi[, k, ], c(0, 0.95)))

j <- 3  # variable
k <- 1  # component
plot(Phi[, k, j], type = "l",
     ylab = paste0("phi_jkd, j=", j, ", (k=", k, ")"))

## ----fig.height = 5, fig.width = 7--------------------------------------------
plot(K, type = "l", ylim = c(0, max(K)), 
     xlab = "iteration", main = "", ylab = "count",
     lwd = 0.5, col = "grey")
points(Kplus, type = "l", col = "red3",
       lwd = 2, lty = 1)
legend("topright", legend = c("K", "K+"), col = c("grey", "red3"),
       lwd = 2)

## ----fig.height = 5, fig.width = 7--------------------------------------------
k <- 3   # component
matplot(Nl[seq(100, M, 10), ((k-1)*L+1):(k*L)],
        type = "l", ylab = "Nl")

## ----fig.height = 5, fig.width = 7--------------------------------------------
p_Kplus <- tabulate(Kplus, nbins = max(Kplus))/M 
barplot(p_Kplus/sum(p_Kplus), xlab = expression(K["+"]), names = 1:length(p_Kplus),
        ylab = expression("p(" ~ K["+"] ~ "|" ~ bold(y) ~ ")"))

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
Nk_Kplus <- (Nk[Kplus == Kplus_hat, ] > 0)

## -----------------------------------------------------------------------------
index <- Kplus == Kplus_hat
Nk[is.na(Nk)] <- 0
Nk_Kplus <- (Nk[Kplus == Kplus_hat, ] > 0)

Pi_k_inter <- Pi_k[index,,] 
Pi_k_Kplus <- array(0, dim = c(M0, Kplus_hat, sum(cat)))  
for (j in 1:sum(cat)) {
  Pi_k_Kplus[, , j] <- Pi_k_inter[, , j][Nk_Kplus]
}

Mu_inter <- Mu[index, , ]       
Mu_Kplus <- array(0, dim = c(M0, Kplus_hat, sum(cat)))  
for (j in 1:sum(cat)) {
  Mu_Kplus[, , j] <- Mu_inter[, , j][Nk_Kplus]
}

Phi_inter <- Phi[index, , ]
Phi_Kplus <- array(0, dim = c(M0, Kplus_hat, r))
for (j in 1:r) {
  Phi_Kplus[, , j] <- Phi_inter[, , j][Nk_Kplus]
}

Eta_inter <- Eta[index, ]
Eta_Kplus <- matrix(Eta_inter[Nk_Kplus], ncol = Kplus_hat)

v <- which(index)
S_Kplus <- matrix(0, M0, N)
for (i in seq_along(v)) {
  m <- v[i]
  perm_S <- rep(0, Kmax)
  perm_S[Nk[m, ] != 0] <- 1:Kplus_hat
  S_Kplus[i, ] <- perm_S[S[m, ]]
}


## -----------------------------------------------------------------------------
Func_init <- nonnormpost_mode[[Kplus_hat]]$pi_k
identified_Kplus <- identifyLCAMixture(
    Pi_k_Kplus, Mu_Kplus, Phi_Kplus, Eta_Kplus, S_Kplus, Func_init)

## -----------------------------------------------------------------------------
identified_Kplus$non_perm_rate

## -----------------------------------------------------------------------------
colMeans(identified_Kplus$Eta)

## -----------------------------------------------------------------------------
z_sp <- apply(identified_Kplus$S, 2,
              function(x) which.max(tabulate(x, Kplus_hat)))
table(z_sp)
library("mclust")
classError(z_sp, z)$errorRate
adjustedRandIndex(z, z_sp)

## ----fig.height = 5, fig.width = 7--------------------------------------------
k <- 1  # component
boxplot(identified_Kplus$Mu[, k, seq(2, 2*r, 2)],
        ylab = "mu_j")
boxplot(identified_Kplus$Phi[, k, ], ylab = "phi_j",
        ylim = quantile(identified_Kplus$Phi[, k, ], c(0, 0.95)))

