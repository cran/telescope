## -----------------------------------------------------------------------------
library("telescope")

## -----------------------------------------------------------------------------
y <- c( 9.172,  9.350,  9.483,  9.558,  9.775, 10.227,
       10.406, 16.084, 16.170, 18.419, 18.552, 18.600,
       18.927, 19.052, 19.070, 19.330, 19.343, 19.349,
       19.440, 19.473, 19.529, 19.541, 19.547, 19.663,
       19.846, 19.856, 19.863, 19.914, 19.918, 19.973,
       19.989, 20.166, 20.175, 20.179, 20.196, 20.215,
       20.221, 20.415, 20.629, 20.795, 20.821, 20.846,
       20.875, 20.986, 21.137, 21.492, 21.701, 21.814,
       21.921, 21.960, 22.185, 22.209, 22.242, 22.249,
       22.314, 22.374, 22.495, 22.746, 22.747, 22.888,
       22.914, 23.206, 23.241, 23.263, 23.484, 23.538,
       23.542, 23.666, 23.706, 23.711, 24.129, 24.285,
       24.289, 24.368, 24.717, 24.990, 25.633, 26.960,
       26.995, 32.065, 32.789, 34.279)
N <- length(y)

## ----fig.height = 5, fig.width = 7--------------------------------------------
hist(y, breaks = 50, main = "", col = "lightgray")   

## -----------------------------------------------------------------------------
Mmax <- 10000
thin <- 1
burnin <- 100

## -----------------------------------------------------------------------------
M <- Mmax/thin

## -----------------------------------------------------------------------------
Kmax <- 50  
Kinit <- 10

## -----------------------------------------------------------------------------
G <- "MixStatic" 
priorOnE0 <- priorOnE0_spec("e0const", 1)

## -----------------------------------------------------------------------------
priorOnK <- priorOnK_spec("Unif", 30)

## -----------------------------------------------------------------------------
r <- 1    # dimension
R <- diff(range(y))
c0 <- 2 + (r-1)/2
C0 <- diag(c(0.02*(R^2)), nrow = r)
g0 <- 0.2 + (r-1) / 2
G0 <- diag(10/(R^2), nrow = r)
B0 <- diag((R^2), nrow = r)
b0 <- as.matrix((max(y) + min(y))/2, ncol = 1)  

## -----------------------------------------------------------------------------
set.seed(1234)
cl_y <- kmeans(y, centers = Kinit, nstart = 30)
S_0 <- cl_y$cluster
mu_0 <- t(cl_y$centers)

## -----------------------------------------------------------------------------
eta_0 <- rep(1/Kinit, Kinit)
sigma2_0 <- array(0, dim = c(1, 1, Kinit))
sigma2_0[1, 1, ] <- 0.5 * C0

## -----------------------------------------------------------------------------
estGibbs <- sampleUniNormMixture(
    y, S_0, mu_0, sigma2_0, eta_0,
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
nonnormpost_mode_list <- estGibbs$nonnormpost_mode_list
acc <- estGibbs$acc
e0 <- estGibbs$e0
alpha <- estGibbs$alpha  

## ----fig.height = 5, fig.width = 7--------------------------------------------
plot(seq_along(K), K, type = "l", ylim = c(0, 30),  
     xlab = "iteration", main = "",
     ylab = expression("K" ~ "/" ~ K["+"]),
     lwd = 0.5, col = "grey")
lines(seq_along(Kplus), Kplus, col = "red3", lwd = 2, lty = 1)
legend("topright", legend = c("K", "K+"),
       col = c("grey", "red3"), lwd = 2)

## ----fig.height = 5, fig.width = 7--------------------------------------------
Kplus <- rowSums(Nk != 0, na.rm = TRUE)  
p_Kplus <- tabulate(Kplus, nbins = max(Kplus))/M  
barplot(p_Kplus/sum(p_Kplus), names = 1:length(p_Kplus), 
        col = "red3", xlab = expression(K["+"]),
        ylab = expression("p(" ~ K["+"] ~ "|" ~ bold(y) ~ ")"))

## -----------------------------------------------------------------------------
quantile(Kplus, probs = c(0.25, 0.5, 0.75))

## -----------------------------------------------------------------------------
Kplus_hat <- which.max(p_Kplus)
Kplus_hat
M0 <- sum(Kplus == Kplus_hat)
M0          

## -----------------------------------------------------------------------------
p_K <- tabulate(K, nbins = max(K, na.rm = TRUE))/M
round(p_K[1:20], digits = 2)

## ----fig.height = 5, fig.width = 7--------------------------------------------
barplot(p_K/sum(p_K), names = 1:length(p_K), xlab = "K", 
        ylab = expression("p(" ~ K ~ "|" ~ bold(y) ~ ")"))

## -----------------------------------------------------------------------------
which.max(tabulate(K, nbins = max(K)))
mean(K)
quantile(K, probs = c(0.25, 0.5, 0.75))

## -----------------------------------------------------------------------------
index <- Kplus == Kplus_hat
Nk[is.na(Nk)] <- 0
Nk_Kplus <- (Nk[index, ] > 0)

## -----------------------------------------------------------------------------
Mu_inter <- Mu[index, , , drop = FALSE]
Mu_Kplus <- array(0, dim = c(M0, 1, Kplus_hat)) 
Mu_Kplus[, 1, ] <- Mu_inter[, 1, ][Nk_Kplus]

Eta_inter <- Eta[index, ]
Eta_Kplus <- matrix(Eta_inter[Nk_Kplus], ncol = Kplus_hat)
Eta_Kplus <- sweep(Eta_Kplus, 1, rowSums(Eta_Kplus), "/")

w <- which(index)
S_Kplus <- matrix(0, M0, length(y))
for (i in seq_along(w)) {
    m <- w[i]
    perm_S <- rep(0, Kmax)
    perm_S[Nk[m, ] != 0] <- 1:Kplus_hat
    S_Kplus[i, ] <- perm_S[S[m, ]]
}

## -----------------------------------------------------------------------------
Func_init <- nonnormpost_mode_list[[Kplus_hat]]$mu  
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

## ----fig.height = 5, fig.width = 7--------------------------------------------
hist(y, breaks = 50)
points(cbind(y, rep(0, N)), col = z_sp, lwd = 2)

