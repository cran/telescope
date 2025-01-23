## -----------------------------------------------------------------------------
library("telescope")

## -----------------------------------------------------------------------------
y <- c(34, 24, 22, 17, 17, 15, 15, 14, 12, 12, 11, 11, 10, 10, 9, 9,
       8, 7, 7, 7, 6, 6, 6, 5, 5, 5, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2,
       2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0)
N <- length(y)

## ----fig.height = 5, fig.width = 7--------------------------------------------
barplot(table(y), col = "lightgray", xlab = "count", ylab = "Frequency")   

## -----------------------------------------------------------------------------
Mmax <- 20000
thin <- 1
burnin <- 1000

## -----------------------------------------------------------------------------
M <- Mmax/thin

## -----------------------------------------------------------------------------
Kmax <- 50  
Kinit <- 10

## -----------------------------------------------------------------------------
G <- "MixDynamic"

## -----------------------------------------------------------------------------
priorOnAlpha <- priorOnAlpha_spec("gam_1_2")

## -----------------------------------------------------------------------------
priorOnK <- priorOnK_spec("BNB_143")

## -----------------------------------------------------------------------------
a0 <- 0.1 
h0 <- 0.5 
b0 <- a0/mean(y) 
H0 <- h0/b0

## -----------------------------------------------------------------------------
set.seed(1234) 
cl_y <- kmeans(y, centers = Kinit, nstart = 30)
S_0 <- cl_y$cluster
mu_0 <- t(cl_y$centers)
eta_0 <- rep(1/Kinit, Kinit)

## -----------------------------------------------------------------------------
estGibbs <- samplePoisMixture(
    y, S_0, mu_0, eta_0, 
    a0, b0, h0, H0,
    M, burnin, thin, Kmax, 
    G, priorOnK, priorOnAlpha)

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
plot(1:length(alpha), alpha,
     type = "l", ylim = c(0, max(alpha)),
     xlab = "iterations", ylab = expression(alpha))

## ----fig.height = 5, fig.width = 7--------------------------------------------
hist(alpha, freq = FALSE, breaks = 50)

## ----fig.height = 5, fig.width = 7--------------------------------------------
plot(1:length(e0), e0,
     type = "l", ylim = c(0, max(e0)),
     xlab = "iterations", ylab = expression(e["0"]))
hist(e0, freq = FALSE, breaks = 50)
mean(e0)
quantile(e0, probs = c(0.25, 0.5, 0.75))

## ----fig.height = 5, fig.width = 7--------------------------------------------
plot(1:length(K), K, type = "l", ylim = c(0, 30), 
     xlab = "iteration", main = "", ylab = "number",
     lwd = 0.5, col = "grey")
points(1:length(K), Kplus, type = "l", col = "red3",
       lwd = 2, lty = 1)
legend("topright", legend = c("K", expression(K["+"])),
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
Func_init <- matrix(nonnormpost_mode[[Kplus_hat]]$mean_muk,
                    nrow = Kplus_hat)
identified_Kplus <- identifyMixture(
    Mu_Kplus, Mu_Kplus, Eta_Kplus, S_Kplus, Func_init)
identified_Kplus$non_perm_rate

