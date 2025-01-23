## -----------------------------------------------------------------------------
library("telescope")

## -----------------------------------------------------------------------------
freq <- c(5, 15, 3, 2, 4, 4, 3, 1, 1, 2, 4, 2, 0, 2, 0, 0, 1, 3, 2, 1,
          2, 1, 3, 3, 2, 4, 1, 0, 0, 4, 1, 3, 2, 2, 7, 3)
pattern <- cbind(F = rep(rep(1:3, each = 4), 3),
                 C = rep(1:3, each = 3 * 4),
                 M = rep(1:4, 9))
fear <- pattern[rep(seq_along(freq), freq),]
pi_stern <- matrix(c(0.74, 0.26, 0.0, 0.71, 0.08, 0.21, 0.22, 0.6, 0.12, 0.06, 
                     0.00, 0.32, 0.68, 0.28, 0.31, 0.41, 0.14, 0.19, 0.40, 0.27), 
                   ncol = 10, byrow = TRUE)

## -----------------------------------------------------------------------------
N <- nrow(fear)
r <- ncol(fear)

## ----fig.width=8, fig.height=8------------------------------------------------
plotBubble(fear)

## -----------------------------------------------------------------------------
Mmax <- 5000
thin <- 10
burnin <- 1000

## -----------------------------------------------------------------------------
M <- Mmax/thin

## -----------------------------------------------------------------------------
Kmax <- 50  
Kinit <- 10

## -----------------------------------------------------------------------------
G <- "MixDynamic"

## -----------------------------------------------------------------------------
priorOnAlpha <- priorOnAlpha_spec("F_6_3")

## -----------------------------------------------------------------------------
priorOnK <- priorOnK_spec("BNB_143")

## -----------------------------------------------------------------------------
cat <- apply(fear, 2, max)
a0 <- rep(1, sum(cat))

## -----------------------------------------------------------------------------
set.seed(1234)
if (requireNamespace("klaR", quietly = TRUE)) {
    cl_y <- klaR::kmodes(data = fear, modes = Kinit,
                         iter.max = 20, weighted = FALSE)
} else {
    cl_y <- kmeans(fear, centers = Kinit, iter.max = 20)
}
S_0 <- cl_y$cluster
eta_0 <- cl_y$size/N

## -----------------------------------------------------------------------------
pi_0 <- do.call("cbind", lapply(1:r, function(j) {
    prop.table(table(S_0, fear[, j]), 1)
}))
round(pi_0, 2)

## -----------------------------------------------------------------------------
result <- sampleLCA(
    fear, S_0, pi_0, eta_0, a0, 
    M, burnin, thin, Kmax, 
    G, priorOnK, priorOnAlpha)

## -----------------------------------------------------------------------------
Eta <- result$Eta
Pi <- result$Pi
Nk <- result$Nk
S <- result$S
K <- result$K
Kplus <- result$Kplus   
nonnormpost_mode <- result$nonnormpost_mode
e0 <- result$e0
alpha <- result$alpha
acc <- result$acc

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
plot(1:length(e0), e0, type = "l", 
     ylim = c(0, max(e0)), 
     xlab = "iterations", ylab = expression(e["0"]))

## ----fig.height = 5, fig.width = 7--------------------------------------------
hist(e0, freq = FALSE, breaks = 50)
mean(e0)
quantile(e0, probs = c(0.25, 0.5, 0.75))

## ----fig.height = 5, fig.width = 7--------------------------------------------
plot(1:length(K), K, type = "l", ylim = c(0, 30), 
     xlab = "iteration", main = "", ylab = "number",
     lwd = 0.5, col = "grey")
points(1:length(K), Kplus, type = "l", col = "red3",
       lwd = 2, lty = 1)
legend("topright", legend = c("K", "K+"), col = c("grey", "red3"),
       lwd = 2)

## ----fig.height = 5, fig.width = 7--------------------------------------------
Kplus <- rowSums(Nk != 0, na.rm = TRUE)  
p_Kplus <- tabulate(Kplus, nbins = max(Kplus))/M 
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
Nk_Kplus <- (Nk[Kplus == Kplus_hat, ] > 0)

## -----------------------------------------------------------------------------
Pi_inter <- Pi[index, , ]
Pi_Kplus <- array(0, dim = c(M0, sum(cat), Kplus_hat))
for (j in 1:sum(cat)) {
  Pi_Kplus[, j, ] <- Pi_inter[, , j][Nk_Kplus]
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
Func_init <- nonnormpost_mode[[Kplus_hat]]$pi
identified_Kplus <- identifyMixture(
    Pi_Kplus, Pi_Kplus, Eta_Kplus, S_Kplus, Func_init)

## -----------------------------------------------------------------------------
identified_Kplus$non_perm_rate

## -----------------------------------------------------------------------------
colMeans(identified_Kplus$Eta)

## -----------------------------------------------------------------------------
Pi_identified <- colMeans(identified_Kplus$Mu)
round(t(Pi_identified), digits = 2)
pi_stern

## -----------------------------------------------------------------------------
z_sp <- apply(identified_Kplus$S, 2,
              function(x) which.max(tabulate(x, Kplus_hat)))
table(z_sp)

