# install.packages("MASS")  # if missing
library(MASS)

# --- Load image matrix (CSV with header) ---
# Replace the path with your own if needed
A <- as.matrix(read.csv("C:/Users/Karol Dziedziul/Dropbox/matrix/lecture5/img_small.csv", header = TRUE))

# Basic sanity check: expected size 128x128
stopifnot(nrow(A) == 128, ncol(A) == 128)

# Display the original image (flip rows for correct orientation)
image(apply(A, 2, rev), axes = FALSE, col = gray.colors(256), asp = 1)

# Target ranks
ks <- c(4, 8, 16, 32)

# --- SVD rank-k approximation ---
svd_k <- function(M, k) {
  s <- svd(M)
  s$u[, 1:k, drop = FALSE] %*% diag(s$d[1:k]) %*% t(s$v[, 1:k, drop = FALSE])
}

par(mfrow = c(1, length(ks)), mar = c(1, 1, 2, 1))
for (k in ks) {
  A_svd <- svd_k(A, k)
  image(A_svd[nrow(A_svd):1, ], axes = FALSE, col = gray.colors(256), asp = 1,
        main = paste("SVD: k =", k))
}

# --- Simple NMF (works for nonnegative matrices) ---
# Note: NMF requires A >= 0. If your data has negatives, shift/scale before using NMF.
library(NMF)

par(mfrow = c(1, length(ks)), mar = c(1, 1, 2, 1))
for (k in ks) {
  fit <- nmf(A, rank = k, method = "brunet", nrun = 5, seed = 123)
  W <- basis(fit)
  H <- coef(fit)
  Ak <- W %*% H
  image(apply(Ak, 2, rev), axes = FALSE, col = gray.colors(256), asp = 1,
        main = paste("NMF: k =", k))
}

# --- CUR with leverage-score style (deterministic top-c/top-r) ---
cur <- function(M, k, c = 4 * k, r = 4 * k) {
  s  <- svd(M)
  Uk <- s$u[, 1:k, drop = FALSE]
  Vk <- s$v[, 1:k, drop = FALSE]
  J  <- order(rowSums(Vk^2), decreasing = TRUE)[1:c]  # select columns
  I  <- order(rowSums(Uk^2), decreasing = TRUE)[1:r]  # select rows
  C  <- M[, J, drop = FALSE]
  R  <- M[I, , drop = FALSE]
  W  <- M[I, J, drop = FALSE]
  C %*% ginv(W) %*% R
}

par(mfrow = c(1, length(ks)), mar = c(1, 1, 2, 1))
for (k in ks) {
  A_cur <- cur(A, k)
  image(apply(A_cur, 2, rev), axes = FALSE, col = gray.colors(256), asp = 1,
        main = paste("CUR: k =", k))
}
