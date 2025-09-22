# --- Face demo: SVD, NMF, CUR on a 128x128 grayscale matrix ------------------
# Expected input: a numeric 128x128 matrix in one of:
#   images/img_small.rds      (contains the matrix)
#   images/img_small.RData    (object named A or img_small)
#   images/img_small.csv      (comma-separated)
#   images/img_small.txt      (whitespace-separated)
# Output: figures/face_svd_k.pdf, figures/face_nmf_k.pdf, figures/face_cur_k.pdf,
#         tables/face_errors.csv

# -------------------------- setup --------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("images",  showWarnings = FALSE, recursive = TRUE)

# helper: load matrix from multiple formats
load_img_matrix <- function(base = "images/img_small") {
  if (file.exists(paste0(base, ".rds"))) {
    A <- readRDS(paste0(base, ".rds"))
  } else if (file.exists(paste0(base, ".RData"))) {
    e <- new.env()
    load(paste0(base, ".RData"), envir = e)
    # try common names
    nm <- intersect(ls(e), c("A","img_small","img","face","X"))
    if (length(nm) == 0) stop("RData found but no object named A/img_small/img/face/X")
    A <- get(nm[1], envir = e)
  } else if (file.exists(paste0(base, ".csv"))) {
    A <- as.matrix(read.csv(paste0(base, ".csv"), header = FALSE))
  } else if (file.exists(paste0(base, ".txt"))) {
    A <- as.matrix(read.table(paste0(base, ".txt"), header = FALSE))
  } else {
    stop("No image matrix found. Put one of: .rds/.RData/.csv/.txt at ", base, ".*")
  }
  if (!is.matrix(A)) A <- as.matrix(A)
  storage.mode(A) <- "double"
  if (any(is.na(A))) stop("Matrix contains NA values.")
  A
}

A <- load_img_matrix()

# enforce 128x128 (warn if not)
if (!all(dim(A) == c(128,128))) {
  warning(sprintf("Matrix is %dx%d; continuing anyway.", nrow(A), ncol(A)))
}

# normalize to [0,1] for display (does not change relative errors)
rng <- range(A)
if (diff(rng) > 0) {
  A01 <- (A - rng[1]) / (rng[2] - rng[1])
} else {
  A01 <- A
}

# plotting helper (base R)
plot_matrix <- function(M, fname = NULL, main = "", zlim = c(0,1)) {
  if (!is.null(fname)) {
    png(filename = fname, width = 960, height = 960, res = 150)
  }
  par(mar = c(1.2,1,2,1))
  image(t(apply(M, 2, rev)),
        axes = FALSE, zlim = zlim, col = gray((0:255)/255), asp = 1,
        main = main)
  if (!is.null(fname)) dev.off()
}

# -------------------------- SVD ----------------------------------------------
sv <- svd(A, nu = 128, nv = 128)
rank_k <- function(sv, k) sv$u[,1:k,drop=FALSE] %*% (sv$d[1:k] * t(sv$v[,1:k,drop=FALSE]))
fro <- function(M) sqrt(sum(M*M))

k_list <- c(4,8,16,32,64)
err_svd <- data.frame(method="SVD", k=k_list, fro_error=NA_real_, rel_error=NA_real_)

# multi-panel PDF with reconstructions
pdf("figures/face_svd_k.pdf", width = 8, height = 8)
par(mfrow = c(2,3), mar = c(1.2,1,2,1))
plot_matrix(A01, main = "Original")
for (i in seq_along(k_list)) {
  k <- k_list[i]
  Ak <- rank_k(sv, k)
  # scale to [0,1] for display
  rngk <- range(Ak)
  Ak01 <- if (diff(rngk) > 0) (Ak - rngk[1])/(rngk[2]-rngk[1]) else Ak
  plot_matrix(Ak01, main = paste0("SVD k=", k))
  err <- fro(A - Ak)
  err_svd$fro_error[i] <- err
}
par(mfrow = c(1,1))
dev.off()
err_svd$rel_error <- err_svd$fro_error / fro(A)

# -------------------------- NMF (Lee-Seung) ----------------------------------
nmf_available <- requireNamespace("NMF", quietly = TRUE)
err_nmf <- data.frame(method="NMF", k=k_list, fro_error=NA_real_, rel_error=NA_real_)
if (nmf_available) {
  # ensure nonnegativity (shift if necessary)
  minA <- min(A)
  Ashift <- if (minA < 0) A - minA else A
  for (i in seq_along(k_list)) {
    k <- k_list[i]
    fit <- NMF::nmf(Ashift, rank = k, method = "lee", .options = "t") # deterministic-ish
    W <- NMF::basis(fit); H <- NMF::coef(fit)
    Ak <- W %*% H
    if (minA < 0) Ak <- Ak + minA
    # display
    rngk <- range(Ak)
    Ak01 <- if (diff(rngk) > 0) (Ak - rngk[1])/(rngk[2]-rngk[1]) else Ak
    if (i == 1) {
      pdf("figures/face_nmf_k.pdf", width = 8, height = 8)
      par(mfrow = c(2,3), mar = c(1.2,1,2,1))
      plot_matrix(A01, main = "Original")
    }
    plot_matrix(Ak01, main = paste0("NMF k=", k))
    if (i == length(k_list)) { par(mfrow = c(1,1)); dev.off() }
    err_nmf$fro_error[i] <- fro(A - Ak)
  }
  err_nmf$rel_error <- err_nmf$fro_error / fro(A)
} else {
  warning("Package 'NMF' not installed. Skipping NMF section.")
  err_nmf <- err_nmf[0,]  # empty
}

# -------------------------- CUR (simple randomized) --------------------------
set.seed(123)
cur_one <- function(A, k) {
  # leverage-like probs via squared column/row norms
  cn <- colSums(A^2); rn <- rowSums(A^2)
  pc <- cn / sum(cn); pr <- rn / sum(rn)
  cidx <- sample(ncol(A), k, prob = pc, replace = FALSE)
  ridx <- sample(nrow(A), k, prob = pr, replace = FALSE)
  C <- A[, cidx, drop=FALSE]
  R <- A[ridx, , drop=FALSE]
  U <- MASS::ginv(A[ridx, cidx, drop=FALSE])  # pseudo-inverse of intersection
  list(Ahat = C %*% U %*% R, cidx=cidx, ridx=ridx)
}

err_cur <- data.frame(method="CUR", k=k_list, fro_error=NA_real_, rel_error=NA_real_)
pdf("figures/face_cur_k.pdf", width = 8, height = 8)
par(mfrow = c(2,3), mar = c(1.2,1,2,1))
plot_matrix(A01, main = "Original")
for (i in seq_along(k_list)) {
  k <- k_list[i]
  fit <- cur_one(A, k)
  Ak <- fit$Ahat
  rngk <- range(Ak)
  Ak01 <- if (diff(rngk) > 0) (Ak - rngk[1])/(rngk[2]-rngk[1]) else Ak
  plot_matrix(Ak01, main = paste0("CUR k=", k))
  err_cur$fro_error[i] <- fro(A - Ak)
}
par(mfrow = c(1,1))
dev.off()
err_cur$rel_error <- err_cur$fro_error / fro(A)

# -------------------------- summary table ------------------------------------
errs <- rbind(err_svd, err_nmf, err_cur)
write.csv(errs, "tables/face_errors.csv", row.names = FALSE)
message("Done. Figures in 'figures/', errors in 'tables/face_errors.csv'.")
