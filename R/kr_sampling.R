#' KernelReSample: sampling probabilities for each test point
#'
#' Implements ONLY the sampling (selection) part used in \code{kr_boot_plus}:
#' (1) clustering on training X to get centers,
#' (2) compute center labels ystar (cluster mean of y),
#' (3) tune lambda (and optionally bw_scale) on centers by BIC,
#' (4) for each test point x0: compute center weights and induce train-point probs.
#'
#' This function does NOT fit any estimator (no KRR / no FALKON). It returns
#' sampling probabilities (and optionally centers-level weights) that can be
#' used by any downstream estimator.
#'
#' @param N Integer. Number of training points (must equal nrow(x)).
#' @param x Matrix. Training covariates (N x p).
#' @param y Numeric. Training responses (length N). Used to compute center labels for BIC tuning.
#' @param x0s Matrix. Test covariates (m x p).
#' @param lambda0 Numeric vector. Candidate lambdas for BIC tuning on centers.
#'
#' @param kernel Kernel family: "gaussian" (default), "laplace", "matern".
#' @param matern_nu Matern smoothness: 0.5, 1.5, or 2.5 (only used if kernel="matern").
#'
#' @param clustering Clustering for centers: "kmeans" (default), "kmeans_pp",
#'   "kmeans_rp", or "kmedoids".
#' @param n_centers Integer. Number of centers/clusters used for selection weights (default = 200).
#' @param kmeans_maxiter Integer. Max iterations for kmeans backends.
#' @param rp_dim Integer. Random projection dimension for "kmeans_rp".
#'
#' @param bw_centers_method Bandwidth rule for *centers* (selection weights): "median" (default) or "detcov".
#' @param bw_scope Scope for centers bandwidth when bw_centers_method is used:
#'   "centers" (default) uses centers; "train_sample" uses a subsample of X.
#' @param bw_train_sample Integer. Subsample size if bw_scope="train_sample".
#' @param tune_bw Logical. If TRUE, allow bw_scales for centers bandwidth.
#' @param bw_scales Numeric. Multiplicative scales applied to centers bandwidth in tuning.
#'
#' @param joint_tune Logical. If TRUE and tune_bw=TRUE, joint tune (bw_scale, lambda) by BIC.
#'   If FALSE (default), tune lambda only (bw fixed at base value).
#'
#' @param normalize_cols Logical. If TRUE (default), normalize each column of the returned
#'   sampling matrix \code{finalprobs} to sum to 1 (so each column is a valid probability vector).
#'
#' @param return Character. What to return:
#'   - "probs": N x m sampling probabilities (default).
#'   - "weights": n_centers x m center weights W (normalized by column).
#'   - "both": both probs and weights.
#'   - "prep": return a prep object containing centers/idx/A_W/bw/lambda, useful if you want to reuse for new x0s.
#' @param verbose Logical. Print progress messages.
#'
#' @return A list containing requested outputs and metadata:
#'   - finalprobs (N x m) if return includes probs
#'   - W (n_centers x m) if return includes weights
#'   - idx, centers, lambda_opt, bw_centers, kernel, clustering, etc.
#'   - timing (elapsed seconds for cluster/ystar/tune/weights/total)
#'
#' @examples
#' set.seed(1)
#' N <- 2000
#' p <- 3
#' x <- matrix(rnorm(N * p), N, p)
#' y <- x[, 1] - 2 * x[, 2] + rnorm(N)
#' m <- 10
#' x0s <- matrix(rnorm(m * p), m, p)
#' lambda0 <- 10^seq(-2, 2, length.out = 10)
#'
#' # 1) Return sampling probabilities (N x m)
#' samp <- kr_sampling(
#'   N = N, x = x, y = y, x0s = x0s, lambda0 = lambda0,
#'   kernel = "gaussian",
#'   clustering = "kmeans",
#'   n_centers = 50,
#'   bw_centers_method = "median",
#'   normalize_cols = TRUE,
#'   return = "probs"
#' )
#' dim(samp$finalprobs)               # N x m
#' range(colSums(samp$finalprobs))    # ~ 1
#'
#' # Draw one bootstrap-resample index set for the first test point
#' probs1 <- samp$finalprobs[, 1]
#' idx_sub <- sample.int(N, size = 200, replace = TRUE, prob = probs1)
#' head(idx_sub)
#'
#' # 2) Return both probs and center weights
#' samp2 <- kr_sampling(
#'   N = N, x = x, y = y, x0s = x0s, lambda0 = lambda0,
#'   kernel = "laplace",
#'   clustering = "kmeans_pp",
#'   n_centers = 50,
#'   normalize_cols = TRUE,
#'   return = "both"
#' )
#' dim(samp2$W)          # n_centers x m
#' dim(samp2$finalprobs) # N x m
#'
#' # 3) Return a prep object (centers/tuning results) for reuse
#' prep <- kr_sampling(
#'   N = N, x = x, y = y, x0s = x0s, lambda0 = lambda0,
#'   n_centers = 50,
#'   return = "prep"
#' )
#' names(prep)
#'
#' @export
kr_sampling <- function(
    N, x, y, x0s, lambda0,
    kernel = c("gaussian","laplace","matern"),
    matern_nu = c(0.5, 1.5, 2.5),
    clustering = c("kmeans","kmeans_pp","kmeans_rp","kmedoids"),
    n_centers = 200,
    kmeans_maxiter = 50,
    rp_dim = 24,
    bw_centers_method = c("median", "detcov"),
    bw_scope = c("centers","train_sample"),
    bw_train_sample = 2000,
    tune_bw = FALSE,
    bw_scales = c(0.25, 0.5, 1, 2),
    joint_tune = FALSE,
    normalize_cols = TRUE,
    return = c("probs","weights","both","prep"),
    verbose = FALSE
) {
  # ---- coerce types ----
  x <- as.matrix(x); x0s <- as.matrix(x0s)
  y <- as.numeric(y)

  # ---- checks ----
  if (nrow(x) != N) stop("N must equal nrow(x).")
  if (length(y) != N) stop("y must have length N.")
  if (ncol(x0s) != ncol(x)) stop("x0s must have same number of columns as x.")
  if (!is.numeric(lambda0) || length(lambda0) < 1) stop("lambda0 must be a numeric vector.")
  if (!is.numeric(n_centers) || length(n_centers) != 1 || n_centers < 2) stop("n_centers must be an integer >= 2.")

  kernel <- match.arg(kernel)
  matern_nu <- as.numeric(match.arg(as.character(matern_nu), as.character(c(0.5,1.5,2.5))))
  clustering <- match.arg(clustering)
  bw_centers_method <- match.arg(bw_centers_method)
  bw_scope <- match.arg(bw_scope)
  return <- match.arg(return)

  n_centers <- min(as.integer(n_centers), N)

  timing <- list()
  t_all <- proc.time()

  # =========================
  # Helpers (copied from kr_boot_plus)
  # =========================
  detcov_safe <- function(Z) {
    Z <- as.matrix(Z)
    if (nrow(Z) < 2) return(1)
    C <- stats::cov(Z)
    d <- suppressWarnings(det(C))
    if (!is.finite(d) || d <= 0) return(1)
    d
  }

  median_dist <- function(Z, max_m = 2000L) {
    Z <- as.matrix(Z)
    m <- nrow(Z)
    if (m <= 1) return(1)
    if (m > max_m) Z <- Z[sample.int(m, max_m), , drop = FALSE]
    stats::median(as.numeric(stats::dist(Z)))
  }

  bw_from_points <- function(Z, method, kernel, max_m = 2000L) {
    Z <- as.matrix(Z)
    if (method == "detcov") {
      d <- detcov_safe(Z)
      if (kernel == "gaussian") {
        return(max(d, 1e-12))
      } else {
        p <- ncol(Z)
        return(max(d^(1/(2 * p)), 1e-12))
      }
    } else {
      md <- median_dist(Z, max_m = max_m)
      if (kernel == "gaussian") return(max(md^2, 1e-12))
      return(max(md, 1e-12))
    }
  }

  sqdist <- function(A, B) {
    A <- as.matrix(A); B <- as.matrix(B)
    AA <- rowSums(A^2)
    BB <- rowSums(B^2)
    D2 <- outer(AA, BB, "+") - 2 * (A %*% t(B))
    D2[D2 < 0] <- 0
    D2
  }

  kernel_mat <- function(A, B, bw) {
    D2 <- sqdist(A, B)
    if (kernel == "gaussian") {
      return(exp(-D2 / max(bw, 1e-12)))
    }
    D <- sqrt(D2)
    ell <- max(bw, 1e-12)
    if (kernel == "laplace") {
      return(exp(-D / ell))
    }
    if (matern_nu == 0.5) {
      return(exp(-D / ell))
    } else if (matern_nu == 1.5) {
      s <- sqrt(3) * D / ell
      return((1 + s) * exp(-s))
    } else {
      s <- sqrt(5) * D / ell
      return((1 + s + (s^2)/3) * exp(-s))
    }
  }

  chol_solve_spd <- function(A, b) {
    R <- chol(A)
    backsolve(R, forwardsolve(t(R), b))
  }

  # ---------------------------
  # Clustering backends (copied)
  # ---------------------------
  kmeanspp_seed <- function(X, k) {
    X <- as.matrix(X)
    N0 <- nrow(X)
    centers <- matrix(0, k, ncol(X))
    centers[1, ] <- X[sample.int(N0, 1), ]
    d2 <- rowSums((X - matrix(centers[1, ], N0, ncol(X), byrow = TRUE))^2)
    for (j in 2:k) {
      s <- sum(d2)
      if (!is.finite(s) || s <= 0) {
        idx <- sample.int(N0, 1)
      } else {
        probs <- d2 / s
        idx <- sample.int(N0, 1, prob = probs)
      }
      centers[j, ] <- X[idx, ]
      d2_new <- rowSums((X - matrix(centers[j, ], N0, ncol(X), byrow = TRUE))^2)
      d2 <- pmin(d2, d2_new)
    }
    centers
  }

  do_clustering <- function(X, k) {
    t0 <- proc.time()
    if (clustering == "kmeans") {
      km <- stats::kmeans(X, centers = k, iter.max = kmeans_maxiter)
      out <- list(idx = km$cluster, centers = km$centers)
    } else if (clustering == "kmeans_pp") {
      initC <- kmeanspp_seed(X, k)
      km <- stats::kmeans(X, centers = initC, iter.max = kmeans_maxiter)
      out <- list(idx = km$cluster, centers = km$centers)
    } else if (clustering == "kmeans_rp") {
      p <- ncol(X)
      q <- min(as.integer(rp_dim), p)
      Rr <- matrix(stats::rnorm(p * q), p, q) / sqrt(q)
      Z <- X %*% Rr
      km <- stats::kmeans(Z, centers = k, iter.max = kmeans_maxiter)
      idx <- km$cluster
      centers <- matrix(0, k, p)
      for (j in 1:k) {
        sel <- which(idx == j)
        centers[j, ] <- if (length(sel) == 0) X[sample.int(nrow(X), 1), ] else colMeans(X[sel, , drop = FALSE])
      }
      out <- list(idx = idx, centers = centers)
    } else { # kmedoids
      if (!requireNamespace("cluster", quietly = TRUE)) {
        stop("Package 'cluster' is required for clustering='kmedoids'. Please install it.")
      }
      pam <- cluster::pam(X, k = k)
      out <- list(idx = pam$clustering, centers = X[pam$id.med, , drop = FALSE])
    }
    timing$cluster <<- (proc.time() - t0)[3]
    out
  }

  # ---------------------------
  # Center-stage tuning via BIC (copied)
  # ---------------------------
  tune_lambda_bw_centers <- function(centers, ystar, base_bw, lambda_grid, scales, joint) {
    best <- list(score = Inf, lam = lambda_grid[1], bw = base_bw)

    for (s in scales) {
      bw_try <- base_bw * s
      KC <- kernel_mat(centers, centers, bw_try)
      KC <- (KC + t(KC)) / 2 + 1e-10 * diag(nrow(KC))

      eg <- eigen(KC, symmetric = TRUE)
      evals <- pmax(eg$values, 0)
      Q <- eg$vectors
      Vy <- drop(t(Q) %*% ystar)

      for (lam in lambda_grid) {
        frac <- evals / (evals + lam)
        fc <- Q %*% (frac * Vy)
        resid <- ystar - drop(fc)
        res2 <- max(sum(resid^2), .Machine$double.eps)
        df <- sum(frac)
        bic <- nrow(centers) * log(res2) + log(nrow(centers)) * df
        if (bic < best$score) best <- list(score = bic, lam = lam, bw = bw_try)
      }

      if (!joint) break
    }
    best
  }

  # =========================
  # 1) clustering
  # =========================
  if (verbose) message("[kr_sampling] clustering...")
  cl <- do_clustering(x, n_centers)
  idx <- as.integer(cl$idx)
  centers <- as.matrix(cl$centers)

  # =========================
  # 2) ystar per center (cluster-mean y)
  # =========================
  if (verbose) message("[kr_sampling] center labels (cluster means)...")
  t0 <- proc.time()

  counts <- tabulate(idx, nbins = n_centers)
  sumy <- tapply(y, idx, sum)
  sumy <- as.numeric(sumy)
  sumy[is.na(sumy)] <- 0
  cnt <- counts
  cnt[cnt == 0] <- 1
  ystar <- sumy / cnt

  timing$ystar <- (proc.time() - t0)[3]

  # =========================
  # 3) choose centers bandwidth + lambda (BIC)
  # =========================
  if (verbose) message("[kr_sampling] tuning lambda (and optionally bw_scale) on centers...")
  t0 <- proc.time()

  base_bw <- if (bw_scope == "centers") {
    bw_from_points(centers, bw_centers_method, kernel)
  } else {
    S <- min(bw_train_sample, nrow(x))
    bw_from_points(x[sample.int(nrow(x), S), , drop = FALSE], bw_centers_method, kernel)
  }

  scales <- if (isTRUE(tune_bw)) bw_scales else 1
  tuned <- tune_lambda_bw_centers(
    centers = centers,
    ystar = ystar,
    base_bw = base_bw,
    lambda_grid = lambda0,
    scales = scales,
    joint = isTRUE(joint_tune) && isTRUE(tune_bw)
  )

  lambda_opt <- tuned$lam
  bw_centers <- tuned$bw
  timing$tune <- (proc.time() - t0)[3]

  # =========================
  # 4) weights/probs (centers bw only)
  # =========================
  if (verbose) message("[kr_sampling] computing sampling probabilities...")
  t0 <- proc.time()

  KC <- kernel_mat(centers, centers, bw_centers)
  KC <- (KC + t(KC)) / 2 + 1e-10 * diag(n_centers)

  tau <- 1 / lambda_opt
  A_W <- diag(n_centers) + tau * KC

  Kxc <- kernel_mat(x0s, centers, bw_centers)     # m x n_centers
  Kx0s <- chol_solve_spd(A_W, t(Kxc))             # n_centers x m

  W <- abs(Kx0s)
  csW <- colSums(W); csW[csW <= 0] <- 1
  W <- sweep(W, 2, csW, "/")

  finalweights <- (n_centers / N) * W
  finalprobs <- finalweights[idx, , drop = FALSE]  # N x m

  if (isTRUE(normalize_cols)) {
    cs <- colSums(finalprobs)
    cs[!is.finite(cs) | cs <= 0] <- 1
    finalprobs <- sweep(finalprobs, 2, cs, "/")
  }

  timing$weights <- (proc.time() - t0)[3]
  timing$total <- (proc.time() - t_all)[3]

  # ---- return control ----
  if (return == "prep") {
    return(list(
      idx = idx,
      centers = centers,
      ystar = ystar,
      kernel = kernel,
      matern_nu = if (kernel == "matern") matern_nu else NULL,
      clustering = clustering,
      n_centers = n_centers,
      bw_centers = bw_centers,
      bw_centers_method = bw_centers_method,
      bw_scope = bw_scope,
      lambda_opt = lambda_opt,
      tau = tau,
      A_W = A_W,
      normalize_cols = normalize_cols,
      timing = timing
    ))
  }

  out <- list(
    idx = idx,
    centers = centers,
    kernel = kernel,
    matern_nu = if (kernel == "matern") matern_nu else NULL,
    clustering = clustering,
    n_centers = n_centers,
    bw_centers = bw_centers,
    bw_centers_method = bw_centers_method,
    bw_scope = bw_scope,
    lambda_opt = lambda_opt,
    tau = tau,
    normalize_cols = normalize_cols,
    timing = timing
  )

  if (return == "weights") {
    out$W <- W
    return(out)
  }
  if (return == "both") {
    out$W <- W
    out$finalprobs <- finalprobs
    return(out)
  }
  out$finalprobs <- finalprobs
  out
}
