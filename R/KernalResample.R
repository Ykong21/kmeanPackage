#' KernelReSample
#'
#' Advanced kernel resampling estimator.
#' Provides configurable (i) clustering to build selection weights and
#' (ii) kernel family + bandwidth rules for both center-stage selection and
#' bootstrap-resample KRR estimation.
#'
#' Key design
#' - Compute centers once (clustering on full training X).
#' - Compute a center-stage bandwidth (rho/ell) for *selection weights* only.
#' - Tune lambda by default (BIC on centers); optionally joint tune (bw_scale, lambda).
#' - For each test point x0: compute sampling probs using centers.
#' - For each bootstrap draw: resample indices with replacement and recompute
#'   the *subset* bandwidth on that resample, then do accurate KRR via Cholesky.
#'
#' Bandwidth units are handled cleanly by kernel type:
#' - Gaussian: exp(-||x-y||^2 / rho)     -> rho is squared-distance scale
#' - Laplace:  exp(-||x-y|| / ell)       -> ell is distance scale
#' - Matern:   uses distance/ell         -> ell is distance scale
#'
#' @param N Integer. Number of training points (must equal nrow(x)).
#' @param n Integer. Bootstrap resample size per draw.
#' @param B Integer. Number of bootstrap draws.
#' @param x Matrix. Training covariates (N x p).
#' @param y Numeric. Training responses (length N).
#' @param x0s Matrix. Test covariates (m x p).
#' @param y0s Numeric. Test responses (length m), used only for MSE reporting.
#' @param lambda0 Numeric vector. Candidate lambdas for BIC tuning on centers.
#'
#' @param kernel Kernel family: "gaussian" (default), "laplace", "matern".
#' @param matern_nu Matern smoothness: 0.5, 1.5, or 2.5 (only used if kernel="matern").
#'
#' @param clustering Clustering for centers: "kmeans" (default), "kmeans_pp",
#'   "kmeans_rp", or "kmedoids".
#' @param n_centers Integer. Number of centers/clusters used for selection weights (default = n).
#' @param kmeans_maxiter Integer. Max iterations for kmeans backends.
#' @param rp_dim Integer. Random projection dimension for "kmeans_rp".
#'
#' @param bw_centers_method Bandwidth rule for *centers* (selection weights):
#'   "detcov" (default) or "median".
#' @param bw_scope Scope for centers bandwidth when bw_centers_method is used:
#'   "centers" (default) uses centers; "train_sample" uses a subsample of X.
#' @param bw_train_sample Integer. Subsample size if bw_scope="train_sample".
#' @param tune_bw Logical. If TRUE, allow bw_scales for centers bandwidth.
#' @param bw_scales Numeric. Multiplicative scales applied to centers bandwidth in tuning.
#'
#' @param bw_subset_method Bandwidth rule for *subset* bandwidth on each bootstrap resample:
#'   "detcov" (default) or "median".
#'
#' @param joint_tune Logical. If TRUE and tune_bw=TRUE, joint tune (bw_scale, lambda) by BIC.
#'   If FALSE (default), tune lambda only (bw fixed at base value).
#' @param verbose Logical. Print progress messages.
#'
#' @return A list with fields:
#'   - y0s, est (m x 1), MSE (length m), Ave_MSE
#'   - lambda_opt, bw_centers, bw_centers_method, bw_subset_method
#'   - clustering, kernel, matern_nu (if applicable)
#'   - timing (named list of elapsed seconds)
#'
#' @examples
#' set.seed(1)
#' N <- 2000
#' p <- 3
#' x <- matrix(rnorm(N * p), N, p)
#' y <- x[, 1] - 2 * x[, 2] + rnorm(N)
#' x0s <- matrix(rnorm(10 * p), 10, p)
#' y0s <- x0s[, 1] - 2 * x0s[, 2]
#' lambda0 <- 10^seq(-2, 2, length.out = 10)
#' out <- kr_boot_plus(
#'   N = N, n = 200, B = 5,
#'   x = x, y = y,
#'   x0s = x0s, y0s = y0s,
#'   lambda0 = lambda0,
#'   kernel = "gaussian",
#'   bw_centers_method = "median",
#'   bw_subset_method = "median"
#' )
#' out$Ave_MSE
#'
#' @export
kr_boot_plus <- function(
    N, n, B, x, y, x0s, y0s, lambda0,
    kernel = c("gaussian","laplace","matern"),
    matern_nu = c(0.5, 1.5, 2.5),
    clustering = c("kmeans","kmeans_pp","kmeans_rp","kmedoids"),
    n_centers = n,
    kmeans_maxiter = 50,
    rp_dim = 24,

    # bandwidth (CENTERS): used only for selection weights
    bw_centers_method = c("median", "detcov"),
    bw_scope = c("centers","train_sample"),
    bw_train_sample = 2000,
    tune_bw = FALSE,
    bw_scales = c(0.25, 0.5, 1, 2),

    # bandwidth (SUBSET): recomputed on each bootstrap resample
    bw_subset_method = c("median","detcov"),

    # tuning
    joint_tune = FALSE,
    verbose = FALSE
) {
  # ---- coerce types ----
  x <- as.matrix(x); x0s <- as.matrix(x0s)
  y <- as.numeric(y); y0s <- as.numeric(y0s)

  # ---- checks ----
  if (nrow(x) != N) stop("N must equal nrow(x).")
  if (length(y) != N) stop("y must have length N.")
  if (n <= 0 || n > N) stop("n must be in {1,...,N}.")
  if (B <= 0) stop("B must be positive.")
  if (ncol(x0s) != ncol(x)) stop("x0s must have same number of columns as x.")
  if (length(y0s) != nrow(x0s)) stop("y0s must have length nrow(x0s).")
  if (!is.numeric(lambda0) || length(lambda0) < 1) stop("lambda0 must be a numeric vector.")
  if (!is.numeric(n_centers) || length(n_centers) != 1 || n_centers < 2) stop("n_centers must be an integer >= 2.")

  kernel <- match.arg(kernel)
  matern_nu <- as.numeric(match.arg(as.character(matern_nu), as.character(c(0.5,1.5,2.5))))
  clustering <- match.arg(clustering)
  bw_centers_method <- match.arg(bw_centers_method)
  bw_subset_method <- match.arg(bw_subset_method)
  bw_scope <- match.arg(bw_scope)

  n_centers <- min(as.integer(n_centers), N)

  timing <- list()
  t_all <- proc.time()

  # =========================
  # Helpers
  # =========================

  # det(cov) with fallback to avoid NA/Inf/non-positive
  detcov_safe <- function(Z) {
    Z <- as.matrix(Z)
    if (nrow(Z) < 2) return(1)
    C <- stats::cov(Z)
    d <- suppressWarnings(det(C))
    if (!is.finite(d) || d <= 0) return(1)
    d
  }

  # median of pairwise distances (distance units)
  median_dist <- function(Z, max_m = 2000L) {
    Z <- as.matrix(Z)
    m <- nrow(Z)
    if (m <= 1) return(1)
    if (m > max_m) Z <- Z[sample.int(m, max_m), , drop = FALSE]
    stats::median(as.numeric(stats::dist(Z)))
  }

  # Bandwidth (rho/ell) from points:
  # - Gaussian needs squared-distance scale rho
  # - Laplace/Matern need distance scale ell
  bw_from_points <- function(Z, method, kernel, max_m = 2000L) {
    Z <- as.matrix(Z)
    if (method == "detcov") {
      d <- detcov_safe(Z)
      if (kernel == "gaussian") {
        # original logic: rho in squared-distance scale
        return(max(d, 1e-12))
      } else {
        # convert volume-like det(cov) to distance scale:
        # det(cov) ~ (scale^2)^p -> scale ~ det(cov)^(1/(2p))
        p <- ncol(Z)
        return(max(d^(1/(2 * p)), 1e-12))
      }
    } else {
      md <- median_dist(Z, max_m = max_m)
      if (kernel == "gaussian") return(max(md^2, 1e-12))
      return(max(md, 1e-12))
    }
  }

  # Pairwise squared distances via ||a||^2 + ||b||^2 - 2 a'b
  sqdist <- function(A, B) {
    A <- as.matrix(A); B <- as.matrix(B)
    AA <- rowSums(A^2)
    BB <- rowSums(B^2)
    D2 <- outer(AA, BB, "+") - 2 * (A %*% t(B))
    D2[D2 < 0] <- 0
    D2
  }

  # Kernel matrix with proper bandwidth semantics per kernel
  kernel_mat <- function(A, B, bw) {
    D2 <- sqdist(A, B)
    if (kernel == "gaussian") {
      # exp(-||x-y||^2 / rho)
      return(exp(-D2 / max(bw, 1e-12)))
    }
    D <- sqrt(D2)
    ell <- max(bw, 1e-12)
    if (kernel == "laplace") {
      # exp(-||x-y|| / ell)
      return(exp(-D / ell))
    }
    # matern
    if (matern_nu == 0.5) {
      return(exp(-D / ell))
    } else if (matern_nu == 1.5) {
      s <- sqrt(3) * D / ell
      return((1 + s) * exp(-s))
    } else { # 2.5
      s <- sqrt(5) * D / ell
      return((1 + s + (s^2)/3) * exp(-s))
    }
  }

  # SPD solve via Cholesky: (t(R)R)x=b (chol returns upper R)
  chol_solve_spd <- function(A, b) {
    R <- chol(A)
    backsolve(R, forwardsolve(t(R), b))
  }

  # ---------------------------
  # Clustering backends
  # ---------------------------

  # Simple kmeans++ seeding
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
      R <- matrix(stats::rnorm(p * q), p, q) / sqrt(q)
      Z <- X %*% R
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
  # Center-stage tuning via BIC
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
  if (verbose) message("[kr_boot_plus] clustering...")
  cl <- do_clustering(x, n_centers)
  idx <- as.integer(cl$idx)
  centers <- as.matrix(cl$centers)

  # =========================
  # 2) ystar per center (cluster-mean y)
  # =========================
  if (verbose) message("[kr_boot_plus] center labels (cluster means)...")
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
  if (verbose) message("[kr_boot_plus] tuning lambda (and optionally bw_scale) on centers...")
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
  if (verbose) message("[kr_boot_plus] computing sampling probabilities...")
  t0 <- proc.time()

  KC <- kernel_mat(centers, centers, bw_centers)
  KC <- (KC + t(KC)) / 2 + 1e-10 * diag(n_centers)

  tau <- 1 / lambda_opt
  A_W <- diag(n_centers) + tau * KC

  # Kxc is m x n_centers
  Kxc <- kernel_mat(x0s, centers, bw_centers)
  # Solve A_W^{-1} * Kxc^T -> n_centers x m
  Kx0s <- chol_solve_spd(A_W, t(Kxc))

  W <- abs(Kx0s)
  cs <- colSums(W); cs[cs <= 0] <- 1
  W <- sweep(W, 2, cs, "/")

  # cluster-to-point mapping
  finalweights <- (n_centers / N) * W
  finalprobs <- finalweights[idx, , drop = FALSE]  # N x m

  timing$weights <- (proc.time() - t0)[3]

  # =========================
  # 5) bootstrap estimation (subset bw recomputed per resample)
  # =========================
  if (verbose) message("[kr_boot_plus] bootstrap KRR...")
  t0 <- proc.time()

  m <- nrow(x0s)
  fx0s <- matrix(0, B, m)

  for (e in 1:m) {
    probs <- finalprobs[, e]
    probs[!is.finite(probs) | probs < 0] <- 0
    s <- sum(probs)
    if (!is.finite(s) || s <= 0) probs <- rep(1 / N, N) else probs <- probs / s

    for (b in 1:B) {
      # WITH replacement (matches your original)
      XIND <- sample.int(N, n, replace = TRUE, prob = probs)
      xsub <- x[XIND, , drop = FALSE]
      ysub <- y[XIND]

      # subset bandwidth recomputed on THIS resample
      bw_sub <- bw_from_points(xsub, bw_subset_method, kernel)

      Kss <- kernel_mat(xsub, xsub, bw_sub)
      Kss <- (Kss + t(Kss)) / 2 + 1e-10 * diag(n)

      A <- diag(n) + tau * Kss
      rhs <- tau * ysub
      a <- chol_solve_spd(A, rhs)

      k0 <- kernel_mat(xsub, matrix(x0s[e, ], 1, ncol(xsub)), bw_sub) # n x 1
      fx0s[b, e] <- drop(t(k0) %*% a)
    }
  }

  est <- matrix(colMeans(fx0s), m, 1)
  mse <- (y0s - est[, 1])^2
  Ave_MSE <- mean(mse)

  timing$est <- (proc.time() - t0)[3]
  timing$total <- (proc.time() - t_all)[3]

  list(
    y0s = y0s,
    est = est,
    MSE = mse,
    Ave_MSE = Ave_MSE,
    lambda_opt = lambda_opt,
    bw_centers = bw_centers,
    bw_centers_method = bw_centers_method,
    bw_subset_method = bw_subset_method,
    clustering = clustering,
    kernel = kernel,
    matern_nu = if (kernel == "matern") matern_nu else NULL,
    timing = timing
  )
}


