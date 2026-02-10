#' K-means kernel estimator with resampling
#' Implements a K-means based kernel estimation procedure with a resampling
#' scheme to obtain predictions at query points. The function follows the
#' original implementation, including tuning over lambda0 and the
#' weighted resampling step.
#'
#' @param N Integer. Number of rows in `x` (usually `nrow(x)`).
#' @param n Integer. Subsample size and number of k-means clusters.
#' @param B Integer. Number of bootstrap/resampling replicates.
#' @param x Numeric matrix/data.frame of size N x p. Training covariates.
#' @param y Numeric vector of length N. Training response.
#' @param x0s Numeric matrix of size m x p. Query points where to predict.
#' @param y0s Numeric vector of length m. True responses at `x0s` (for MSE reporting).
#' @param lambda0 Numeric vector. Candidate lambda values for BIC tuning.
#'
#' @return A list with elements:
#' \describe{
#'   \item{coeffi}{Last fitted coefficient vector `avec` (from the last resample).}
#'   \item{y0s}{The provided `y0s`.}
#'   \item{est}{m x 1 matrix of estimated mean predictions at `x0s`.}
#'   \item{MSEKMean}{Vector of squared errors (y0s - est)^2.}
#'   \item{Ave_MSE}{Mean squared error across x0s.}
#'   \item{optimal_lam}{Selected lambda (scalar).}
#'   \item{Converted_Time_Table}{Character vector of elapsed times for Total/Tuning/Est.}
#'   \item{x0s}{The provided `x0s`.}
#' }
#'
#' @examples
#' set.seed(1)
#' N <- 200
#' p <- 3
#' x <- matrix(rnorm(N * p), N, p)
#' y <- x[, 1] - 2 * x[, 2] + rnorm(N)
#' x0s <- matrix(rnorm(10 * p), 10, p)
#' y0s <- x0s[, 1] - 2 * x0s[, 2]
#' out <- kr_boot(
#'   N = N, n = 20, B = 5,
#'   x = x, y = y,
#'   x0s = x0s, y0s = y0s,
#'   lambda0 = 10^seq(-2, 2, length.out = 10)
#' )
#' out$Ave_MSE
#'
#' @export
kr_boot <- function(N, n, B, x, y, x0s, y0s, lambda0) {
  
  # ---- coerce types (non-invasive) ----
  x <- as.matrix(x)
  x0s <- as.matrix(x0s)
  y <- as.numeric(y)
  y0s <- as.numeric(y0s)
  
  # ---- basic checks (as requested) ----
  if (nrow(x) != N) stop("N must equal nrow(x).")
  if (length(y) != N) stop("y must have length N.")
  if (n <= 0 || n > N) stop("n must be in {1,...,N}.")
  if (B <= 0) stop("B must be positive.")
  if (ncol(x0s) != ncol(x)) stop("x0s must have same number of columns as x.")
  if (length(y0s) != nrow(x0s)) stop("y0s must have length nrow(x0s).")
  if (!is.numeric(lambda0) || length(lambda0) < 1) stop("lambda0 must be a numeric vector of candidate lambdas.")
  
  # local trace helper (original code uses tr())
  tr <- function(A) sum(diag(A))
  
  GaussianKernel <- function(x1, x2, rho) {
    Kxy <- exp(-sum((x1 - x2)^2) / (rho)) # Euclidean Norm
    return(Kxy)
  }
  
  New_MatrixCal <- function(x, y, rho) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    xxt <- matrix(rep(diag(x %*% t(x)), dim(y)[1]),
                  dim(x)[1], dim(y)[1], byrow = FALSE)
    yyt <- matrix(rep(diag(y %*% t(y)), dim(x)[1]),
                  dim(x)[1], dim(y)[1], byrow = TRUE)
    Kmatsub0 <- exp(-((xxt + yyt - 2 * x %*% t(y)) / rho))
    return(Kmatsub0)
  }
  
  KernelEst <- function(x, y, lambda, rho, x0) {
    N <- dim(x)[1]
    Kmat <- New_MatrixCal(x, x, rho)
    tau <- 1 / lambda
    avec <- tau * solve(diag(rep(1, N)) + tau * Kmat) %*% y
    Kx0 <- New_MatrixCal(x0, x, rho)
    fx0 <- Kx0 %*% avec
    return(fx0)
  }
  
  ptm3 <- proc.time()
  ## Time of KMeans Center
  ## Proposed approach
  KMcluster <- stats::kmeans(x, n)
  centers <- as.matrix(KMcluster$centers)
  rhocenters <- det(stats::cov(centers))
  
  KCmat <- New_MatrixCal(centers, centers, rhocenters)
  KxcMat <- New_MatrixCal(x0s, centers, rhocenters)
  
  ## Time of Nearest Y
  # Nearest response y for centers
  ystar <- rep(0, n)
  distMat <- as.matrix(stats::dist(rbind(x, centers), diag = TRUE))[1:N, (N + 1):(n + N)]
  
  for (a in 1:n) {
    ystar[a] <- y[which.min(distMat[, a])]
  }
  
  # Time lambda0
  ptm_33 <- proc.time()
  
  eig <- eigen(KCmat)
  eigenvalues <- eig$values
  Q <- eig$vectors

  A1 <- t(KCmat) %*% Q
  A2 <- t(Q) %*% ystar
  trmat <- t(Q) %*% t(KCmat) %*% Q
  
  KM_BIC0 <- c(rep(0, length(lambda0)))
  for (i in 1:length(lambda0)) {
    diagmat <- diag(1 / (eigenvalues + lambda0[i]))
    fx <- A1 %*% diagmat %*% A2
    trv <- tr(diagmat %*% trmat)
    KM_BIC0[i] <- n * log(t(c(ystar - fx)) %*% c(ystar - fx)) + log(n) * trv
  }
  Lambda0 <- lambda0[which.min(KM_BIC0)]
  lambda_paths <- Lambda0
  TimeKMeans_L0s <- (proc.time() - ptm_33)
  
  # weights
  lambdaOri <- Lambda0
  tauOri <- 1 / lambdaOri
  Kx0s <- solve(diag(rep(1, n)) + tauOri * KCmat) %*% t(KxcMat) ### Here n*r matrix
  weights <- as.matrix(abs(Kx0s) / colSums(abs(Kx0s)))
  finalweights <- n * weights / N
  finalprobs <- finalweights[KMcluster$cluster, ]
  
  fx0s <- matrix(0, B, dim(x0s)[1])
  newfx0ss <- rep(0, dim(x0s)[1])
  
  ptm_36_0 <- proc.time()
  for (e in 1:(dim(x0s)[1])) {
    # Time Estimation
    for (b in 1:B) {
      XIND <- sample(c(1:N), n, prob = finalprobs[, e], replace = TRUE)
      xsub <- x[XIND, ]
      ysub <- y[XIND]
      rho_sub <- det(stats::cov(xsub))
      ###
      Kmatsub <- New_MatrixCal(xsub, xsub, rho_sub)
      
      avec <- tauOri * solve(diag(rep(1, n)) + tauOri * Kmatsub) %*% ysub
      
      Kx0 <- matrix(0, dim(x0s)[1], n)
      
      k0 <- apply(xsub, 1, GaussianKernel, x2 = x0s[e, ], rho = rho_sub)
      
      fx0s[b, e] <- t(k0) %*% avec
    }
    newfx0ss[e] <- mean(fx0s[, e])
  } # e
  
  TimeKMeans_Ests <- proc.time() - ptm_36_0
  TimeKMeans_Totals <- proc.time() - ptm3
  
  est <- matrix(newfx0ss, dim(x0s)[1], 1)
  
  MSEKMean <- apply(est, 2, function(x) (y0s - est)^2)
  Ave_MSE <- mean(MSEKMean)
  
  Time_Table <- rbind(TimeKMeans_Totals,
                      TimeKMeans_L0s,
                      TimeKMeans_Ests)[, 1:3]
  
  t_convt <- paste(Time_Table %/% 3600, ":",
                   (Time_Table %% 3600) %/% 60, "'",
                   round((Time_Table %% 3600) %% 60, 0), sep = "")
  Converted_Time_Table <- matrix(t_convt, 3, 3, byrow = FALSE)
  rownames(Converted_Time_Table) <- c("Total", "Tuning", "Est")
  colnames(Converted_Time_Table) <- c("user.self", "sys.self", "elapsed")
  
  return(list(
    coeffi = avec,
    y0s = y0s,
    est = est,
    MSEKMean = MSEKMean,
    Ave_MSE = Ave_MSE,
    optimal_lam = lambda_paths,
    Converted_Time_Table = Converted_Time_Table[, 1],
    x0s = x0s
  ))
}
