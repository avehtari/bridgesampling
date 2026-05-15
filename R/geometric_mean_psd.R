## Matrix geometric mean of two positive (semi-)definite covariance
## matrices: A # B = A^{1/2} (A^{-1/2} B A^{-1/2})^{1/2} A^{1/2}.
##
## Reference: Seyboldt, Carlson, Carpenter (2026), "Preconditioning
## Hamiltonian Monte Carlo by minimizing Fisher Divergence". The
## dense Fisher-divergence-optimal preconditioner is the matrix
## geometric mean of the sample covariance of the draws and the
## inverse covariance of the scores; the helper here is the
## bridge-sampling analogue used for `proposal_fit = "hybrid"`.
##
## The geometric mean is the (affine-invariant Fisher-Rao)
## Riemannian midpoint of A and B on the positive-definite cone.
##
## Implementation: symmetric eigendecomposition of A, then of the
## inner sandwich M = A^{-1/2} B A^{-1/2}, with PSD clipping on tiny
## (round-off) negative eigenvalues and a pseudo-inverse guard
## (1/sqrt(lambda) -> 0) for exactly-zero eigenvalues.
##
## Inputs are assumed positive semi-definite and already
## nearPD-corrected by the caller (sample/score covariance
## estimates always are). On any un-recoverable input
## (eigendecomposition failure, non-finite eigenvalues, non-finite
## output) the function errors.
##
## Arguments:
##   A, B: p x p symmetric matrices, assumed positive semi-definite.
##         A and B should already be nearPD-corrected by the caller
##         when they are sample/score covariance estimates.
##
## Returns: a p x p symmetric matrix.
.geometric_mean_psd <- function(A, B) {
  if (!is.matrix(A) || !is.matrix(B))
    stop(".geometric_mean_psd(): A and B must be matrices.", call. = FALSE)
  if (!identical(dim(A), dim(B)) || nrow(A) != ncol(A))
    stop(".geometric_mean_psd(): A and B must be square and the same size.",
         call. = FALSE)

  ## Force exact symmetry up front -- both eigen() and the
  ## sandwiches below assume it.
  A <- (A + t(A)) / 2
  B <- (B + t(B)) / 2

  eigA <- tryCatch(eigen(A, symmetric = TRUE),
                   error = function(e) NULL)
  if (is.null(eigA) || any(!is.finite(eigA$values)))
    stop(".geometric_mean_psd(): eigendecomposition of A failed ",
         "(non-finite eigenvalues or eigen() error).", call. = FALSE)

  ## Clip tiny negative eigenvalues from round-off. We do NOT bail
  ## out on small positive eigenvalues: the caller pre-regularises
  ## A via Matrix::nearPD(), so A is PD by construction; small
  ## eigenvalues are real and the geometric mean is mathematically
  ## well-defined via the pseudo-inverse guard on invsL below.
  lam <- pmax(eigA$values, 0)
  U     <- eigA$vectors
  sqrtL <- sqrt(lam)
  invsL <- ifelse(sqrtL > 0, 1 / sqrtL, 0)    # pseudo-inverse for true zeros
  A_half     <- U %*% (sqrtL * t(U))
  A_neg_half <- U %*% (invsL * t(U))

  M <- A_neg_half %*% B %*% A_neg_half
  M <- (M + t(M)) / 2

  eigM <- tryCatch(eigen(M, symmetric = TRUE),
                   error = function(e) NULL)
  if (is.null(eigM) || any(!is.finite(eigM$values)))
    stop(".geometric_mean_psd(): eigendecomposition of ",
         "A^{-1/2} B A^{-1/2} failed.", call. = FALSE)

  mu <- pmax(eigM$values, 0)
  V    <- eigM$vectors
  M_half <- V %*% (sqrt(mu) * t(V))

  G <- A_half %*% M_half %*% A_half
  G <- (G + t(G)) / 2

  if (any(!is.finite(G)))
    stop(".geometric_mean_psd(): non-finite entries in A # B; ",
         "caller's nearPD preprocessing should make this unreachable.",
         call. = FALSE)

  G
}
