context("Hybrid score-matched proposal covariance (matrix geometric mean)")

# Helper: simulate iid samples + analytic gradients from N(mu, Sigma).
.sim_gaussian <- function(n, mu, Sigma, seed = 1L) {
  set.seed(seed)
  L <- chol(Sigma)
  Z <- matrix(rnorm(n * length(mu)), nrow = n)
  X <- Z %*% L + matrix(mu, nrow = n, ncol = length(mu), byrow = TRUE)
  G <- -t(solve(Sigma, t(X) - mu))
  list(samples = X, gradients = G)
}

# -----------------------------------------------------------------------
# 1. .geometric_mean_psd(): matrix geometric mean invariants
# -----------------------------------------------------------------------
test_that(".geometric_mean_psd() satisfies Sigma # Sigma == Sigma", {
  set.seed(1)
  for (p in c(1L, 2L, 5L, 20L)) {
    S <- crossprod(matrix(rnorm(p * p), p)) + diag(p)
    G <- .geometric_mean_psd(S, S)
    expect_equal(G, S, tolerance = 1e-9,
                 info = sprintf("p = %d", p))
  }
})

test_that(".geometric_mean_psd() is AIRM-symmetric: A # B == B # A", {
  set.seed(2)
  for (p in c(2L, 5L, 10L)) {
    A <- crossprod(matrix(rnorm(p * p), p)) + diag(p)
    B <- crossprod(matrix(rnorm(p * p), p)) + diag(p)
    GAB <- .geometric_mean_psd(A, B)
    GBA <- .geometric_mean_psd(B, A)
    expect_equal(GAB, GBA, tolerance = 1e-9,
                 info = sprintf("p = %d", p))
  }
})

test_that(".geometric_mean_psd() agrees with sqrt(Sigma) when A = I", {
  set.seed(3)
  p <- 4L
  S <- crossprod(matrix(rnorm(p * p), p)) + diag(p)
  es <- eigen(S, symmetric = TRUE)
  S_sqrt <- es$vectors %*% diag(sqrt(es$values)) %*% t(es$vectors)
  S_sqrt <- (S_sqrt + t(S_sqrt)) / 2
  G <- .geometric_mean_psd(diag(p), S)
  expect_equal(G, S_sqrt, tolerance = 1e-9)
})

test_that(".geometric_mean_psd() handles rank-deficient A via pseudo-inverse", {
  ## The helper does not bail on rank-deficient inputs: an exactly-zero
  ## eigenvalue is handled by the sqrt-pseudo-inverse guard
  ## (1/sqrt(lambda) -> 0) so the geometric mean is still defined
  ## (zero-eigenvalue directions are preserved as zero in the output,
  ## which the caller's downstream nearPD() then handles).
  A <- matrix(0, 2, 2); A[1, 1] <- 1
  out <- .geometric_mean_psd(A, diag(2))
  expect_true(all(is.finite(out)))
  expect_equal(dim(out), c(2L, 2L))
})

test_that(".geometric_mean_psd() errors on non-finite A", {
  ## Un-recoverable input: NaN in A makes the eigendecomposition
  ## return non-finite eigenvalues, which the helper rejects with an
  ## explicit error.
  A_bad <- matrix(c(1, NaN, NaN, 1), 2, 2)
  expect_error(.geometric_mean_psd(A_bad, diag(2)),
               "eigendecomposition of A failed")
})

# -----------------------------------------------------------------------
# 2. .score_covariance() recovers the sample covariance on Gaussian draws.
# -----------------------------------------------------------------------
test_that(".score_covariance() on Gaussian draws matches sample covariance", {
  ps <- c(1L, 5L, 20L)
  for (p in ps) {
    Sigma <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
    sim <- .sim_gaussian(n = 5000L, mu = rnorm(p), Sigma = Sigma, seed = 42L)
    sc <- .score_covariance(sim$gradients)
    expect_true(isTRUE(sc$ok))
    fro <- function(A) sqrt(sum(A^2))
    rel_err <- fro(sc$Sigma_score - Sigma) / fro(Sigma)
    expect_lt(rel_err, 0.1)
  }
})

# -----------------------------------------------------------------------
# 3. bridge_sampler.matrix(proposal_fit = "hybrid") on multivariate-Gaussian
#    targets matches the analytic logml = 0. The matrix geometric mean of
#    V_sample and Sigma_score (both consistent estimators of Sigma) is
#    again a consistent estimator of Sigma, so the proposal density is
#    near-optimal and the bridge estimate is within a few MCSE of zero.
# -----------------------------------------------------------------------
test_that("bridge_sampler.matrix hybrid matches analytic logml on Gaussian", {
  skip_if_not_installed("mvtnorm")
  for (p in c(1L, 5L, 20L)) {
    Sigma <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
    mu <- rep(0, p)
    sim <- .sim_gaussian(n = 4000L, mu = mu, Sigma = Sigma, seed = 7L)
    parameters <- paste0("x", seq_len(p))
    colnames(sim$samples) <- parameters
    colnames(sim$gradients) <- parameters
    lb <- setNames(rep(-Inf, p), parameters)
    ub <- setNames(rep(Inf, p), parameters)
    log_post <- function(theta, data) {
      mvtnorm::dmvnorm(theta, mean = data$mu, sigma = data$Sigma, log = TRUE)
    }
    bs <- bridge_sampler(
      samples = sim$samples,
      log_posterior = log_post,
      data = list(mu = mu, Sigma = Sigma),
      lb = lb,
      ub = ub,
      gradients = sim$gradients,
      proposal_fit = "hybrid",
      silent = TRUE
    )
    expect_true(is.finite(bs$logml))
    expect_lt(abs(bs$logml - 0), 5 * bs$mcse_logml + 0.1)
  }
})

# -----------------------------------------------------------------------
# 4. bridge_sampler.CmdStanMCMC end-to-end smoke test.
# -----------------------------------------------------------------------
test_that("bridge_sampler.CmdStanMCMC accepts proposal_fit = 'hybrid'", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  if (!file.exists(cmdstanr::cmdstan_path())) {
    skip("CmdStan is not installed in the expected path for cmdstanr.")
  }

  stan_code <- "
  parameters { real mu; }
  model { mu ~ normal(0, 1); }
  "
  tmp_stan <- tempfile(fileext = ".stan")
  writeLines(stan_code, tmp_stan)
  mod <- cmdstanr::cmdstan_model(tmp_stan, compile_model_methods = TRUE)
  fit <- mod$sample(
    iter_warmup = 200,
    iter_sampling = 400,
    chains = 2,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
  fit$init_model_methods()
  bs <- bridge_sampler(fit, proposal_fit = "hybrid", silent = TRUE)
  expect_true(is.finite(bs$logml))
})

# -----------------------------------------------------------------------
# 5. bridge_sampler.stanfit mirror.
# -----------------------------------------------------------------------
test_that("bridge_sampler.stanfit accepts proposal_fit = 'hybrid'", {
  skip_on_cran()
  skip_if_not_installed("callr")

  ## rstan and cmdstanr using model methods in the same R process may
  ## lead to segmentation faults. Both bundle the Stan C++ library and
  ## ship overlapping global symbols (stan::math::*, boost::*, Eigen
  ## lazy-init tables); Linux's dyn.load resolves the second backend's
  ## symbols into incompatible memory laid out for the first, and
  ## segfaults at load time. Aligning Stan major.minor versions does
  ## not always fix this because the two backends may have been
  ## compiled with different preprocessor flags and produce different
  ## template instantiations. We therefore run the cmdstanr branch in
  ## the current process and the rstan branch in a fresh subprocess
  ## via callr::r().
  ##
  ## IMPORTANT: do NOT use skip_if_not_installed("rstan") or
  ## requireNamespace("rstan") in this parent process. By the time this
  ## test_that block runs, the CmdStanMCMC-hybrid test_that block above
  ## has already loaded a Stan-symbol-carrying .so into the process via
  ## fit$init_model_methods(); attempting to load rstan into the same
  ## process segfaults at library.dynam(). Detect rstan via
  ## system.file() (pure filesystem lookup, no dyn.load).
  if (!nzchar(system.file(package = "rstan")))
    skip("rstan not installed")
  pkg_dir <- pkgload::pkg_path()
  result <- callr::r(
    function(pkg_dir) {
      suppressPackageStartupMessages(devtools::load_all(pkg_dir, quiet = TRUE))
      stan_code <- "
        parameters { real mu; }
        model { mu ~ normal(0, 1); }
      "
      fit <- tryCatch(
        rstan::stan(
          model_code = stan_code,
          iter       = 600,
          warmup     = 200,
          chains     = 2,
          refresh    = 0
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) return(list(skip = TRUE, logml = NA_real_))
      bs <- bridge_sampler(fit, proposal_fit = "hybrid", silent = TRUE)
      list(skip = FALSE, logml = bs$logml)
    },
    args = list(pkg_dir = pkg_dir)
  )
  if (isTRUE(result$skip)) skip("rstan compile/sample unavailable in subprocess")
  expect_true(is.finite(result$logml))
})

# -----------------------------------------------------------------------
# 6. NaN gradient rows -> warning + fallback to sample covariance.
# -----------------------------------------------------------------------
test_that("NaN gradient rows trigger fallback warning", {
  skip_if_not_installed("mvtnorm")
  p <- 3L
  Sigma <- diag(p)
  sim <- .sim_gaussian(n = 1000L, mu = rep(0, p), Sigma = Sigma, seed = 23L)
  parameters <- paste0("x", seq_len(p))
  colnames(sim$samples) <- parameters
  colnames(sim$gradients) <- parameters
  # Poison MOST gradient rows so n_finite < p + 1 (forces fallback in helper).
  G_bad <- sim$gradients
  G_bad[1:998, ] <- NaN
  lb <- setNames(rep(-Inf, p), parameters)
  ub <- setNames(rep(Inf, p), parameters)
  log_post <- function(theta, data) {
    mvtnorm::dmvnorm(theta, mean = rep(0, p), sigma = diag(p), log = TRUE)
  }
  expect_warning(
    bridge_sampler(
      samples = sim$samples,
      log_posterior = log_post,
      data = NULL,
      lb = lb,
      ub = ub,
      gradients = G_bad,
      proposal_fit = "hybrid",
      silent = TRUE
    ),
    "score covariance unavailable",
    perl = FALSE
  )
})

# -----------------------------------------------------------------------
# 7. Bounded-parameter guard fires for the matrix dispatcher.
# -----------------------------------------------------------------------
test_that("Bounded parameters trigger fallback warning under hybrid", {
  skip_if_not_installed("mvtnorm")
  p <- 2L
  set.seed(31)
  X <- matrix(rexp(p * 1000), ncol = p) # strictly positive samples
  G <- -X # arbitrary, never used because of bounded-fallback
  parameters <- paste0("x", seq_len(p))
  colnames(X) <- parameters
  colnames(G) <- parameters
  lb <- setNames(rep(0, p), parameters) # bounded below at 0
  ub <- setNames(rep(Inf, p), parameters)
  log_post <- function(theta, data) {
    sum(dexp(theta, rate = 1, log = TRUE))
  }
  expect_warning(
    bridge_sampler(
      samples = X,
      log_posterior = log_post,
      data = NULL,
      lb = lb,
      ub = ub,
      gradients = G,
      proposal_fit = "hybrid",
      silent = TRUE
    ),
    "bounded parameters",
    perl = FALSE
  )
})

# -----------------------------------------------------------------------
# 8. warp3 silently ignores proposal_fit = "hybrid".
# -----------------------------------------------------------------------
test_that("method = 'warp3' ignores proposal_fit = 'hybrid'", {
  skip_if_not_installed("mvtnorm")
  p <- 3L
  Sigma <- diag(p)
  sim <- .sim_gaussian(n = 1500L, mu = rep(0, p), Sigma = Sigma, seed = 41L)
  parameters <- paste0("x", seq_len(p))
  colnames(sim$samples) <- parameters
  colnames(sim$gradients) <- parameters
  lb <- setNames(rep(-Inf, p), parameters)
  ub <- setNames(rep(Inf, p), parameters)
  log_post <- function(theta, data) {
    mvtnorm::dmvnorm(theta, mean = rep(0, p), sigma = diag(p), log = TRUE)
  }
  bs <- bridge_sampler(
    samples = sim$samples,
    log_posterior = log_post,
    data = NULL,
    lb = lb,
    ub = ub,
    method = "warp3",
    gradients = sim$gradients,
    proposal_fit = "hybrid",
    silent = TRUE
  )
  expect_true(is.finite(bs$logml))
  expect_equal(bs$method, "warp3")
})
