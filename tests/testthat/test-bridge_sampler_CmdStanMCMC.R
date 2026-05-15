context("test bridge_sampler cmdstanmcmc method")

testthat::test_that("bridge_sampler() works for CmdStanMCMC and basic sanity checks", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not_installed("bridgesampling")
  testthat::skip_if_not_installed("posterior")

  # Require a working CmdStan toolchain
  if (!file.exists(cmdstanr::cmdstan_path())) {
    testthat::skip("CmdStan is not installed in the expected path for cmdstanr.")
  }

  set.seed(123)

  N     <- 60L
  sigma <- 1
  y     <- rnorm(N, mean = 0.5, sd = sigma)

  data_list <- list(N = N, y = y, sigma = sigma)

  stan_code <- "
  data {
    int<lower=1> N;
    vector[N] y;
    real<lower=0> sigma;
  }
  parameters {
    real mu;
  }
  model {
    mu ~ normal(0, 1);
    y ~ normal(mu, sigma);
  }"

  tf <- tempfile(fileext = ".stan")
  on.exit(unlink(tf), add = TRUE)
  writeLines(stan_code, tf)
  mod <- cmdstanr::cmdstan_model(tf, quiet = TRUE, force_recompile = TRUE)

  fit <- mod$sample(
    data = data_list,
    seed = 202,
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 2000,
    iter_sampling = 4000,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )

  bs <- bridgesampling::bridge_sampler(fit, silent = TRUE, use_neff = FALSE)

  testthat::expect_s3_class(fit, "CmdStanMCMC")
  testthat::expect_true(is.list(bs))
  testthat::expect_true(is.finite(bs$logml))

  fit2 <- mod$sample(
    data = data_list,
    seed = 203,
    chains = 4,
    parallel_chains = 2,
    iter_warmup = 750,
    iter_sampling = 2000,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
  bs2 <- bridgesampling::bridge_sampler(fit2, silent = TRUE, use_neff = FALSE)
  testthat::expect_true(is.finite(bs2$logml))
})

testthat::test_that("CmdStanMCMC bridge estimate roughly agrees with rstan", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not_installed("bridgesampling")
  testthat::skip_if_not_installed("callr")

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
  ## via callr::r(), then compare the two scalar bridge-sampling
  ## estimates here.
  ##
  ## IMPORTANT: do NOT use requireNamespace("rstan") or
  ## skip_if_not_installed("rstan") in this parent process. By the time
  ## this test_that block runs, the cmdstanr-using test_that block above
  ## has already loaded a Stan-symbol-carrying .so into the process via
  ## bridge_sampler.CmdStanMCMC(); attempting to load rstan into the
  ## same process segfaults at library.dynam(). Detect rstan via
  ## system.file() (pure filesystem lookup, no dyn.load).
  if (!nzchar(system.file(package = "rstan")))
    testthat::skip("rstan not installed")

  if (!file.exists(cmdstanr::cmdstan_path())) {
    testthat::skip("CmdStan is not installed in the expected path for cmdstanr.")
  }

  set.seed(456)

  # Same model/data as above
  N     <- 60L
  sigma <- 1
  y     <- rnorm(N, mean = 0.25, sd = sigma)
  data_list <- list(N = N, y = y, sigma = sigma)

  stan_code <- "
  data {
    int<lower=1> N;
    vector[N] y;
    real<lower=0> sigma;
  }
  parameters {
    real mu;
  }
  model {
    mu ~ normal(0, 1);
    y ~ normal(mu, sigma);
  }"

  # --- CmdStanR fit (current process) ---
  tf <- tempfile(fileext = ".stan")
  on.exit(unlink(tf), add = TRUE)
  writeLines(stan_code, tf)
  mod_cs <- cmdstanr::cmdstan_model(tf, quiet = TRUE)
  fit_cs <- mod_cs$sample(
    data = data_list,
    seed = 777,
    chains = 4,
    parallel_chains = 2,
    iter_warmup = 10000,
    iter_sampling = 3000,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
  bs_cmd <- bridgesampling::bridge_sampler(fit_cs, silent = TRUE, use_neff = FALSE)
  testthat::expect_true(is.finite(bs_cmd$logml))

  # --- RStan fit (separate subprocess to avoid in-process Stan symbol
  #     collision with cmdstanr's already-loaded compiled model) ---
  pkg_dir <- pkgload::pkg_path()
  rstan_result <- callr::r(
    function(pkg_dir, stan_code, data_list) {
      suppressPackageStartupMessages(devtools::load_all(pkg_dir, quiet = TRUE))
      sm <- rstan::stan_model(model_code = stan_code)
      fit_rs <- rstan::sampling(
        sm, data = data_list, seed = 777,
        chains = 4, iter = 10000, warmup = 3000,
        refresh = 0, cores = 2
      )
      bs_rstan <- bridge_sampler(fit_rs, silent = TRUE, use_neff = FALSE)
      list(logml = bs_rstan$logml)
    },
    args = list(pkg_dir = pkg_dir, stan_code = stan_code,
                data_list = data_list)
  )
  testthat::expect_true(is.finite(rstan_result$logml))

  # Compare the two bridge estimates. Tolerance accounts for MC/bridge variance.
  testthat::expect_lt(abs(bs_cmd$logml - rstan_result$logml), 0.5)
})
