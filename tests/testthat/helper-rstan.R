## Test helper: detect whether it is safe to load rstan (or any rstan-
## dependent package, e.g. rstanarm) in the current R process.
##
## rstan and cmdstanr using model methods in the same R process may
## lead to segmentation faults. Both bundle the Stan C++ library and
## ship overlapping global symbols (stan::math::*, boost::*, Eigen
## lazy-init tables); Linux's dyn.load resolves the second backend's
## symbols into incompatible memory laid out for the first, and
## segfaults at load time. Aligning Stan major.minor versions does not
## always fix this because the two backends may have been compiled
## with different preprocessor flags and produce different template
## instantiations.
##
## Practical consequence in this test suite: any test file that runs
## AFTER a cmdstanr-using test in the same `devtools::test()`
## invocation cannot safely load rstan; even `requireNamespace("rstan")`
## or `skip_if_not_installed("rstan")` will segfault because they
## internally call library.dynam() on rstan.so.
##
## This helper checks rstan installation via the filesystem only
## (system.file()) and refuses to attempt the load if cmdstanr has
## been loaded into the process. Callers run the rstan branch in a
## callr subprocess if they really need the assertion.
skip_if_stan_unsafe <- function(pkg = "rstan") {
  if (!nzchar(system.file(package = pkg))) {
    testthat::skip(paste0(pkg, " not installed"))
  }
  if ("cmdstanr" %in% loadedNamespaces()) {
    testthat::skip(paste0(
      "cmdstanr already loaded in this R process; loading ", pkg,
      " here would segfault (rstan + cmdstanr in-process Stan C++ ",
      "symbol collision). Run this test in isolation via ",
      "testthat::test_file() or in a callr::r() subprocess."))
  }
  invisible(TRUE)
}
