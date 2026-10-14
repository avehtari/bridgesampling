# rstan (and other packages built against StanHeaders) and the cmdstanr
# model methods each contain a copy of Stan Math. On Linux, GCC emits the
# static data of stan::math::AutodiffStackSingleton as STB_GNU_UNIQUE
# symbols, which the dynamic loader binds once per process. So when one of
# the two is loaded, loading the other in the same R process segfaults in
# dyn.load(), in either order (stan-dev/cmdstanr#1177). These helpers skip a
# test that would load the second copy.

# The mangled name of AutodiffStackSingleton<vari_base, chainable_alloc>::instance_
.stan_math_symbol <- paste0(
  "_ZN4stan4math22AutodiffStackSingletonINS0_9vari_baseENS0_15",
  "chainable_allocEE9instance_E"
)

# Is a library loaded by Rcpp::sourceCpp() that contains Stan Math? That is
# the cmdstanr model methods (and $expose_functions()); their DLLs are named
# sourceCpp_<n>. rstan model DLLs (file<hex>) also contain Stan Math, but
# they come from the same Stan Math as rstan.so, so they are not counted.
# rstan.so itself registers its routines and turns dynamic lookup off.
.cmdstanr_methods_loaded <- function() {
  any(vapply(getLoadedDLLs(), function(dll) {
    startsWith(dll[["name"]], "sourceCpp_") && !is.null(tryCatch(
      getNativeSymbolInfo(.stan_math_symbol, PACKAGE = dll),
      error = function(e) NULL
    ))
  }, logical(1)))
}

# Is the DLL of a package loaded that is StanHeaders or is built against it
# (rstan, rstanarm, ...)?
.rstan_loaded <- function() {
  any(vapply(names(getLoadedDLLs()), function(name) {
    if (identical(name, "StanHeaders")) {
      return(TRUE)
    }
    linking_to <- tryCatch(
      utils::packageDescription(name, fields = "LinkingTo"),
      error = function(e) NA, warning = function(w) NA
    )
    !is.na(linking_to) && grepl("StanHeaders", linking_to, fixed = TRUE)
  }, logical(1)))
}

# Use before a test loads rstan, also before skip_if_not_installed("rstan"),
# which loads rstan with requireNamespace().
skip_if_cmdstanr_methods_loaded <- function() {
  if (.cmdstanr_methods_loaded()) {
    testthat::skip(paste(
      "cmdstanr model methods are loaded in this R process;",
      "loading rstan would segfault (stan-dev/cmdstanr#1177)"
    ))
  }
}

# Use before a test calls cmdstanr model methods ($init_model_methods(),
# $log_prob(), $grad_log_prob(), $unconstrain_draws()).
skip_if_rstan_loaded <- function() {
  if (.rstan_loaded()) {
    testthat::skip(paste(
      "rstan is loaded in this R process;",
      "cmdstanr model methods would segfault (stan-dev/cmdstanr#1177)"
    ))
  }
}
