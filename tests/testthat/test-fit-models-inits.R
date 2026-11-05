# Regression tests for user-defined initial values in MLE models

# Prepare example data
lung_test <- survival::lung
lung_test$event <- as.integer(lung_test$status == 2)

form <- survival::Surv(time, event) ~ 1

# Starting values
wei_inits <- c(1.3, 420)
exp_inits <- c(0.003)
lnorm_inits <- c(5.5, 0.9)
genf_inits <- c(5.543222, 1, 0, 1)

# Helper: check initial values and convergence
expect_valid_fit <- function(model, expected_inits = NULL) {
  expect_equal(model$opt$convergence, 0)

  if (is.null(expected_inits)) {
    expect_null(model$call$inits)
  } else {
    expect_equal(
      as.numeric(model$call$inits),
      expected_inits
    )
  }
}

test_that("single distribution accepts numeric and list inits", {
  m1 <- fit.models(
    form, lung_test, distr = "wei",
    method = "mle", inits = wei_inits
  )

  m2 <- fit.models(
    form, lung_test, distr = "wei",
    method = "mle", inits = list(wei_inits)
  )

  expect_valid_fit(m1$models[[1]], wei_inits)
  expect_valid_fit(m2$models[[1]], wei_inits)

  expect_equal(
    m1$models[[1]]$res[, "est"],
    m2$models[[1]]$res[, "est"],
    tolerance = 1e-4
  )
})

test_that("multiple distributions receive corresponding inits", {
  dists <- c("exp", "wei", "lnorm", "genf")

  inits <- list(
    exp_inits,
    wei_inits,
    lnorm_inits,
    genf_inits
  )

  m <- fit.models(
    form, lung_test, distr = dists,
    method = "mle", inits = inits
  )

  expect_length(m$models, 4)

  for (i in seq_along(dists)) {
    expect_valid_fit(m$models[[i]], inits[[i]])
  }
})

test_that("NULL entries use default initial values", {
  m <- fit.models(
    form, lung_test,
    distr = c("exp", "wei"),
    method = "mle",
    inits = list(NULL, wei_inits)
  )

  expect_valid_fit(m$models[[1]])
  expect_valid_fit(m$models[[2]], wei_inits)
})

test_that("GenF accepts a numeric initial-value vector", {
  m <- fit.models(
    form, lung_test,
    distr = "genf",
    method = "mle",
    inits = genf_inits
  )

  expect_valid_fit(m$models[[1]], genf_inits)
})

test_that("default fitting remains unchanged", {
  m_default <- fit.models(
    form, lung_test,
    distr = c("exp", "wei"),
    method = "mle"
  )

  m_custom <- fit.models(
    form, lung_test,
    distr = "wei",
    method = "mle",
    inits = wei_inits
  )

  expect_valid_fit(m_default$models[[1]])
  expect_valid_fit(m_default$models[[2]])

  expect_equal(
    m_default$models[[2]]$res[, "est"],
    m_custom$models[[1]]$res[, "est"],
    tolerance = 1e-4
  )
})

test_that("incorrect inits list length is rejected", {
  expect_error(
    fit.models(
      form, lung_test,
      distr = c("exp", "wei"),
      method = "mle",
      inits = list(wei_inits)
    ),
    "one element per"
  )
})



test_that("bc example supports multiple distributions", {

  data("bc", package = "flexsurv", envir = environment())
  f <- survival::Surv(recyrs, censrec) ~ group

  # Original single-distribution example
  m_single <- fit.models(
    formula = f,
    data = bc,
    distr = "wei",
    inits = c(2, 5),
    method = "mle"
  )

  expect_equal(
    as.numeric(m_single$models[[1]]$call$inits),
    c(2, 5)
  )
  expect_equal(m_single$models[[1]]$opt$convergence, 0)

  # Multiple distributions: default for exp,
  # user-defined starting values for Weibull
  m_multi <- fit.models(
    formula = f,
    data = bc,
    distr = c("exp", "wei"),
    inits = list(NULL, c(2, 5)),
    method = "mle"
  )

  expect_length(m_multi$models, 2)
  expect_null(m_multi$models[[1]]$call$inits)
  expect_equal(
    as.numeric(m_multi$models[[2]]$call$inits),
    c(2, 5)
  )

  expect_equal(m_multi$models[[1]]$opt$convergence, 0)
  expect_equal(m_multi$models[[2]]$opt$convergence, 0)
})


test_that("bc models receive separate initial-value vectors", {

  data("bc", package = "flexsurv", envir = environment())
  f <- survival::Surv(recyrs, censrec) ~ group

  # Obtain valid initial values for each model
  m_default <- fit.models(
    formula = f,
    data = bc,
    distr = c("exp", "wei"),
    method = "mle"
  )

  initial_values <- lapply(
    m_default$models,
    function(m) unname(m$res[, "est"])
  )

  # Refit both distributions using their own vectors
  m_custom <- fit.models(
    formula = f,
    data = bc,
    distr = c("exp", "wei"),
    inits = initial_values,
    method = "mle"
  )

  expect_length(m_custom$models, 2)

  for (i in seq_along(initial_values)) {
    expect_equal(
      as.numeric(m_custom$models[[i]]$call$inits),
      as.numeric(initial_values[[i]])
    )
    expect_equal(m_custom$models[[i]]$opt$convergence, 0)
  }
})
