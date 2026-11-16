test_that("state occupancy restarts each simulation", {
  p <- function(x) tibble::tibble(
    profile = 1, time = 0:1, lambda1 = x, lambda2 = x
  )
  x <- make_state_occupancy(
    2, p(c(1, .5)), p(c(0, 0)), p(c(0, .5)),
    p(c(1, 1)), p(c(0, 0)), c(100, 0, 0)
  )
  # Each draw starts at 100; one cycle moves half to death.
  expect_equal(x$`Pre-progressed`, c(100, 50, 100, 50))
  expect_equal(x$Progressed, rep(0, 4))
  expect_equal(x$Death, c(0, 50, 0, 50))
})
