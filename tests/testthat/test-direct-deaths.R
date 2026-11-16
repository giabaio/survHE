test_that("direct deaths use the progression indicator", {
  d <- data.frame(id = 1:3, prog = c(0, 0, 1), death = c(1, 0, 1),
                  prog_t = c(5, 4, 2), death_t = c(5, 4, 6), treat = 0)
  x <- make_data_multi_state(d)
  expect_equal(x$status[x$id == 1 & x$trans == 2], 1)
  expect_equal(x$status[x$trans == 1], d$prog)
  expect_equal(x$status[x$id == 2], c(0, 0))
  expect_equal(x$status[x$id == 3 & x$trans %in% c(2, 3)], c(0, 1))
  expect_equal(x$Tstop[x$trans == 2], c(5, 4, 2))
})
