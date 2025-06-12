
test_that("output plots don't throw an error and return ggplot2 list", {
  
  data(bc)
  
  mle <- fit.models(formula=Surv(recyrs,censrec) ~ group, data=bc,
                    distr=c("exponential", "weibull"), method="mle")
  
  p.mle <- make.surv(mle)
  p.mle2 <- make.surv(mle, nsim = 2)
  p.mle3 <- make.surv(mle, nsim = 2, mod = 2)
  
  ## Need to change the test here from using expect_type(...,"list"), which 
  ## now breaks given the move to s7 for the newer version of 'ggplot2'
  expect_true(is_ggplot(plot(mle, nsim=10)))
  
  expect_true(is_ggplot(psa.plot(p.mle)))
  expect_true(is_ggplot(psa.plot(p.mle)))
  
  expect_true(is_ggplot(psa.plot(p.mle2)))
  
  expect_true(is_ggplot(psa.plot(p.mle3)))
  
  expect_true(is_ggplot(plot(mle, add.km=TRUE)))
  expect_true(is_ggplot(plot(mle, add.km=TRUE, sim = 10)))
})
