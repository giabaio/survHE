
test_that("output plots don't throw an error and return ggplot2 list", {
  
  data(bc)
  
  mle <- fit.models(formula=Surv(recyrs,censrec) ~ group, data=bc,
                    distr=c("exponential", "weibull"), method="mle")
  
  p.mle <- make.surv(mle)
  p.mle2 <- make.surv(mle, nsim = 2)
  p.mle3 <- make.surv(mle, nsim = 2, mod = 2)
  
  expect_type(plot(mle, nsim=10), "list")
  
  expect_type(psa.plot(p.mle), "list")
  expect_type(psa.plot(p.mle), "list")
  
  expect_type(psa.plot(p.mle2), "list")
  
  expect_type(psa.plot(p.mle3), "list")
  
  expect_type(plot(mle, add.km=TRUE), "list")
  expect_type(plot(mle, add.km=TRUE, sim = 10), "list")
})
