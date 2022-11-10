
test_that("plots don't throw an error and return ggplot2 list", {
  
  data(bc)
  
  mle <- fit.models(formula=Surv(recyrs,censrec) ~ group, data=bc,
                    distr="exp", method="mle")
  
  p.mle <- make.surv(mle)
  p.mle2 <- make.surv(mle, nsim = 2)
  
  expect_type(psa.plot(p.mle), "list")
  expect_type(plot(mle, nsim=10), "list")
})
