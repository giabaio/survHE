test_that("MLE profiles retain shared Weibull shape draws", {
  d <- data.frame(time=c(1,2,4,8,1.2,2.4,4.8,9.6), event=1,
                  arm=factor(rep(c("A","B"), each=4)))
  f <- fit.models(Surv(time,event)~arm, data=d, distr="weibull")
  for (nd in list(NULL, list(list(arm="A"), list(arm="B")))) {
    set.seed(1)
    s <- suppressMessages(make.surv(f, nsim=20, t=c(0,2), newdata=nd))
    expect_identical(s$sim[[1]][,"shape"], s$sim[[2]][,"shape"])
  }
})
