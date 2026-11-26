test_that("make.surv keeps the order of the newdata profiles", {
  d <- data.frame(time=c(1,2,4,8,3,6,12,24), event=1,
                  arm=factor(rep(c("A","B"), each=4)))
  f <- fit.models(Surv(time,event)~arm, data=d, distr="exp")
  set.seed(1)
  s_ab <- make.surv(f, t=5, newdata=list(list(arm="A"), list(arm="B")))
  set.seed(1)
  s_ba <- make.surv(f, t=5, newdata=list(list(arm="B"), list(arm="A")))
  expect_equal(unname(s_ba$des.mat[,"armB"]), c(1,0))
  expect_equal(s_ba$des.mat, s_ab$des.mat[2:1,,drop=FALSE])
  expect_equal(s_ba$S, s_ab$S[2:1])
})
