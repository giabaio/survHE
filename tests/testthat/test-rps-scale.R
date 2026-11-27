test_that("make.surv uses the scale and timescale of a fitted rps model", {
  data(bc, package="flexsurv")
  nd <- list(list(group="Good"), list(group="Poor"))
  for (opt in list(list(scale="odds"), list(scale="normal"), list(timescale="identity"))) {
    f <- do.call(fit.models, c(list(formula=Surv(recyrs,censrec)~group, data=bc,
                                    distr="rps", k=1), opt))
    set.seed(1)
    s <- make.surv(f, t=c(2,5), newdata=nd)
    ref <- summary(f$models[[1]], t=c(2,5), ci=FALSE,
                   newdata=data.frame(group=c("Good","Poor")))
    expect_equal(c(s$S[[1]]$S, s$S[[2]]$S), c(ref[[1]]$est, ref[[2]]$est),
                 tolerance=0.01)
  }
})
