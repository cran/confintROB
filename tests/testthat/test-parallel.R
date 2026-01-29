# Tests for parallel processing
# These tests require parallel and doParallel (in Suggests)

test_that("parallel and sequential results are identical", {
  skip_if_not_installed("parallel")
  skip_if_not_installed("doParallel")
  library(lme4)
  library(foreach)
  library(parallel)
  library(doParallel)

  control <- lmerControl(check.conv.grad = "ignore")

  model.ML <-
    lmer(Yield ~ (1 | Batch),
         Dyestuff,
         REML = FALSE,
         control = control)

  # Run in parallel
  cl <- makeCluster(2, type = "PSOCK")
  doParallel::registerDoParallel(cl)

  set.seed(123)
  parallelResults <- confintROB(model.ML, nsim = 3)

  stopCluster(cl)
  registerDoSEQ()

  # Run sequentially
  set.seed(123)
  singleThreadedResults <- confintROB(model.ML, nsim = 3)

  expect_equal(parallelResults, singleThreadedResults)
})
