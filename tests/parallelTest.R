require(confintROB)
require(lme4)
require(robustlmm)
require(robustvarComp)
require(foreach)
require(parallel)
require(doParallel)

cl <- makeCluster(2, type = "PSOCK")
doParallel::registerDoParallel(cl)

control <- lmerControl(check.conv.grad = "ignore")

model.ML <-
  lmer(Yield ~ (1 | Batch),
       Dyestuff,
       REML = FALSE,
       control = control)

set.seed(123)
parallelResults <- confintROB(model.ML, nsim = 3)

stopCluster(cl)
registerDoSEQ()

set.seed(123)
singleThreadedResults <- confintROB(model.ML, nsim = 3)
stopifnot(all.equal(parallelResults, singleThreadedResults))