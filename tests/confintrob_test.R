require(confintROB)
require(lme4)
require(robustlmm)

test <- function(object, ...) {
  doTest <- function(method, ...) {
    cat(
      "Running test for object of class ",
      class(object),
      " with arguments method = ",
      method,
      "\n"
    )
    set.seed(1234)
    result <-
      confintROB(
        object = object,
        level = .95,
        method = method,
        nsim = 10,
        boot.type = "parametric",
        clusterID = "Subject",
        ...
      )
    print(result, digits = 2)
  }
  
  for (method in c("boot", "BCa", "Wald")) {
    doTest(method, ...)
  }
}

test.wild <- function(object, ...) {
  cat("Running test.wild for object of class ",  class(object), "\n")
  set.seed(123)
  y 		 <- getME(object, "y")
  X      <- as.matrix(getME(object, "X"))
  id     <- getME(object, "flist")[[1]]
  bet    <- unname(fixef(object))
  result <-
    confintROB:::createWildSampleFunction(y = y,
                                          X = X,
                                          id = id,
                                          bet = bet)(1)
  print(result, digits = 5)
}


control <- lmerControl(check.conv.grad = "ignore")

model.ds.ML <-
  lmer(Yield ~ (1 | Batch),
       Dyestuff,
       REML = FALSE,
       control = control)
print(summary(model.ds.ML), digits = 2)
test(model.ds.ML)
test.wild(model.ds.ML, .export = "control")

model.ds.DAStau <-
  rlmer(
    Yield ~ (1 | Batch),
    Dyestuff,
    rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
    rho.b = chgDefaults(smoothPsi, k = 5.14, s = 10),
    rho.sigma.b = chgDefaults(smoothPsi, k = 5.14, s = 10),
    init = function(...)
      lmer(..., control = control)
  )
print(summary(model.ds.DAStau), digits = 2)
test(model.ds.DAStau)
test.wild(model.ds.DAStau, .export = "control")

model.ss.ML <-
  lmer(
    Reaction ~ Days + (Days | Subject),
    data = sleepstudy,
    REML = FALSE,
    control = control
  )
print(summary(model.ss.ML), digits = 2)
test.wild(model.ss.ML, .export = "control")

model.ss.DAStau <-
  rlmer(
    Reaction ~ Days + (Days | Subject),
    data = sleepstudy,
    rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
    rho.b = chgDefaults(smoothPsi, k = 5.14, s = 10),
    rho.sigma.b = chgDefaults(smoothPsi, k = 5.14, s = 10),
    init = function(...)
      lmer(..., control = control)
  )
print(summary(model.ss.DAStau), digits = 2)
test.wild(model.ss.DAStau, control = control, .export = "control")