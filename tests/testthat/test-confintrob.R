# Tests for lmerMod and rlmerMod objects
# These tests require robustlmm (in Suggests)

test_that("confintROB works with lmerMod objects", {
  skip_if_not_installed("robustlmm")
  library(lme4)
  library(robustlmm)

  # Use control variable to test that bakeCall() properly handles this
  control <- lmerControl(check.conv.grad = "ignore")

  model.ds.ML <-
    lmer(Yield ~ (1 | Batch),
         Dyestuff,
         REML = FALSE,
         control = control)

  expect_snapshot({
    print(summary(model.ds.ML), digits = 2)
  })

  # Test boot method
  expect_snapshot({
    cat("Running test for object of class lmerMod with method = boot\n")
    set.seed(1234)
    result <- confintROB(
      object = model.ds.ML,
      level = .95,
      method = "boot",
      nsim = 10,
      boot.type = "parametric"
    )
    print(result, digits = 2)
  })

  # Test BCa method - this tests that bakeCall() works correctly
  expect_snapshot({
    cat("Running test for object of class lmerMod with method = BCa\n")
    set.seed(1234)
    result <- confintROB(
      object = model.ds.ML,
      level = .95,
      method = "BCa",
      nsim = 10,
      boot.type = "parametric"
    )
    print(result, digits = 2)
  })

  # Test Wald method
  expect_snapshot({
    cat("Running test for object of class lmerMod with method = Wald\n")
    set.seed(1234)
    result <- confintROB(
      object = model.ds.ML,
      level = .95,
      method = "Wald",
      nsim = 10,
      boot.type = "parametric"
    )
    print(result, digits = 2)
  })

  # Test wild bootstrap sampling
  expect_snapshot({
    cat("Running test.wild for object of class lmerMod\n")
    set.seed(123)
    y <- getME(model.ds.ML, "y")
    X <- as.matrix(getME(model.ds.ML, "X"))
    id <- getME(model.ds.ML, "flist")[[1]]
    bet <- unname(fixef(model.ds.ML))
    result <- confintROB:::createWildSampleFunction(
      y = y, X = X, id = id, bet = bet
    )(1)
    print(result, digits = 5)
  })
})

test_that("confintROB works with rlmerMod objects (Dyestuff)", {
  skip_if_not_installed("robustlmm")
  library(lme4)
  library(robustlmm)

  control <- lmerControl(check.conv.grad = "ignore")

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

  expect_snapshot({
    print(summary(model.ds.DAStau), digits = 2)
  })

  # Test boot method
  expect_snapshot({
    cat("Running test for object of class rlmerMod with method = boot\n")
    set.seed(1234)
    result <- confintROB(
      object = model.ds.DAStau,
      level = .95,
      method = "boot",
      nsim = 10,
      boot.type = "parametric",
      clusterID = "Subject"
    )
    print(result, digits = 2)
  })

  # Test BCa method
  expect_snapshot({
    cat("Running test for object of class rlmerMod with method = BCa\n")
    set.seed(1234)
    result <- confintROB(
      object = model.ds.DAStau,
      level = .95,
      method = "BCa",
      nsim = 10,
      boot.type = "parametric",
      clusterID = "Subject"
    )
    print(result, digits = 2)
  })

  # Test Wald method
  expect_snapshot({
    cat("Running test for object of class rlmerMod with method = Wald\n")
    set.seed(1234)
    result <- confintROB(
      object = model.ds.DAStau,
      level = .95,
      method = "Wald",
      nsim = 10,
      boot.type = "parametric",
      clusterID = "Subject"
    )
    print(result, digits = 2)
  })

  # Test wild bootstrap sampling
  expect_snapshot({
    cat("Running test.wild for object of class rlmerMod\n")
    set.seed(123)
    y <- getME(model.ds.DAStau, "y")
    X <- as.matrix(getME(model.ds.DAStau, "X"))
    id <- getME(model.ds.DAStau, "flist")[[1]]
    bet <- unname(fixef(model.ds.DAStau))
    result <- confintROB:::createWildSampleFunction(
      y = y, X = X, id = id, bet = bet
    )(1)
    print(result, digits = 5)
  })
})

test_that("wild bootstrap works with sleepstudy models", {
  skip_if_not_installed("robustlmm")
  library(lme4)
  library(robustlmm)

  control <- lmerControl(check.conv.grad = "ignore")

  model.ss.ML <-
    lmer(
      Reaction ~ Days + (Days | Subject),
      data = sleepstudy,
      REML = FALSE,
      control = control
    )

  expect_snapshot({
    print(summary(model.ss.ML), digits = 2)
  })

  expect_snapshot({
    cat("Running test.wild for object of class lmerMod (sleepstudy)\n")
    set.seed(123)
    y <- getME(model.ss.ML, "y")
    X <- as.matrix(getME(model.ss.ML, "X"))
    id <- getME(model.ss.ML, "flist")[[1]]
    bet <- unname(fixef(model.ss.ML))
    result <- confintROB:::createWildSampleFunction(
      y = y, X = X, id = id, bet = bet
    )(1)
    print(result, digits = 5)
  })

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

  expect_snapshot({
    print(summary(model.ss.DAStau), digits = 2)
  })

  expect_snapshot({
    cat("Running test.wild for object of class rlmerMod (sleepstudy)\n")
    set.seed(123)
    y <- getME(model.ss.DAStau, "y")
    X <- as.matrix(getME(model.ss.DAStau, "X"))
    id <- getME(model.ss.DAStau, "flist")[[1]]
    bet <- unname(fixef(model.ss.DAStau))
    result <- confintROB:::createWildSampleFunction(
      y = y, X = X, id = id, bet = bet
    )(1)
    print(result, digits = 5)
  })
})
