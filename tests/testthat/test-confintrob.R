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

  expect_snapshot(transform = strip_trailing_ws, {
    print(summary(model.ds.ML), digits = 2)
  })

  # Test boot method
  expect_snapshot(transform = strip_trailing_ws, {
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
  expect_snapshot(transform = strip_trailing_ws, {
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
  expect_snapshot(transform = strip_trailing_ws, {
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
  expect_snapshot(transform = strip_trailing_ws, {
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

  expect_snapshot(transform = strip_trailing_ws, {
    print(summary(model.ds.DAStau, df = "none"), digits = 2)
  })

  # Test boot method
  expect_snapshot(transform = strip_trailing_ws, {
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
  expect_snapshot(transform = strip_trailing_ws, {
    cat("Running test for object of class rlmerMod with method = BCa\n")
    set.seed(1234)
    result <- confintROB(
      object = model.ds.DAStau,
      level = .95,
      method = "BCa",
      nsim = 10,
      boot.type = "parametric",
      clusterID = "Batch"
    )
    print(result, digits = 2)
  })

  # Test Wald method
  expect_snapshot(transform = strip_trailing_ws, {
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
  expect_snapshot(transform = strip_trailing_ws, {
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

  expect_snapshot(transform = strip_trailing_ws, {
    print(summary(model.ss.ML), digits = 2)
  })

  expect_snapshot(transform = strip_trailing_ws, {
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

  expect_snapshot(transform = strip_trailing_ws, {
    print(summary(model.ss.DAStau, df = "none"), digits = 2)
  })

  expect_snapshot(transform = strip_trailing_ws, {
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

test_that("parametric bootstrap redraws the random effects", {
  library(lme4)

  model.ss.ML <-
    lmer(Reaction ~ Days + (Days | Subject),
         data = sleepstudy,
         REML = FALSE)

  set.seed(1234)
  result <- suppressMessages(confintROB(
    object = model.ss.ML,
    level = .95,
    method = "boot",
    nsim = 200,
    boot.type = "parametric"
  ))

  # If the random effects were held fixed, the bootstrap standard deviation of
  # the fixed effects would be less than half the model-based standard error.
  estimates <- attr(result, "fullResults")$bootstrap_estimates
  bootstrap.sd <- apply(estimates[, 1:2], 2, sd)
  model.se <- sqrt(diag(as.matrix(vcov(model.ss.ML))))
  ratio <- unname(bootstrap.sd / model.se)
  expect_gt(min(ratio), 0.7)
  expect_lt(max(ratio), 1.4)
})

test_that("wild bootstrap and jackknife use independent clusters", {
  library(lme4)

  # Nested: the clusters are the 10 batches, not the 30 casks
  model.nested <- lmer(strength ~ (1 | batch/cask), data = Pastes)
  expect_identical(names(getME(model.nested, "flist"))[1], "cask:batch")
  expect_identical(confintROB:::getClusters(model.nested), Pastes$batch)
  expect_silent(confintROB:::checkClusterID(model.nested, method = "boot",
                                            boot.type = "wild"))

  # clusterID overrides the default
  expect_identical(confintROB:::getClusters(model.nested, "cask:batch"),
                   getME(model.nested, "flist")[["cask:batch"]])
  expect_identical(confintROB:::getClusters(model.nested, Pastes$sample),
                   factor(Pastes$sample))
  expect_error(confintROB:::getClusters(model.nested, "Subject"),
               "not a grouping factor")
  expect_error(confintROB:::getClusters(model.nested, 1:3),
               "one value per observation")

  # The wild bootstrap gives the same weight to all casks of a batch
  y <- getME(model.nested, "y")
  X <- as.matrix(getME(model.nested, "X"))
  bet <- unname(fixef(model.nested))
  set.seed(1)
  ystar <- confintROB:::createWildSampleFunction(y, X, Pastes$batch, bet)(1)[[1]]
  ratio <- (ystar - drop(X %*% bet)) / (y - drop(X %*% bet))
  expect_true(all(tapply(ratio, Pastes$batch, function(r) diff(range(r))) <
                    1e-8))

  set.seed(1)
  result <- suppressWarnings(suppressMessages(
    confintROB(model.nested, method = "BCa", nsim = 5, boot.type = "wild")
  ))
  expect_named(attr(result, "fullResults"),
               c("BCa", "Percentile", "bootstrap_estimates", "biasBCa", "acc"))

  # Crossed: no independent clusters, clusterID is required
  model.crossed <- lmer(diameter ~ (1 | plate) + (1 | sample), data = Penicillin)
  expect_error(confintROB(model.crossed, nsim = 2, boot.type = "wild"),
               "crossed")
  expect_error(confintROB(model.crossed, method = "BCa", nsim = 2,
                          boot.type = "parametric"),
               "crossed")
  expect_message(suppressWarnings(
    confintROB(model.crossed, nsim = 2, boot.type = "parametric",
               clusterID = "plate")
  ), "Ignoring argument 'clusterID'")
  set.seed(1)
  result <- suppressWarnings(suppressMessages(
    confintROB(model.crossed, nsim = 2, boot.type = "wild", clusterID = "plate")
  ))
  expect_identical(dim(result), c(4L, 2L))
})

test_that("rlmer nested fits use the coarsest grouping factor", {
  skip_if_not_installed("robustlmm")
  library(robustlmm)

  model.nested <- rlmer(strength ~ (1 | batch/cask), data = Pastes)
  expect_identical(confintROB:::getClusters(model.nested), Pastes$batch)
  expect_identical(nlevels(confintROB:::getClusters(model.nested)), 10L)
})

test_that("a single grouping factor keeps the clusters and draws", {
  library(lme4)

  model.ss.ML <- lmer(Reaction ~ Days + (Days | Subject),
                      data = sleepstudy, REML = FALSE)
  expect_identical(confintROB:::getClusters(model.ss.ML),
                   getME(model.ss.ML, "flist")[[1]])
  expect_identical(confintROB:::getClusters(model.ss.ML, "Subject"),
                   getME(model.ss.ML, "flist")[[1]])
})

test_that("the jackknife leaves out one whole cluster at a time", {
  library(lme4)

  model <- lmer(strength ~ (1 | batch/cask), data = Pastes, REML = FALSE)
  coefs <- t(sapply(levels(Pastes$batch), function(b) {
    fit <- suppressMessages(update(model, data = Pastes[Pastes$batch != b, ]))
    confintROB:::getEstimates.lmerMod(fit)
  }))
  expected <- confintROB:::computeAcc(coefs)
  names(expected) <- confintROB:::createNames(model)
  acc <- suppressMessages(confintROB:::clusjack(model, .export = NULL))
  expect_equal(acc, expected)
})
