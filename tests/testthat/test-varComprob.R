# Tests for varComprob objects from robustvarComp package
# These tests require robustvarComp (in Suggests)

test_that("createParamSampleFunction works with varComprob S-estimator", {
  skip_if_not_installed("robustvarComp")
  library(lme4)
  library(robustvarComp)

  participant <- sleepstudy$Subject
  within <- sleepstudy$Days

  # Build the argument "groups" of the varComprob() function
  n <- length(unique(participant))
  J <- length(unique(within))
  groups <- cbind(rep(1:J, each = n), rep((1:n), J))

  # Build the argument "varcov" of the varComprob() function
  z1 <- rep(1, J)
  z2 <- unique(within)

  K <- list(
    sigma2_u0 = tcrossprod(z1, z1),
    Covariance = tcrossprod(z1, z2) + tcrossprod(z2, z1),
    sigma2_u1 = tcrossprod(z2, z2)
  )

  # Estimation with S-estimator
  suppressWarnings(
    model.S <- varComprob(
      Reaction ~ 1 + Days,
      groups = groups,
      data = sleepstudy,
      varcov = K,
      control = varComprob.control(
        lower = c(0, -Inf, 0),
        method = "S",
        psi = "rocke",
        max.it = 1,
        init = list(
          beta = c("(Intercept)" = 253.835569743834, Days = 10.7736608268214),
          gamma = c(
            sigma2_u0 = 1.59549700005736,
            Covariance = -0.0711447985744645,
            sigma2_u1 = 0.0765023178239254
          ),
          eta0 = c("error variance" = 692.556625895202),
          scale = 10752.1432565101
        )
      )
    )
  )

  expect_snapshot({
    print(summary(model.S), digits = 2)
  })

  # Test createParamSampleFunction
  expect_snapshot({
    cat("Running test for object of class", paste(class(model.S), collapse = " "), "\n")
    sample <- confintROB:::createParamSampleFunction(model = model.S, data = sleepstudy)
    set.seed(123)
    result11 <- c(sample(1), sample(1))
    set.seed(123)
    result2 <- sample(2)
    names(result11) <- names(result2)
    stopifnot(all.equal(result11, result2))
    print(head(result2[[1]]), digits = 5)
    print(head(result2[[2]]), digits = 5)
  })
})

test_that("createParamSampleFunction works with varComprob composite-TAU estimator", {
  skip_if_not_installed("robustvarComp")
  library(lme4)
  library(robustvarComp)

  participant <- sleepstudy$Subject
  within <- sleepstudy$Days

  n <- length(unique(participant))
  J <- length(unique(within))
  groups <- cbind(rep(1:J, each = n), rep((1:n), J))

  z1 <- rep(1, J)
  z2 <- unique(within)

  K <- list(
    sigma2_u0 = tcrossprod(z1, z1),
    Covariance = tcrossprod(z1, z2) + tcrossprod(z2, z1),
    sigma2_u1 = tcrossprod(z2, z2)
  )

  control <- varComprob.control(
    lower = c(0, -Inf, 0),
    max.it = 1,
    init = list(
      beta = c("(Intercept)" = 250.945321738908, Days = 10.2320816031076),
      gamma = c(
        sigma2_u0 = 2.17362686604633,
        Covariance = -0.0704396118106077,
        sigma2_u1 = 0.132062984417908
      ),
      eta0 = c("error variance" = 376.800691794604),
      scales = c(
        293.57715136143, 358.95262673052, 465.547583256656,
        561.3346991483, 692.21765047862, 932.623947285384,
        641.528419359161, 846.716921562313, 924.543567137878,
        365.994312558323, 481.953914967322, 585.564052671342,
        697.829285167833, 1009.71707572247, 672.461886751178,
        982.606142686251, 936.132126983003, 248.037407578449,
        374.605889784185, 536.450389280523, 854.773265534817,
        632.866330961722, 855.224511580672, 962.333779104256,
        391.221328441633, 629.884894368671, 834.926952170133,
        882.869865599689, 1022.24447287146, 1168.56340641807,
        575.172734225926, 715.931584462354, 671.517853836347,
        949.863650035998, 1052.4253043978, 760.626391277738,
        523.076365944673, 681.762701791185, 943.357505068095,
        914.246654077684, 856.56616457374, 1309.32923881337,
        717.252457844454, 685.620374481247, 781.788840625603
      )
    )
  )

  suppressWarnings(
    model.cTAU <- varComprob(
      Reaction ~ 1 + Days,
      groups = groups,
      data = sleepstudy,
      varcov = K,
      control = control
    )
  )

  expect_snapshot({
    print(summary(model.cTAU), digits = 2)
  })

  expect_snapshot({
    cat("Running test for object of class", paste(class(model.cTAU), collapse = " "), "\n")
    sample <- confintROB:::createParamSampleFunction(model = model.cTAU, data = sleepstudy)
    set.seed(123)
    result11 <- c(sample(1), sample(1))
    set.seed(123)
    result2 <- sample(2)
    names(result11) <- names(result2)
    stopifnot(all.equal(result11, result2))
    print(head(result2[[1]]), digits = 5)
    print(head(result2[[2]]), digits = 5)
  })
})

test_that("createParamSampleFunction works with permuted data", {
  skip_if_not_installed("robustvarComp")
  library(lme4)
  library(robustvarComp)

  participant <- sleepstudy$Subject
  within <- sleepstudy$Days

  n <- length(unique(participant))
  J <- length(unique(within))
  groups <- cbind(rep(1:J, each = n), rep((1:n), J))

  z1 <- rep(1, J)
  z2 <- unique(within)

  K <- list(
    sigma2_u0 = tcrossprod(z1, z1),
    Covariance = tcrossprod(z1, z2) + tcrossprod(z2, z1),
    sigma2_u1 = tcrossprod(z2, z2)
  )

  control <- varComprob.control(
    lower = c(0, -Inf, 0),
    max.it = 1,
    init = list(
      beta = c("(Intercept)" = 250.945321738908, Days = 10.2320816031076),
      gamma = c(
        sigma2_u0 = 2.17362686604633,
        Covariance = -0.0704396118106077,
        sigma2_u1 = 0.132062984417908
      ),
      eta0 = c("error variance" = 376.800691794604),
      scales = c(
        293.57715136143, 358.95262673052, 465.547583256656,
        561.3346991483, 692.21765047862, 932.623947285384,
        641.528419359161, 846.716921562313, 924.543567137878,
        365.994312558323, 481.953914967322, 585.564052671342,
        697.829285167833, 1009.71707572247, 672.461886751178,
        982.606142686251, 936.132126983003, 248.037407578449,
        374.605889784185, 536.450389280523, 854.773265534817,
        632.866330961722, 855.224511580672, 962.333779104256,
        391.221328441633, 629.884894368671, 834.926952170133,
        882.869865599689, 1022.24447287146, 1168.56340641807,
        575.172734225926, 715.931584462354, 671.517853836347,
        949.863650035998, 1052.4253043978, 760.626391277738,
        523.076365944673, 681.762701791185, 943.357505068095,
        914.246654077684, 856.56616457374, 1309.32923881337,
        717.252457844454, 685.620374481247, 781.788840625603
      )
    )
  )

  # Original model
  suppressWarnings(
    model.cTAU <- varComprob(
      Reaction ~ 1 + Days,
      groups = groups,
      data = sleepstudy,
      varcov = K,
      control = control
    )
  )

  sample_original <- confintROB:::createParamSampleFunction(model = model.cTAU, data = sleepstudy)
  set.seed(123)
  result_original <- sample_original(2)

  # Permuted dataset
  set.seed(1)
  permutation <- sample.int(nrow(sleepstudy))
  groups_permuted <- groups[permutation, ]
  data_permuted <- sleepstudy[permutation, ]

  expect_snapshot({
    cat("Permutation head:", head(permutation), "\n")
  })

  suppressWarnings(
    model.cTAU_permuted <- varComprob(
      Reaction ~ 1 + Days,
      groups = groups_permuted,
      data = data_permuted,
      varcov = K,
      control = control
    )
  )

  expect_snapshot({
    print(summary(model.cTAU_permuted), digits = 2)
  })

  expect_snapshot({
    cat("Running test for permuted data\n")
    sample_permuted <- confintROB:::createParamSampleFunction(
      model = model.cTAU_permuted,
      data = data_permuted
    )
    set.seed(123)
    result_permuted <- sample_permuted(2)
    print(head(result_permuted[[1]]), digits = 5)
    print(head(result_permuted[[2]]), digits = 5)
  })

  # Verify permutation equivalence
  result_expected <- lapply(result_original, `[`, permutation)
  sample_permuted <- confintROB:::createParamSampleFunction(
    model = model.cTAU_permuted,
    data = data_permuted
  )
  set.seed(123)
  result_permuted <- sample_permuted(2)
  expect_equal(result_expected, result_permuted)
})
