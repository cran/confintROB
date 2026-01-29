# createParamSampleFunction works with varComprob S-estimator

    Code
      print(summary(model.S), digits = 2)
    Output
      Method:  S

      Fixed effects:
                  Estimate Std. Error z value Pr(>|z|)
      (Intercept)    253.8        6.6    38.2   <2e-16 ***
      Days            10.8        1.6     6.9    7e-12 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

       Random effect variances:
       sigma2_u0 Covariance  sigma2_u1
      1104.97202  -49.27180   52.98219

       Residual variance:
      error variance
            692.5566

       Value of the objective function at convergence:
      [1] 10752.14


      Robustness weights:
       6 observations c(1,6,7,11,15,18) are outliers with |weight| = 0 ( < 0.005556);
       The remaining 12 ones are summarized as
           Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
      0.3215799 0.4935559 0.9459881 0.8379223 1.0917200 1.2086610

---

    Code
      cat("Running test for object of class", paste(class(model.S), collapse = " "),
      "\n")
    Output
      Running test for object of class varComprob.S varComprob.fit varComprob
    Code
      sample <- confintROB:::createParamSampleFunction(model = model.S, data = sleepstudy)
      set.seed(123)
      result11 <- c(sample(1), sample(1))
      set.seed(123)
      result2 <- sample(2)
      names(result11) <- names(result2)
      stopifnot(all.equal(result11, result2))
      print(head(result2[[1]]), digits = 5)
    Output
           1      2      3      4      5      6
      251.71 322.29 206.80 325.15 274.34 329.96
    Code
      print(head(result2[[2]]), digits = 5)
    Output
           1      2      3      4      5      6
      220.15 261.71 366.48 295.06 290.76 366.08

# createParamSampleFunction works with varComprob composite-TAU estimator

    Code
      print(summary(model.cTAU), digits = 2)
    Output
      Method:  compositeTau

      Fixed effects:
                  Estimate Std. Error z value Pr(>|z|)
      (Intercept)    250.9        7.8    32.2   <2e-16 ***
      Days            10.2        1.5     6.7    3e-11 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

       Random effect variances:
       sigma2_u0 Covariance  sigma2_u1
       819.02411  -26.54169   49.76142

       Residual variance:
      error variance
            376.8007

       Value of the objective function at convergence:
      [1] 7675.007

---

    Code
      cat("Running test for object of class", paste(class(model.cTAU), collapse = " "),
      "\n")
    Output
      Running test for object of class varComprob.compositeTau varComprob.fit varComprob
    Code
      sample <- confintROB:::createParamSampleFunction(model = model.cTAU, data = sleepstudy)
      set.seed(123)
      result11 <- c(sample(1), sample(1))
      set.seed(123)
      result2 <- sample(2)
      names(result11) <- names(result2)
      stopifnot(all.equal(result11, result2))
      print(head(result2[[1]]), digits = 5)
    Output
           1      2      3      4      5      6
      251.66 308.23 213.82 315.15 274.40 321.11
    Code
      print(head(result2[[2]]), digits = 5)
    Output
           1      2      3      4      5      6
      225.00 257.90 344.42 289.68 289.04 346.64

# createParamSampleFunction works with permuted data

    Code
      cat("Permutation head:", head(permutation), "\n")
    Output
      Permutation head: 68 167 129 162 43 14

---

    Code
      print(summary(model.cTAU_permuted), digits = 2)
    Output
      Method:  compositeTau

      Fixed effects:
                  Estimate Std. Error z value Pr(>|z|)
      (Intercept)    250.9        7.8    32.2   <2e-16 ***
      Days            10.2        1.5     6.7    3e-11 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

       Random effect variances:
       sigma2_u0 Covariance  sigma2_u1
       819.02411  -26.54169   49.76142

       Residual variance:
      error variance
            376.8007

       Value of the objective function at convergence:
      [1] 7675.007

---

    Code
      cat("Running test for permuted data\n")
    Output
      Running test for permuted data
    Code
      sample_permuted <- confintROB:::createParamSampleFunction(model = model.cTAU_permuted,
        data = data_permuted)
      set.seed(123)
      result_permuted <- sample_permuted(2)
      print(head(result_permuted[[1]]), digits = 5)
    Output
          68    167    129    162     43     14
      315.19 319.75 307.52 269.58 272.30 318.29
    Code
      print(head(result_permuted[[2]]), digits = 5)
    Output
          68    167    129    162     43     14
      323.84 268.18 335.08 214.43 293.81 309.89

