## Coverage of the 95% parametric percentile bootstrap interval when the
## bootstrap samples are generated with simulate(..., use.u = TRUE) (random
## effects held fixed at their predicted values, only the errors are redrawn)
## and with simulate(..., use.u = FALSE) (random effects and errors are both
## redrawn, as in Mason, Cantoni & Ghisletta, 2021, section "Parametric
## Bootstrap", step 1).
##
## Data are simulated from a known linear mixed model without contamination
## and fitted with lmer(). Both settings are applied to the same data sets.
## Everything but the sample function is confintROB's own machinery
## (bootstrap(), createFitFunction.lmerMod()), so the only difference between
## the two settings is the value of use.u.
##
## Run from the package root with
##   Rscript inst/simulations/parametric-use.u-coverage.R
## Optional environment variables: NDATA (number of data sets per design,
## default 500), NSIM (bootstrap samples per data set, default 500), CORES.
## Results are written next to this script. The results do not depend on the
## number of cores.

suppressPackageStartupMessages({
  library(lme4)
  library(parallel)
  library(confintROB)
})

nData <- as.integer(Sys.getenv("NDATA", 500))
nsim <- as.integer(Sys.getenv("NSIM", 500))
cores <- as.integer(Sys.getenv("CORES", max(1, detectCores() - 1)))
outDir <- Sys.getenv("OUTDIR", "inst/simulations")
level <- 0.95
options(width = 200)
seed <- 20261010

## true parameters, in the order confintROB returns them
designs <- list(
  intercept = list(
    formula = y ~ x + (1 | g),
    truth = c("(Intercept)" = 1, x = 0.5, "SD (Intercept)" = 1, sigma = 1)
  ),
  slope = list(
    formula = y ~ x + (x | g),
    truth = c("(Intercept)" = 1, x = 0.5, "SD (Intercept)" = 1, "SD x" = 0.5,
              "Cor (Intercept), x" = 0.3, sigma = 1)
  )
)
nGroups <- c(10, 30)
x <- 0:5 # observations per group

generate <- function(design, nGroup) {
  truth <- designs[[design]]$truth
  d <- data.frame(g = factor(rep(seq_len(nGroup), each = length(x))),
                  x = rep(x, nGroup))
  b0 <- rnorm(nGroup)
  b1 <- numeric(nGroup)
  if (design == "slope") {
    rho <- truth[["Cor (Intercept), x"]]
    b1 <- truth[["SD x"]] * (rho * b0 + sqrt(1 - rho^2) * rnorm(nGroup))
  }
  b0 <- truth[["SD (Intercept)"]] * b0
  d$y <- truth[["(Intercept)"]] + b0[d$g] + (truth[["x"]] + b1[d$g]) * d$x +
    rnorm(nrow(d), sd = truth[["sigma"]])
  d
}

## local copy of confintROB:::parametric.lmerMod with use.u as a switch
percentileCI <- function(model, use.u) {
  sample <- function(n) {
    simulate(model, n, use.u = use.u)
  }
  fit <- confintROB:::createFitFunction.lmerMod(model)
  results <- confintROB:::bootstrap(model, nsim, 100, NULL, sample, fit)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  t(apply(results, 2, quantile, probs = probs, na.rm = TRUE))
}

runOne <- function(seed, design, nGroup) {
  set.seed(seed)
  d <- generate(design, nGroup)
  model <- lmer(designs[[design]]$formula, d, REML = FALSE)
  ci <- lapply(c("TRUE" = TRUE, "FALSE" = FALSE), function(use.u)
    tryCatch(percentileCI(model, use.u), error = function(e) NULL))
  list(estimates = confintROB:::getEstimates.lmerMod(model), ci = ci)
}

runCell <- function(design, nGroup, seeds) {
  res <- mclapply(seeds, function(s)
    suppressMessages(suppressWarnings(runOne(s, design, nGroup))),
    mc.cores = cores, mc.preschedule = FALSE)
  truth <- designs[[design]]$truth
  do.call(rbind, lapply(names(res[[1]]$ci), function(use.u) {
    ci <- lapply(res, function(r) r$ci[[use.u]])
    ok <- !vapply(ci, is.null, NA)
    lower <- vapply(ci[ok], function(m) m[, 1], truth)
    upper <- vapply(ci[ok], function(m) m[, 2], truth)
    coverage <- rowMeans(lower <= truth & truth <= upper)
    data.frame(design = design, nGroup = nGroup, use.u = use.u,
               parameter = names(truth), truth = unname(truth),
               nData = sum(ok), coverage = coverage,
               mcse = sqrt(coverage * (1 - coverage) / sum(ok)),
               meanLength = rowMeans(upper - lower), row.names = NULL)
  }))
}

set.seed(seed)
cells <- expand.grid(nGroup = nGroups, design = names(designs),
                     stringsAsFactors = FALSE)
cellSeeds <- lapply(seq_len(nrow(cells)), function(i)
  sample.int(.Machine$integer.max, nData))
summary <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  cat("design", cells$design[i], "with", cells$nGroup[i], "groups\n")
  runCell(cells$design[i], cells$nGroup[i], cellSeeds[[i]])
}))

attr(summary, "settings") <- list(
  nData = nData, nsim = nsim, level = level, seed = seed, x = x,
  designs = designs, R = R.version.string,
  lme4 = as.character(packageVersion("lme4")))
saveRDS(summary, file.path(outDir, "parametric-use.u-coverage.rds"))
write.csv(summary, file.path(outDir, "parametric-use.u-coverage.csv"),
          row.names = FALSE)
print(format(summary, digits = 3), row.names = FALSE)
