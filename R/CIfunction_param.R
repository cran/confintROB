# Title
# @title Parametric Confidence Interval
#
# @param model an object of class lmerMod, rlmerMod or varComprob
# @param nsim number of bootstrap samples, positive integer
# @param .export passed on to \code{\link[foreach]{foreach}}
#
# @return Returns a parametric Confidence Interval
#' @importFrom stats simulate formula quantile
#' @importFrom lme4 lmer
parametric <- function(model, nsim, max.tries, .export, data, varComprob.random, ...)
  UseMethod("parametric")

parametric.lmerMod <-
  function(model, nsim, max.tries = 100, .export, data, varComprob.random, ...) {
    class(model) <- "lmerMod"
    sample <- function(n) {
      simulate(model, n, use.u = TRUE)
    }
    fit <- createFitFunction.lmerMod(model)
    bootstrap(model, nsim, max.tries, .export, sample, fit, ...)
  }

parametric.rlmerMod <-
  function(model, nsim, max.tries = 100, .export, data, varComprob.random, ...) {
    parametric.lmerMod(model, nsim, max.tries, .export, data, varComprob.random, ...)
  }