# Title
# @title Confidence Interval for lmerMod, rlmerMod or varComprob objects
#
# @param model an object of class lmerMod, rlmerMod or varComprob
# @param clusterID text variable indicating the clustering variable (only
#   needed to do wild bootstrap for varComprob objects)
# @param method type of Confidence intervals: "parametric" or "wild" are possible
# @param nsim number of bootstrap samples, positive integer
# @param confint.level confidence level < 1
# @param .export passed on to \code{\link[foreach]{foreach}}
# @param test.sampling: previous sample results to be combined with new ones
# @param ... additional options passed on to \code{lmer} (if applicable)
# @return Returns a Confidence Interval
# @export
confintLMM <-
  function(model,
           clusterID,
           method,
           nsim,
           confint.level,
           .export,
           test.sampling,
           varComprob.data,
           varComprob.random,
           ...) {
    if (method == "parametric") {
      results <-
        parametric(model, nsim, .export = .export, data = varComprob.data, varComprob.random = varComprob.random, ...)
    } else if (method == "wild") {
      results <- wild(model, clusterID, nsim, .export = .export, data = varComprob.data, varComprob.random = varComprob.random, ...)
    } else {
      stop("Method not implemented for ", method, " type of bootstrap")
    }
    
    if (!is.null(test.sampling[["bootstrap_estimates"]])) {
      results <- rbind(test.sampling[["bootstrap_estimates"]], results)
    }
    
    if (nrow(results) == 1) {
      CI <- t(rbind(results, results))
    } else{
      probs <- c((1 - confint.level) / 2, 1 - (1 - confint.level) / 2)
      CI <-
        t(apply(results, 2, quantile, probs = probs, na.rm = TRUE))
    }
    bounds <-
      paste(format(
        100 * c((1 - confint.level) / 2, 1 - (1 - confint.level) / 2),
        trim = TRUE,
        scientific = FALSE,
        digit = 3
      ), "%")
    colnames(CI) <- bounds
    return(list(Percentile = CI, bootstrap_estimates = results))
  }