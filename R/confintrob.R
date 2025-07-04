#' Implements the classical Wald-type CI, the parametric and the wild bootstrap (Modugno & Giannerini, 2013) for linear mixed models estimated with the robust estimators in \code{\link[robustlmm]{rlmer}} (Koller, 2016; Koller & Stahel, 2022) and 
#' \code{\link[robustvarComp]{varComprob}} (Agostinelli & Yohai, 2016) and the classical estimators in \code{\link[lme4]{lmer}} (Bates et al., 2013). For bootstrap methods, percentile, Bias-Corrected, and accelerated (BCa) versions are available. All these versions are tested in Mason, Cantoni & Ghisletta (2021, 2024).
#' @description
#' Confidence Intervals (CIs) based on parametric and semi-parametric bootstrap (and Wald-type) for robust and classical Linear Mixed Models estimators.
#' @details \code{confintROB} computes 5 types of CIs based on arguments \code{method} and \code{boot.type}.
#' 
#' \code{method}:
#' 
#' - \code{"Wald"}: CIs computation is based on standard error estimates
#' 
#' - \code{"boot"}: CIs are computed using the respective percentile of estimates obtained on bootstrap sample(s) based on the confidence \code{level}
#' 
#' - \code{"BCa"}: based on the Jacknife method, the Bias-Corrected and accelerated parameters are computed on the original sample to correct estimates obtained from the bootstrap sample(s)
#' 
#' \code{boot.type} (for \code{method}s \code{"boot"} and \code{"BCa"}):
#' 
#' - \code{"parametric"}: the classical parametric bootstrap scheme is used to create bootstrap sample(s)
#' 
#' - \code{"wild"}: the semi-parametric bootstrap scheme is used to create bootstrap sample(s)
#' @references Agostinelli, C., & Yohai, V. J. (2016). Composite robust estimators for linear mixed models. Journal of the American Statistical Association, 111 (516), 1764-1774. doi:10.1080/01621459.2015.1115358
#' 
#' Bates, D., Machler, M., Bolker, B., & Walker, S. (2015). Fitting linear mixed-effects models using lme4. Journal of Statistical Software, 67 (1), 1-48. doi: 10.18637/jss.v067.i01
#' 
#' Koller, M. (2016). robustlmm: An R package for robust estimation of linear mixed-effects models. Journal of Statistical Software, 75 (6), 1-24. doi: 10.18637/jss.v075.i06
#' 
#' Koller, M., & Stahel, W. A. (2022). Robust estimation of general linear mixed effects models. In P. M. Yi & P. K. Nordhausen (Eds.), Robust and multivariate statistical methods. Springer Nature Switzerland AG.
#' 
#' Mason, F., Cantoni, E., & Ghisletta, P. (2021). Parametric and semi-parametric bootstrap-based confidence intervals for robust linear mixed models. Methodology, 17 (4), 271-295. doi: 10.5964/meth.6607
#' 
#' Mason, F., Cantoni, E., & Ghisletta, P. (2024). Linear mixed models and latent growth curve models for group comparison studies contaminated by outliers. Psychological methods. doi: 10.1037/met0000643
#' 
#' Modugno, L., & Giannerini, S. (2013). The wild bootstrap for multilevel models. Communications in Statistics - Theory and Methods, 44 (22), 4812-4825. doi: 10.1080/03610926.2013.80280
#' 
#' @seealso \code{\link[lme4]{lmer}}
#' \code{\link[robustlmm]{rlmer}} 
#' \code{\link[robustvarComp]{varComprob}}
#' 
#' @title Confidence Intervals for Robust and Classical Linear Mixed Model Estimators.
#' @param object an object of class \code{lmerMod}, \code{rlmerMod} or \code{varComprob}
#' @param parm parameters for which intervals are sought. Specified by an integer vector of positions (see example) or a character vector of parameter names. Fixed effects parameters are ordered according their appearance in the formula. For the order of variance components, see argument \code{order} of the \code{varCorr} function from package \code{\link[lme4]{lme4}}. By default, the CIs of all the parameters of the model are computed.
#' @param level confidence level in ]0; 1[
#' @param method type of CIs: \code{"Wald"}, \code{"boot"}, \code{"BCa"}
#' @param nsim number of bootstrap samples, positive integer
#' @param boot.type type of bootstrap: \code{"wild"} or \code{"parametric"}
#' @param clusterID text variable indicating the clustering variable. This is only required for method \code{"BCa"} or for boot.type \code{"wild"} for \code{\link[robustvarComp]{varComprob}} objects. This vector is used to identify the level of the "cluster" to which these resampling methods should be applied and is not included in the \code{\link[robustvarComp]{varComprob}} objects.
#' @param verify.saved check if an existing CI is already computed. Only for the vignette examples
#' @param .export passed on to \code{\link[foreach]{foreach}}
#' @param varComprob.data a data frame object used to fit the original model. This is only required for the \code{varComprob} objects
#' @param varComprob.random text variable describing the random effect structure as it would be written with \code{\link[lme4]{lmer}} from \code{\link[lme4]{lme4}}. This is only required for the \code{varComprob} objects
#' @param ... additional arguments passed on to \code{\link[lme4]{lmer}} (if applicable)
#'
#' @return a numeric table (matrix with column and row names) of CIs.
#' @export
#' @examples
#' if (require(robustlmm)) {
#'   model.RSE <- rlmer(Reaction ~ 1 + Days + (Days|Subject),
#'                     data = sleepstudy)
#'   CI.RSE <- confintROB(model.RSE, level = .95, boot.type = "wild",
#'                     parm = c(1,2), # indicates that only CIs
#'                                    # for the intercept and Days are asked.
#'                     nsim = 10) # small nsim for speed, in practice use, e.g., 5000
#'   print(CI.RSE)
#' }
#' @importFrom lme4 confint.merMod
#' @importFrom utils packageVersion
confintROB <-
  function(object,
           parm,
           level = .95,
           method = c("boot", "BCa", "Wald"),
           nsim = 5000,
           boot.type = c("wild", "parametric"),
           clusterID,
           verify.saved = NULL,
           .export = NULL,
           varComprob.data,
           varComprob.random,
           ...) {
    method <- match.arg(method)
    boot.type <- match.arg(boot.type)
    
    if (method == "Wald") {
      result <- confintWald(object, parm, level)
      return(result)
    }
    
    checkRobustlmmVersion(object)
    checkClusterID(object, clusterID, method, boot.type)
    checkData(object, varComprob.data, method, boot.type)
    checkK(object)
    checkRandom(object, varComprob.random, boot.type)
      
    test.sampling <-
      confintLMM(object,
                 clusterID,
                 boot.type,
                 1,
                 level,
                 .export,
                 test.sampling = NULL,
                 varComprob.data,
                 varComprob.random,
                 ...)
    if (isValidSavedResults(verify.saved, test.sampling, nsim)) {
      return(verify.saved)
    }
    
    fullResults <-
      confintLMM(object,
                 clusterID,
                 boot.type,
                 nsim - 1,
                 level,
                 .export,
                 test.sampling,
                 varComprob.data,
                 varComprob.random,
                 ...)
    
    retVal <- fullResults[[1]]
    if (method == "BCa") {
      ci_BCa <-
        try(confint_BCa(object, nsim, clusterID, fullResults[[2]], level, .export, varComprob.data))
      if (is(ci_BCa, "try-error")) {
        warning("Method 'BCa' failed, returning bootstrap estimates")
      } else {
        retVal <- ci_BCa[[1]]
        fullResults <- c(ci_BCa[1], fullResults, ci_BCa[-1])
      }
    } else if (method != "boot") {
      warning("Argument 'method' was ignored")
    }
    
    class(fullResults) <- "confintROB-fullResults"
    if (!missing(parm)) {
      retVal <- retVal[parm, , drop = FALSE]
    }
    attr(retVal, "fullResults") <- fullResults
    return(retVal)
  }

isValidSavedResults <- function(verify.saved, test.sampling, nsim) {
  if (is.null(verify.saved)) {
    return(FALSE)
  }
  if (nsim != nrow(attr(verify.saved, "fullResults")$bootstrap_estimates)) {
    return(FALSE)
  }
  expected <-
    attr(verify.saved, "fullResults")$bootstrap_estimates[1, ]
  if (is.null(expected)) {
    return(FALSE)
  }
  actual <- test.sampling$bootstrap_estimates[1, ]
  check <-
    all.equal(expected,
              actual,
              check.attributes = FALSE,
              tolerance = 1e-9)
  return(isTRUE(check))
}

#' @export
`print.confintROB-fullResults` <- function(x, ...) {
  cat(
    "  Full results of confintROB, a list with components:\n  ",
    paste("\"", names(x), "\"", collapse = ", ", sep = ""),
    "\n"
  )
  invisible(x)
}
