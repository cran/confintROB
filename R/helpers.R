checkRobustlmmVersion <- function(object) {
  if (is(object, "rlmerMod") &&
      packageVersion("robustlmm") < "3.1-1") {
    stop(
      "Please update package 'robustlmm'. In order for 'confintROB' to work ",
      "correctly, at least version 3.1-1 is needed."
    )
  }
}

checkClusterID <- function(object, clusterID, method, boot.type) {
  if (missing(clusterID)) {
    if (is(object, "varComprob") &&
        (method == "BCa" || boot.type == "wild")) {
      stop(
        "Argument 'clusterID' is needed to do BCa or wild bootstrap",
        " for varComprob objects"
      )
    }
  } else if (getOption("confintROB.checkClusterID.message.ignored", TRUE)) {
    message(
      "Ignoring argument 'clusterID' as it's not needed ",
      "for this combination of arguments"
    )
    options("confintROB.checkClusterID.message.ignored" = FALSE)
  }
}

checkData <- function(object, varComprob.data, method, boot.type) {
  if (missing(varComprob.data)) {
    if (is(object, "varComprob") &&
        (boot.type == "parametric" || boot.type == "wild")) {
      stop("Argument 'varComprob.data' is needed to compute bootstrap Confidence Intervals",
           " for varComprob objects")
    }
  } else if (getOption("confintROB.checkdata.message.ignored", TRUE)) {
    message("Ignoring argument 'varComprob.data' as it's not needed ",
            "for this combination of arguments")
    options("confintROB.checkdata.message.ignored" = FALSE)
  }
}

checkK <- function(object) {
  if (is(object, "varComprob")) {
    if (base::is.null(object$K)) {
      stop("Argument 'K' is needed in the original call",
           " for varComprob objects")
    }
  }
}

checkRandom <- function(object, varComprob.random, boot.type) {
  if (is(object, "varComprob")) {
    if (missing(varComprob.random) &&
      (boot.type == "parametric" || boot.type == "wild")) {
      stop("Argument 'varComprob.random' is needed to compute bootstrap Confidence Intervals",
           " for varComprob objects")
    }
  }
}

getRequiredPackages <- function(object) {
  packages <- "MASS"
  if (is(object, "rlmerMod")) {
    packages <- c(packages, "robustlmm")
  } else if (is(object, "varComprob")) {
    packages <- c(packages, "robustvarComp")
  } else {
    packages <- c(packages, "lme4")
  }
  return(packages)
}

sdcor <- function(object) {
  vc <- as.data.frame(VarCorr(object))
  randoms <- vc[, 5]
  names(randoms) <-
    paste("Sigma",
          replace_na(unlist(vc[, 1]), ""),
          replace_na(unlist(vc[, 2]), ""),
          replace_na(unlist(vc[, 3]), ""))
  return(randoms)
}

getEstimates <- function(object, ...)
  UseMethod("getEstimates")

getEstimates.lmerMod <- function(object, ...) {
  c(fixef(object), sdcor(object))
}

getEstimates.rlmerMod <- function(object, ...) {
  getEstimates.lmerMod(object, ...)
}

getEstimates.varComprob <- function(object, ...) {
  c(object$fixef, sigma2 = object$eta0, object$eta)
}

createNames <- function(object) {
  fixed <- fixef(object)
  randoms <- sdcor(object)
  return(c(names(fixed), names(randoms)))
}

fitLmer <- function(object, bdd, ..., y = formula(object)[2],
                    formulrest = as.character(formula(object))[3],
                    formulboot = paste(y, "~", formulrest)) {
  dots <- list(...)
  if (is.null(dots[["control"]])) {
    control <- getControl(object)
    return(lmer(
      formulboot,
      data = bdd,
      REML = FALSE,
      control = control,
      ...
    ))
  }
  return(lmer(formulboot, data = bdd, REML = FALSE, ...))
}

#' @importFrom stats getCall
#' @importFrom lme4 lmerControl
getControl <- function(object) {
  control <- getCall(object)[["control"]]
  if (!is(control, "lmerControl")) {
    control <- lmerControl()
  }
  return(control)
}

createFitFunction.lmerMod <- function(model) {
  function(yboot, ...) {
    bdd <- model@frame
    bdd$yboot <- yboot
    model.bootr  <- fitLmer(model, bdd, ..., y = "yboot")
    OK <- length(model.bootr@optinfo$conv$lme4$messages) == 0
    c(OK, getEstimates.lmerMod(model.bootr))
  }
}

createFitFunction.varComprob <- function(model, data, random) {
  formulboot <- combineFormulas(formula(model), random, "yboot")
  function(yboot, ...) {
    bdd <- data
    bdd$yboot <- yboot
    model.bootr <- fitLmer(model, bdd, ..., formulboot = formulboot)
    OK <- length(model.bootr@optinfo$conv$lme4$messages) == 0
    c(OK, getEstimates.lmerMod(model.bootr))
  }
}

combineFormulas <- function(fixefFormula,
                            randomFormula,
                            response = as.character(fixef)[2]) {
  fixef <- lastCharacterElement(fixefFormula)
  random <- lastCharacterElement(randomFormula)
  combinedFormula <- paste(response, "~", fixef, "+", random)
  return(combinedFormula)
}

lastCharacterElement <- function(formula) {
  char <- as.character(formula)
  return(char[length(char)])
}

bootstrap <-
  function(model,
           nsim,
           max.tries,
           .export,
           sample,
           fit,
           ...) {
    result <- NULL
    it <- 0
    remaining.nsim <- nsim
    
    while (remaining.nsim > 0 && (it <- it + 1) < max.tries) {
      samples  <- sample(remaining.nsim)
      itresult <-
        bootstrap.iteration(model, samples, .export, fit, ...)
      result <- rbind(result, itresult)
      remaining.nsim <- nsim - NROW(result)
    }
    
    if (remaining.nsim > 0) {
      if (nsim == 1) {
        stop("Failed to produce a valid model fit after ", it, " tries.")
      }
      warning("Failed to produce ",
              nsim,
              " valid model fits after ",
              it,
              " tries.")
    }
    return(result)
  }

bootstrap.iteration <- function(model, samples, .export, fit, ...) {
  `%foreachOp%` <- getForeachOperator()
  yboot <- NULL ## make R CMD CHECK happy
  resultr <- foreach(
    yboot = samples,
    .combine = "rbind",
    .packages = getRequiredPackages(model),
    .export = .export
  ) %foreachOp% fit(yboot, ...)
  
  if (NCOL(resultr) == 1) {
    if (resultr[1]) {
      return(resultr[-1])
    } else {
      return(NULL)
    }
  }
  resultr[resultr[, 1] == 1, -1, drop = FALSE]
}

getForeachOperator <- function() {
  if (foreach::getDoParRegistered()) {
    return(foreach::`%dopar%`)
  } else {
    return(foreach::`%do%`)
  }
}
