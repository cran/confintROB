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
  } else {
    message(
      "Ignoring argument 'clusterID' as it's not needed ",
      "for this combination of arguments"
    )
  }
}

checkData <- function(object, varComprob.data, method, boot.type) {
  if (missing(varComprob.data)) {
    if (is(object, "varComprob") &&
        (boot.type == "parametric" || boot.type == "wild")) {
      stop("Argument 'varComprob.data' is needed to compute bootstrap Confidence Intervals",
           " for varComprob objects")
    }
  } else {
    message("Ignoring argument 'varComprob.data' as it's not needed ",
            "for this combination of arguments")
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
  ) %foreachOp% suppressWarnings(fit(yboot, ...))
  
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

# Bake model call by evaluating variable references
#
# Replaces variable references in a model's call with their actual values.
# This ensures that update() works correctly even when the original variables
# are no longer in scope (e.g., during jackknife resampling).
#
# Uses R's built-in substitute() for robust symbol replacement.
#
# @param model an object of class lmerMod or rlmerMod
# @param envir environment to look up variables (default: formula environment)
# @return the model with its call modified to contain evaluated values
bakeCall <- function(model, envir = NULL) {
  cl <- getCall(model)
  if (is.null(cl)) return(model)

  # Build a substitution environment by collecting variables from relevant scopes
  # Priority: caller_env, formula environment, model frame environment
  subst_env <- buildSubstEnv(model, envir)

  # Arguments to skip (these are special or should remain as-is)
  skip_args <- c("", "data", "formula")

  # Process each argument in the call
  for (i in seq_along(cl)) {
    arg_name <- names(cl)[i]
    if (is.null(arg_name)) arg_name <- ""

    # Skip the function name (first element) and special arguments
    if (i == 1 || arg_name %in% skip_args) next

    arg <- cl[[i]]

    # If argument is a symbol (variable name), try to substitute it
    if (is.symbol(arg)) {
      name <- as.character(arg)
      if (exists(name, envir = subst_env, inherits = FALSE)) {
        cl[[i]] <- subst_env[[name]]
      }
    }
    # If argument is a call that defines a function (stored as expression)
    # e.g., init = function(...) lmer(..., control = control)
    else if (is.call(arg) && identical(arg[[1]], as.symbol("function"))) {
      # Use substitute() on the function body (3rd element)
      if (length(arg) >= 3) {
        arg[[3]] <- do.call(substitute, list(arg[[3]], subst_env))
        cl[[i]] <- arg
      }
    }
    # If argument is an actual function object, substitute in its body
    else if (is.function(arg)) {
      body(arg) <- do.call(substitute, list(body(arg), subst_env))
      cl[[i]] <- arg
    }
  }

  # Update the model's call
  if (is(model, "lmerMod") || is(model, "rlmerMod")) {
    model@call <- cl
  }

  return(model)
}

# Build substitution environment from model and caller context
#
# Collects variables from relevant environments into a single environment
# for use with substitute().
#
# @param model an object of class lmerMod or rlmerMod
# @param envir caller's environment (optional)
# @return environment containing variables for substitution
buildSubstEnv <- function(model, envir = NULL) {
  # Create new environment for substitution (don't inherit from anywhere)
  subst_env <- new.env(parent = emptyenv())

  # Collect source environments in reverse priority order
  # (later ones override earlier ones)
  source_envs <- list()

  # Global environment (lowest priority)
  source_envs <- c(source_envs, list(globalenv()))

  # Model frame environment
  if (is(model, "lmerMod") || is(model, "rlmerMod")) {
    frame_env <- attr(model@frame, ".Environment")
    if (!is.null(frame_env)) {
      source_envs <- c(source_envs, list(frame_env))
    }
  }

  # Formula environment
  formula_env <- environment(formula(model))
  if (!is.null(formula_env)) {
    source_envs <- c(source_envs, list(formula_env))
  }

  # Caller's environment (highest priority)
  if (!is.null(envir)) {
    source_envs <- c(source_envs, list(envir))
  }

  # Copy non-function, non-environment values to substitution environment
  for (src_env in source_envs) {
    for (name in ls(src_env, all.names = TRUE)) {
      tryCatch({
        val <- get(name, envir = src_env, inherits = FALSE)
        # Only substitute data values, not functions or environments
        if (!is.function(val) && !is.environment(val)) {
          subst_env[[name]] <- val
        }
      }, error = function(e) NULL)
    }
  }

  return(subst_env)
}
