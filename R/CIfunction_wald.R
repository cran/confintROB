confintWald <- function(object, parm, level = .95) {
  effects <- getEffects(object)
  delta <- computeDelta(object, level)
  inf <- effects - delta
  sup <- effects + delta
  bounds <-
    paste(format(
      100 * c((1 - level) / 2, 1 - (1 - level) / 2),
      trim = TRUE,
      scientific = FALSE,
      digit = 3
    ), "%")
  result <- cbind(inf,sup)
  colnames(result) <- bounds
  if (!missing(parm)) {
    result <- result[parm, , drop = FALSE]
  }
  return(result)
}

computeDelta <- function(object, level = .95, ...)
  UseMethod("computeDelta")

computeDelta.lmerMod <- function(object, level = .95, ...) {
  summ <- summary(object)
  alpha <- 1 - level
  quantile <- qnorm(level + alpha / 2)
  delta <- summ$coefficients[, 2] * quantile
  return(delta)
}

computeDelta.rlmerMod <- function(object, level = .95, ...) {
  computeDelta.lmerMod(object, level, ...)
}

computeDelta.varComprob <- function(object, level = .95, ...) {
  summ <- summary(object)
  alpha <- 1 - level
  quantile <- qnorm(level + alpha / 2)
  deltaBeta <- summ$zTable[, 2] * quantile
  deltaEta <- sqrt(diag(object$vcov.eta)) * quantile
  delta <- c(deltaBeta, deltaEta)
  return(delta)
}

getEffects <- function(object, ...)
  UseMethod("getEffects")

getEffects.lmerMod <- function(object, ...) {
  effects <- fixef(object)
  return(effects)
}

getEffects.rlmerMod <-function(object, ...) {
  getEffects.lmerMod(object, ...)
}

getEffects.varComprob <- function(object, ...) {
  effects <- c(object$beta, object$eta)
  return(effects)
}