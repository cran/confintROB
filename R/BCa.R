# Title
# @title Jackknife estimates
#
# @param model an object of class lmerMod or rlmerMod
# @param clusterID text variable indicating the clustering variable. Only
#   required for varComprob objects
# @param .export passed on to \code{\link[foreach]{foreach}}
#
# @return Jackknife estimates and accelerated parameter
#' @importFrom lme4 fixef VarCorr
#' @importFrom stats update
#' @importFrom tidyr replace_na
#' @importFrom foreach `%dopar%` foreach
# @export
clusjack <- function (model, clusterID, .export, data, ...)
  UseMethod("clusjack")

clusjack.lmerMod <- function(model, clusterID, .export, data, ...) {
  summ       <- summary(model)
  data       <- model@frame
  nomid      <- names(summ$ngrps)[1]
  idvect     <- data[nomid]
  n          <- nrow(data)
  
  cluster    <- table(idvect)
  clusters   <- unlist(labels(cluster))
  nc         <- length(clusters)
  Obsno      <- unlist(split(1:n, cluster))
  
  coefs      <- jackesampling(nc, Obsno, model, data, .export)
  acc        <- computeAcc(coefs)
  names(acc) <- createNames(model)
  return(acc)
}

clusjack.rlmerMod <- function(model, clusterID, .export, ...) {
  clusjack.lmerMod(model, clusterID, .export, ...)
}

# Title
# @title estimates on Jacknife samples with lmer and rlmer
#
# @param nc length of the cluster ID
# @param Obsno index of the clusers
# @param model an object of class lmerMod or rlmerMod
# @param data an object of class data.frame
# @param .export passed on to \code{\link[foreach]{foreach}}
#
# @return Jackknife lmer estimates
# @export
jackesampling <-
  function(nc, Obsno, model, data, .export) {
    i <- NULL ## Make R CMD check happy
    `%foreachOp%` <- getForeachOperator()
    foreach(
      i = 1:nc,
      .combine = "rbind",
      .packages = getRequiredPackages(model),
      .export = .export
    ) %foreachOp% {
      obs       <- unlist(Obsno[-i])
      modeljack <- update(model, data = data[obs, ])
      jackcoef  <- getEstimates.lmerMod(modeljack)
      
    }
  }

computeAcc <- function(coefs) {
  meanCoefs <- colMeans(coefs, na.rm = TRUE)
  uu        <- meanCoefs - t(coefs)
  uu2       <- uu * uu
  return(rowSums(uu2 * uu, na.rm = TRUE) / (6 * rowSums(uu2, na.rm = TRUE) ^ 1.5))
}

# Title
# @title BCa Confidence Interval
#
# @param model an object of class lmerMod, rlmerMod or varComprob
# @param nsim number of bootstrap samples, positive integer
# @param data an object of class data.frame
# @param clusterID text variable indicating the clustering variable.
# @param coefs matrix of coefficient estimates in the bootstrap samples
# @param confint.level confidence level < 1
# @param .export passed on to \code{\link[foreach]{foreach}}
#
# @return BCa Confidence Interval
#' @importFrom stats qnorm pnorm
# @export
confint_BCa <-
  function(model,
           nsim,
           clusterID,
           coefs,
           confint.level,
           .export,
           data,
           ...) {
    acc <- clusjack(model, clusterID, .export, data, ...)
    estimates <- getEstimates(model)
    p <- length(estimates)
    biascorr <-
      qnorm(rowSums(t(coefs) < estimates, na.rm = TRUE) / nsim) #remplacé nsim par length(unique(clusterID))?
    ci_BCa <- matrix(NA_real_, nrow = p, ncol = 2)
    quantiles <- qnorm((1 - confint.level) / 2) * c(1, -1)
    for (i in 1:p) {
      adjusted.quantiles <-
        biascorr[i] + (biascorr[i] + quantiles) / (1 - acc[i] * (biascorr[i] + quantiles))
      probs <- pnorm(adjusted.quantiles)
      ci_BCa[i, ] <- quantile(coefs[, i], probs, type = 1)
    }
    bounds <-
      paste(format(
        100 * c((1 - confint.level) / 2, 1 - (1 - confint.level) / 2),
        trim = TRUE,
        scientific = FALSE,
        digit = 3
      ), "%")
    colnames(ci_BCa) <- bounds
    row.names(ci_BCa) = names(estimates)
    return(list(
      BCa = ci_BCa,
      biasBCa = biascorr,
      acc = acc
    ))
  }