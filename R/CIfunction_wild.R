# Title
# @title Wild Confidence Interval
#
# @param model an object of class lmerMod, rlmerMod or varComprob
# @param clusterID text variable indicating the clustering variable (only
#   needed to do wild bootstrap for varComprob objects)
# @param nsim number of bootstrap samples, positive integer
#
# @return Returns a wild Confidence Interval
#' @importFrom stats model.matrix terms
#' @importFrom foreach `%dopar%` foreach

wild <-
  function(model,
           clusterID,
           nsim,
           max.tries,
           .export,
           data,
           varComprob.random,
           ...)
    UseMethod("wild")

# @param v1 first possible value for weights in the wild bootstrap
# @param v2 second possible value for weights in the wild bootstrap
# @param TT length of cluster ID
# @param p1 probability associated to v1
# @param p2 probability associated to v2
# @param Int matrix I
# @param Pt orthogonal projection
# @param rt_hat bootstrap residuals
# @param Xt design matrix for fixed effects
# @param bet vector of fixed effects
# @param bdd dataset
# @param model an object of class rlmerMod or lmerMod
# @param rt_b emplacement for resiudals
# @param yt_b emplacement for y bootstrap
#' @importFrom MASS ginv
#' @importFrom lme4 getME
wild.lmerMod <-
  function(model,
           clusterID,
           nsim,
           max.tries = 100,
           .export,
           data,
           varComprob.random,
           ...) {
    y 		 <- getME(model, "y")
    X      <- as.matrix(getME(model, "X"))
    id     <- getME(model, "flist")[[1]]
    bet    <- unname(fixef(model))
    sample <- createWildSampleFunction(y, X, id, bet)
    fit    <- createFitFunction.lmerMod(model)
    bootstrap(model, nsim, max.tries, .export, sample, fit, ...)
  }

createWildSampleFunction <- function(y, X, id, bet) {
  # Define the probability of weights
  v1 		<- -(sqrt(5) - 1) / 2
  v2 		<- (sqrt(5) + 1) / 2
  p1 		<- (sqrt(5) + 1) / (2 * sqrt(5))
  p2 		<- 1 - p1
  
  # Define the disturbances vector for each participant
  ids    <- unlist(unique(id))
  TT 		 <- length(ids)
  X      <- as.matrix(X)
  tXX    <- solve(crossprod(X))
  Xbet   <- drop(X %*% unname(bet))
  rt_hat <- unname(y - Xbet)
  n 		 <- nrow(X)
  rt_hat_scaled <- numeric(n)
  idxMap <- integer(n)
  
  for (i in 1:TT) {
    idx <- id == ids[i]
    Xt  <- X[idx, , drop = FALSE]
    Pt  <- tcrossprod(Xt %*% tXX, Xt)
    Int <- diag(sum(idx))
    rt_hat_scaled[idx] <- sqrt(diag(ginv(Int - Pt))) * rt_hat[idx]
    idxMap[idx] <- i
  }
  
  oneSample <- function() {
    wt <- base::sample(c(v1, v2), TT, replace = TRUE, prob = c(p1, p2))
    Xbet + rt_hat_scaled * wt[idxMap]
  }
  sample <- function(n) {
    replicate(n, oneSample(), simplify = FALSE)
  }
  return(sample)
}

wild.rlmerMod <-
  function(model, clusterID, nsim, max.tries = 100, .export, data, varComprob.random, ...) {
    wild.lmerMod(model, clusterID, nsim, max.tries, .export, data, varComprob.random, ...)
  }