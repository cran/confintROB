# Title
# @title Wild Confidence Interval for varComprob objects
#
# @param model an object of class varComprob
# @param clusterID text variable indicating the clustering variable
# @param nsim number of bootstrap samples, positive integer
# @param max.tries number of times to try to produce a valid model fit
#   before giving up
# @param .export passed on to \code{\link[foreach]{foreach}}
#
# @return Returns a wild Confidence Interval
# @export
wild.varComprob <-
  function(model,
           clusterID,
           nsim,
           max.tries = 100,
           .export,
           data,
           varComprob.random,
           ...) {
    bdd	   <- model$model
    myvar  <- all.vars(model$terms)
    y 		 <- as.matrix(bdd[myvar[1]])
    id     <- data[[clusterID]]
    bet    <- as.vector(unname(model$fixef))
    X 	   <- model.matrix(model, data)
    sample <- createWildSampleFunction(y, X, id, bet)
    fit    <- createFitFunction.varComprob(model,data,varComprob.random)
    bootstrap(model, nsim, max.tries, .export, sample, fit)
  }

