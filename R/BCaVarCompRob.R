# Title
# @title Jackknife estimates
#
# @param model an object of class varComprob
# @param data an object of class data.frame
# @param clusterID text variable indicating the clustering variable
# @param .export passed on to \code{\link[foreach]{foreach}}
#
# @return Jackknife estimates
# @export
clusjack.varComprob <-
  function (model, clusterID, .export, data,...) {
    idvect <- data[clusterID]
    n_obs <- nrow(data)
    estimates <- getEstimates.varComprob(model)
    p <- length(estimates)
    cluster <- table(idvect)
    clusters <- unlist(labels(cluster))
    nc <- length(clusters)
    
    coefs <- matrix(NA, nrow = nc, ncol = p)
    Obsno <- split(1:n_obs, model$model$`(groups)`[, 2])
    ii <-
      NULL ## to avoid warning about no visible binding in R CMD check
    `%foreachOp%` <- getForeachOperator()
    coefs = foreach(
      ii = 1:nc,
      .combine = "rbind",
      .packages = getRequiredPackages(model),
      .export = .export
    ) %foreachOp% {
      obs <- unlist(Obsno[-ii])
      i <- which(model$model$`(groups)`[, 2] == ii)
      groups_i <- model$model$`(groups)`[-i, ]
      data_i <- data[obs,]
      modeljack <-
        eval(substitute(
          update(model, data = data, groups = groups),
          list(data = data_i, groups = groups_i)
        ))
      jackcoef <- getEstimates.varComprob(modeljack)
    }
    colnames(coefs) <- names(estimates)
    acc <- computeAcc(coefs)
    names(acc) <- names(estimates)
    return(acc)
    
  }
