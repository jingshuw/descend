#' Default is to set max.quantile = max.sparse = 0.98. For those that are too sparse, set max.quantile = max.sparse = 0.99 and try again.
#'
#' @inheritParams deconvSingle
#' @param count.matrix the observed UMI count matrix. Each row is a gene and each column is a cell
#' @param n.cores the number of cores used for parallel computing. Default is 1. Used only when parallel computing is done in a single machine. For using multi-machine cores, need to assign \code{cl} explicitly
#' @param cl an object of class "cluster". See more details in \code{\link[parallel]{makeCluster}}
#' 
#' @return a list of DESCEND objects
#'
#' @export

descend <- function(count.matrix,
                    scaling.consts = NULL,
                    Z = NULL,
                    Z0 = NULL,
                    n.cores = 1, cl = NULL,
                    do.LRT.test = F, 
                    family = c("Poisson", "Negative Binomial"),
                    NB.size = 100,
                    verbose = T, 
                    control = list()) {

   control <- do.call("DESCEND.control", control)

  Y <- as.matrix(count.matrix)
  gene.names <- rownames(Y)

  print(paste("DESCEND starts deconvolving distribution of", length(gene.names), "genes!"))

  if (is.null(cl)) {
    require(doParallel)
    registerDoParallel(n.cores)

    results <-  foreach(v = iter(Y, "row"), .noexport = c("Y")) %dopar% {
      if (verbose)
        print("Compute for gene")
      if (mean(v == 0) > control$max.sparse) {
        control$max.quantile <- 0.99
        control$max.sparse <- 0.99
      }

      temp <- try(deconvSingle(v,
                               scaling.consts = scaling.consts,
                               Z = Z, Z0 = Z0,
                               plot.density = F,
                               do.LRT.test = do.LRT.test,
                               family = family, 
                               NB.size = NB.size,
                               control = control, verbose = F))
      if (class(temp) == "try-error")
        return(NA)
      return(temp)
    }

  } else if (!is.null(cl)) {
    Y <- as.list(as.data.frame(t(Y)))
    clusterExport(cl, c("scaling.consts", "Z", "Z0", "do.LRT.test", "family", "NB.size",
                        "control"), environment())

    results <- parLapply(cl, Y, function(v) {

                         source("~/Dropbox/sparse_factor_bic/code/g_model/package_code/deconvSingle.R")
                         source("~/Dropbox/sparse_factor_bic/code/g_model/package_code/g_model.R")

                         if ("logMsg" %in% ls())
                           logMsg(do.LRT.test)

                         if (mean(v == 0) > control$max.sparse) {
                           control$max.quantile <- 0.99
                           control$max.sparse <- 0.99
                         }

                         temp <- try(deconvSingle(v,
                                                  scaling.consts = scaling.consts,
                                                  Z = Z, Z0 = Z0,
                                                  plot.density = F,
                                                  do.LRT.test = do.LRT.test,
                                                  family = family, 
                                                  NB.size = NB.size,
                                                  control = control, verbose = F))
                         if (class(temp) == "try-error")
                           return(NA)
                         return(temp)
                    })
  } 
  names(results) <- gene.names

  return(results)
}

#' Grep the value and standard deviation of the estimated parameters calculated from a list of descend objects 
#'
#' @param descend.list a list of descend objects
#' 
#' @export

getEstimates <- function(descend.list) {
  temp <- lapply(descend.list, function(ll) {
                      if (class(ll) == "DESCEND")
                        return(ll@estimates[, c(1, 4)])
                      else 
                        return(NA)
                    })
  idx <- which(!sapply(temp, function(mat)is.na(mat[1])))[1]
  n.est <- nrow(temp[[idx]])
  est.names <- rownames(temp[[idx]])
  report.names <- colnames(temp[[idx]])

  temp <- lapply(temp, function(mat) {
                      if (is.na(mat[1])) {
                        mm <- matrix(NA, n.est, 2)
                        colnames(mm) <- report.names
                        rownames(mm) <- est.names
                        return(mm)
                      }
                      else
                        return(mat)
                    })

  est.list <- lapply(1:length(est.names), function(i) {
                      result <- t(sapply(temp, function(mat) mat[i, ]))
                      colnames(result) <- report.names
                      rownames(result) <- names(descend.list)
                      return(result)
                    })

  names(est.list) <- est.names
  return(est.list)
}

getPval <- function(descend.list) {
  temp <- lapply(descend.list, function(ll) {
                      if (class(ll) == "DESCEND")
                        return(ll@pval)
                      else 
                        return(NA)
                    })

  idx <- which(!sapply(temp, is.na))[1]
  n.test <- nrow(temp[[idx]])
  if (n.test == 0) 
    return(NA)
  test.names <- rownames(temp[[idx]])
  report.names <- colnames(temp[[idx]])

  temp <- lapply(temp, function(mat) {
                      if (is.na(mat[1])) {
                        mm <- matrix(NA, n.est, 1)
                        colnames(mm) <- report.names
                        rownames(mm) <- test.names
                      }
                      else
                        return(mat)
                    })


  test.result <- t(sapply(temp, function(mat) mat[, 1]))

  rownames(test.result) <- names(descend.list)

  return(test.result)

}
