#' Apply descend to each gene of the count matrix 
#'
#' @inheritParams deconvSingle
#' @param count.matrix the observed UMI count matrix. Each row is a gene and each column is a cell
#' @param n.cores the number of cores used for parallel computing. Default is 1. Used only when parallel computing is done in a single machine. For using multi-machine cores, need to assign \code{cl} explicitly
#' @param cl an object of class "cluster". See more details in \code{\link[parallel]{makeCluster}}
#'
#'
#' @return a list of DESCEND objects
#'
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

#' Grab the value and standard deviation of the estimated parameters calculated from a list of descend objects 
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

#' Grab the likelihood ratio test p-values of the coefficients and active fraction is the tests are performed from a list of descend objects 
#'
#' @param descend.list a list of descend objects
#' 
#' @export

getPval <- function(descend.list) {
  temp <- lapply(descend.list, function(ll) {
                      if (class(ll) == "DESCEND")
                        return(ll@pval)
                      else 
                        return(NA)
                    })

  idx <- which(!sapply(temp, is.na))[1]
  n.test <- nrow(temp[[idx]])
  if (n.test == 0) {
    print("The likelihood ratio tests were not performed")
    return(NA)
  }
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

#' Finding highly variable genes (HVG) based on Gini coefficients or CV
#'
#' The function finds highly variable genes based on the Gini or CV estimates from DESCEND. A quantile natural cubic spline regression of the estimated selected criteria on the log of the estimated mean for the genes is performed as the general trend. HVGs are selected as the genes whose normalized difference between its estimates and the trend is larger than a threshold.
#'
#' @inheritParams getEstimates
#' @param criteria the cirteria of HVG finding, either "Gini" or "CV"
#' @param quantile quantile of the quantile regression. Defined as the \code{tau} parameter in the function \code{\link[quantreg]{rqss}}. A larger quantile yields a more stringent selection of HVG. Default is 0.7.
#' @param threshold threshold of the normalized difference for HVG selection. Default is 3. A larger value results in a more stringent selection
#' @param plot.result whether plot the selection results or not. Default is TRUE.
#'
#' @return A list of elements:
#' \item{score.mat}{A score matrix of the genes. Each row is a gene}
#' \item{HVG.genes}{A vector of the selected HVG gene names}
#'
#' @export 


findHVG <- function(descend.list, 
                    criteria = c("Gini", "CV"),
                    quantile = 0.7,
                    threshold = 3,
                    plot.result = T) {
  require(quantreg)
  require(splines)
  cirteria <- match.arg(criteria, c("Gini", "CV"))
  ests <- getEstimates(descend.list)
  x <- log(ests$Mean[, 1])
  if (cirteria == "Gini") {
    y <- ests$Gini[, 1]
    sd <- ests$Gini[, 2]
  } else {
    y <- ests$CV[, 1]
    sd <- ests$CV[, 2]
  }
  if (plot.result)
    plot(x, y, pch = 20, col = rgb(0, 0, 0, 0.4),
         xlab = "log Mean Relative Expression", ylab = paste("Estimated", criteria), 
         main = "HVG selection by DESCEND", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, font.main = 1)

  idx <- !is.na(x)

  X <- model.matrix(y[idx] ~ bs(x[idx], df = 10))
  temp <- rqss(y[idx] ~ bs(x[idx], df = 10), tau = quantile)
  y.pred <- X %*% temp$coef
  x1 <- x[idx]
  ix <- sort(x1, index.return = T)$ix

  if (plot.result)
    points(x1[ix], y.pred[ix], type = "l", col = "blue", lwd = 2)

  score <- cbind(y[idx] - y.pred, (y[idx] - y.pred)/sd[idx])
  rownames(score) <- names(x[idx])
  colnames(score) <- c("Difference", "Normalized Difference")

  sel.genes <- score[score[, 2] > threshold, ]
  sel.genes <- sel.genes[order(sel.genes[, 2], decreasing = T), ]

  sel.gene.names <- rownames(sel.genes)



  points(x[sel.gene.names], 
         y[sel.gene.names], col = rgb(1, 0, 0, 0.6), pch = 20)

  
  return(list(score.mat = score, HVG.genes = sel.gene.names))

}

