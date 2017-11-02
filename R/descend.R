#' Apply DESCEND to all the genes in the count matrix for one cell population
#'
#' Apply DESCEND to deconvolve the true expression level distribution for every geneand calculate relavant distribution measurements. Parallel computing is allowed. For deconvolution of two or more cell populations, see \code{\link{descendMultiPop}}. For model details, see \code{\link{deconvG}}. 
#'
#' @inheritParams deconvSingle
#' @param count.matrix the observed UMI count matrix. Each row is a gene and each column is a cell. The column sums are used as the input for \code{scaling.consts} when both \code{ercc.matrix} and \code{scaling.consts} are NULL.
#' @param ercc.matrix the ERCC spike-ins are used for computing the cell-specific efficiency constants as \code{scaling.consts} when \code{scaling.consts} is NULL. Each row is a spike-in genes and each column is a cell. The number and order of the columns should be the same as the number and order of the columns of \code{count.matrix}.
#' @param ercc.trueMol the true input number of molecules of the ercc spike-ins when \code{ercc.matrix} is not NULL.
#' @param n.cores the number of cores used for parallel computing. Default is 1. Used only when parallel computing is done in a single machine. For using multi-machine cores, need to assign \code{cl} explicitly. If \code{verbose} is TRUE, then a separated file is created to store the progress of each slave cores.
#' @param cl an object of class "cluster". See more details in \code{\link[parallel]{makeCluster}}
#' @param type Default is "FORK" to save memory. Change it to "PSOCK" if you are using Windows and cl is NULL. More details see \code{\link[parallel]{makeCluster}}
#' @param show.message whether show messages for the computing progresses. Default is TRUE
#'
#'
#' @return a list of DESCEND objects. The length of the list is the same as the number of genes. NA if the gene is too sparse or DESCEND fails to find a solution.
#'
#' @examples
#' \dontrun{
#' ## For a Windows machine add the argument: type = "PSOCK"
#' data(zeisel)
#' result <- runDescend(zeisel$count.matrix.small, 
#'                      scaling.consts = zeisel$library.size, n.cores = 3)
#' hvg <- findHVG(result)
#' hvg$HVG.genes
#' 
#' 
#' 
#' result1 <- runDescend(zeisel$count.matrix.small, 
#'                       zeisel$ercc.matrix, ercc.trueMol = zeisel$trueMol,
#'                       Z0 = log(zeisel$cell.size),
#'                       n.cores = 3)
#' 
#' ests <- getEstimates(result1)
#' ests$CV
#' }
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import iterators
#' @rdname runDescend
#' @export 

runDescend <- function(count.matrix,
                       ercc.matrix = NULL,
                       scaling.consts = NULL,
                       Z = NULL,
                       Z0 = NULL,
                       n.cores = 1, cl = NULL, type = "FORK",
                       do.LRT.test = F, 
                       family = c("Poisson", "Negative Binomial"),
                       NB.size = 100,
                       show.message = T,
                       verbose = T, 
                       ercc.trueMol = NULL,
                       control = list()) {

   control <- do.call("DESCEND.control", control)

  Y <- as.matrix(count.matrix)
  gene.names <- rownames(Y)
  if (sum(is.na(Y)) > 0)
    stop("Missing values in count.matrix are not allowed!")

  if (sum(Y < 0) > 0 || sum(abs(Y - round(Y))) > 0)
    stop("The data input should be raw UMI counts!!")


  if (is.null(scaling.consts)) {
    if (is.null(ercc.matrix) || is.null(ercc.trueMol)) {
      if (!is.null(ercc.matrix))
        print("The input number of molecules for the ERCC spike-in genes are not provided, library size normalization is used instead.")

      scaling.consts <- colSums(Y)
    }
    else {
      if (sum(is.na(ercc.matrix)) > 0)
        stop("Missing values in ercc.matrix are not allowed!")

      if (sum(ercc.matrix < 0) > 0 || 
          sum(abs(ercc.matrix - round(ercc.matrix))) > 0)
        stop("The ERCC input should be raw UMI counts!!")

      if (ncol(Y) != ncol(ercc.matrix))
        stop("The ERCC matrix should have the same number and order of the columns as the count matrix!")
      scaling.consts <- colSums(ercc.matrix) / sum(ercc.trueMol)
    }
  }

  if (show.message) {
    print(paste("DESCEND starts deconvolving distribution of", 
                length(gene.names), "genes!"))

    print("Estimating the time to finish this task ...")
  }

  if (!is.null(cl)) {
    n.cores <- length(parallel::clusterCall(cl, print, 1))
  }

  while(T) {
    idx <- sample(1:nrow(Y), min(3, nrow(Y)))
    tt <- system.time(temp <- apply(Y[idx, ], 1, function(v) {
                                    try(deconvSingle(v,
                                                     scaling.consts = scaling.consts,
                                                     Z = Z, Z0 = Z0,
                                                     plot.density = F,
                                                     do.LRT.test = do.LRT.test,
                                                     family = family, 
                                                     NB.size = NB.size,
                                                     control = control, verbose = F), silent = T)}))
    if (sum(sapply(temp, function(ll)class(ll) != "DESCEND")) != min(3, nrow(Y)))
      break
  }

  if (show.message)
    print(paste("Estimated total time is", 
                ceiling(tt[3]/sum(sapply(temp, 
                                         function(ll)class(ll) == "DESCEND")) * nrow(Y)/60), 
                "minutes without parallization!"))

#  if (is.null(cl)) {
    requireNamespace("doParallel")
    requireNamespace("foreach")
    requireNamespace("iterators")
    outfile <- paste("verbose_log_", as.numeric(Sys.time()), ".txt", sep = "")
    if (is.null(cl)) {
      if (verbose)
        cl <- makeCluster(n.cores, type = type, outfile = outfile)
      else
        cl <- makeCluster(n.cores, type = type)
      if (verbose && show.message)
        print(paste("log outputs are stored in file", outfile))
    }

    if (show.message)
      print(paste(n.cores, "cores are used in parallel."))

    registerDoParallel(cl)
    results <-  foreach(v = iter(Y, "row"),
                        .noexport = c("Y")) %dopar% {
      requireNamespace("descend")

      if (verbose)
        print(paste("Start computing for one gene!", Sys.time()))
    
      temp <- try(deconvSingle(as.vector(v),
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
  names(results) <- gene.names

  return(results)
}

#' Grab the value and standard deviation of the estimated DESCEND elements calculated from a list of descend objects. For examples, see  \code{\link{runDescend}}
#'
#' @param descend.list a list of descend objects computed from {\code{\link{runDescend}}}
#' @return A list of matrices where each element is for a distribution measurement or a coefficient if covariates are presented. For each distribution measurement or coefficient, the matrix contains two columns. Each row is for a gene. The first column is the estimated value and second value is the estimated standard deviation.
#'
#' @rdname getEstimates
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

#' Grab the likelihood ratio test p-values if the tests are performed from a list of descend objects 
#'
#' @inheritParams getEstimates
#'
#' @return A matrix of one column. Each row is for a distribution measurement or a coefficient if covariates are presented. 
#'
#' @rdname getPval
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
                        mm <- matrix(NA, n.test, 1)
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
#' The function finds highly variable genes based on the Gini or CV estimates from DESCEND. A quantile natural cubic spline regression of the estimated selected criteria on the log of the estimated mean for the genes is performed as the general trend. HVGs are selected as the genes whose normalized difference between its estimates and the trend is larger than a threshold. For examples, see \code{\link{runDescend}}
#'
#' @inheritParams getEstimates
#' @param criteria the cirteria of HVG finding, either "Gini" or "CV"
#' @param quantile quantile of the quantile regression. Defined as the \code{tau} parameter in the function \code{\link[quantreg]{rqss}}. A larger quantile yields a more stringent selection of HVG. Default is 0.7.
#' @param threshold threshold of the normalized difference for HVG selection. Default is 3. A larger value results in a more stringent selection
#' @param plot.result whether plot the selection results or not. Default is TRUE.
#' @param spline.df the degree of freedom of the spline functions used to fit the quantile regression curve
#'
#' @return A list of elements:
#' \item{score.mat}{A score matrix of the genes. Each row is a gene}
#' \item{HVG.genes}{A vector of the selected HVG gene names}
#'
#'
#' @import quantreg
#' @import splines
#' @import grDevices
#' @rdname findHVG
#' @export 


findHVG <- function(descend.list, 
                    criteria = c("Gini", "CV"),
                    quantile = 0.7,
                    threshold = 3,
                    plot.result = T, spline.df = 5) {
  requireNamespace("quantreg")
  requireNamespace("splines")
  criteria <- match.arg(criteria, c("Gini", "CV"))
  ests <- getEstimates(descend.list)
  x <- log(ests$Mean[, 1])
  if (criteria == "Gini") {
    y <- ests$Gini[, 1]
    sd <- ests$Gini[, 2]
  } else {
    y <- ests$CV[, 1]
    sd <- ests$CV[, 2]
  }
  if (plot.result)
    plot(x, y, pch = 20, col = grDevices::rgb(0, 0, 0, 0.4),
         xlab = "log Mean Relative Expression", ylab = paste("Estimated", criteria), 
         main = "HVG selection by DESCEND", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, font.main = 1)

  idx <- !is.na(x)

  X <- model.matrix(y[idx] ~ bs(x[idx], df = spline.df))
  temp <- rqss(y[idx] ~ bs(x[idx], df = spline.df), tau = quantile)
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
         y[sel.gene.names], col = grDevices::rgb(1, 0, 0, 0.6), pch = 20)

  
  return(list(score.mat = score, HVG.genes = sel.gene.names))

}

