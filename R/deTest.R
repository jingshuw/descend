#' DESCEND applied to two or more cell populations
#'
#' This function is used when two or more cell populations are compared with each other and is a first step for differential testing between any two of the cell populations. The true expression distribution is deconvolved for each cell population separately while \code{Z0} is scaled to have mean 0 (combining all populations) to compute a meaningful \code{Z0} adjusted nonzero fraction. For deconvolution of a single cell population, see \code{\link{runDescend}}. For model details, see \code{\link{deconvG}}. Depending on the number of cell types, number of cells and the dimension of \code{Z} and \code{Z0}, this function can take a very long time to run even on a cluster and occupy massive memory for the DESCEND results (as we have a DESCEND object for each cell type and each gene). In this scenario, we suggest users to run \code{runDescend} and save the descend result for each cell type separately, then follow the code inside this function for normalization of \code{Z0} and the calculation of \code{Z0} adjusted nonzero Fraction.
#' @inheritParams runDescend
#' @param labels a vector of factors or characters, indicating the cell popluation label of each cell. The length of \code{labels} should be the same as the number of columns of \code{count.matrix}
#' @param center.Z0 whether to center Z0 to make \code{Z0} adjusted nonzero fraction more meaningful. Default is TRUE. Set it to FALSE if \code{Z0} has already been properly centered
#'
#' @return a list with elements
#' \item{descend.list.list}{a list of DESCEND object lists. Each element is a DESCEND object list for one of the cell populations computed from \code{runDescend}.}
#' \item{model}{model parameters, including the actual \code{scaling.consts}, \code{Z}, the rescaled \code{Z0}, \code{control}, \code{family} and \code{NB.size}}
#'
#' @examples
#' \dontrun{
#' data(zeisel)
#'  set.seed(1)
#'  ## For a Windows machine add the argument: 
#'  ## type = "PSOCK" to each of the function that need parallization.
#'  result.multi <- descendMultiPop(zeisel$count.matrix.small,
#'                                  labels = zeisel$labels,
#'                                  scaling.consts = zeisel$library.size,
#'                                  Z0 = log(zeisel$cell.size), verbose = TRUE, show.message = TRUE,
#'                                  n.cores = 3)
#'  ## try 100 null genes first
#'  detest.result <- deTest(result.multi, c("endothelial-mural", "pyramidal CA1"),
#'                          zeisel$count.matrix.small, zeisel$labels,
#'                          verbose = TRUE, show.message = TRUE,
#'                          N.genes.null = 100, n.cores = 3)
#'  
#'  ## 100 null genes may not get small enough p-values
#'  detest.result <- deTest.more(result.multi, detest.result, 
#'                               c("endothelial-mural", "pyramidal CA1"),
#'                               zeisel$count.matrix.small, labels = zeisel$labels, 
#'                               N.more.genes = 200, verbose = TRUE, 
#'                               n.cores = 3)
#'  
#'  layout(matrix(1:4, nrow = 2))
#'  de.scores1 <- plotDeTest(result.multi, c("endothelial-mural", "pyramidal CA1"),
#'                          detest.result, measurement.name = "Gini", alpha = 0.05)
#'  de.scores2 <- plotDeTest(result.multi, c("endothelial-mural", "pyramidal CA1"),
#'                          detest.result, measurement.name = "Nonzero Mean", 
#'                          alpha = 0.05, log = "xy")
#'  de.scores3 <- plotDeTest(result.multi, c("endothelial-mural", "pyramidal CA1"),
#'                          detest.result, measurement.name = "Nonzero Fraction", alpha = 0.1)
#'  de.scores4 <- plotDeTest(result.multi, c("endothelial-mural", "pyramidal CA1"),
#'                          detest.result, measurement.name = "Adjusted Nonzero Fraction", alpha = 0.1)
#' }
#'
#' @rdname descendMultiPop
#' @export

descendMultiPop <- function(count.matrix,
                            labels,
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
                            center.Z0 = T, 
                            control = list()) {
   
  count.matrix <- as.matrix(count.matrix)
  if (ncol(count.matrix) != length(labels)) {
    stop("Every cell need to have a cell population label. The length of labels need to match the number of columns of count.matrix!")
  }

  labels <- as.factor(labels)

  n.pop <- length(levels(labels))

  if (!is.null(Z0)) {
    Z0 <- as.matrix(Z0)
    if (center.Z0)
      Z0 <- apply(Z0, 2, scale, scale = F)
  }

  if (!is.null(Z))
    Z <- as.matrix(Z)

  if (is.null(scaling.consts)) {
    if (is.null(ercc.matrix) || is.null(ercc.trueMol)) {
      if (!is.null(ercc.matrix))
        print("The input number of molecules for the ERCC spike-in genes are not provided, library size normalization is used instead.")

      scaling.consts <- colSums(count.matrix, na.rm = T)
    }
    else {
      if (sum(is.na(ercc.matrix)) > 0)
        stop("Missing values in ercc.matrix are not allowed!")

      if (sum(ercc.matrix < 0) > 0 || 
          sum(abs(ercc.matrix - round(ercc.matrix))) > 0)
        stop("The ERCC input should be raw UMI counts!!")

      if (ncol(count.matrix) != ncol(ercc.matrix))
        stop("The ERCC matrix should have the same number and order of the columns as the count matrix!")
      scaling.consts <- colSums(ercc.matrix) / sum(ercc.trueMol)
    }
  }

  if (!file.exists("temp"))
    dir.create("temp")
  print("Directory ./temp will be used to store DESCEND result for each cell type!")



  result.list <- lapply(levels(labels), function(str) {
                              print(paste("DESCEND is applied to cell population", 
                                          str, "!"))
                              idx <- which(labels == str)
                              scaling.consts.temp <- scaling.consts[idx]
                          
                              if (!is.null(Z)) 
                                Z.temp <- Z[idx, , drop = F]
                              else
                                Z.temp <- NULL

                              if (!is.null(Z0))
                                Z0.temp <- Z0[idx, , drop = F]
                              else
                                Z0.temp <- NULL

							  result <- runDescend(count.matrix[, idx, drop = F],
                                                ercc.matrix = NULL,
                                                scaling.consts = scaling.consts.temp,
                                                Z = Z.temp,
                                                Z0 = Z0.temp,
                                                n.cores = n.cores, cl = cl, type = type,
                                                do.LRT.test = do.LRT.test,
                                                family = family, NB.size = NB.size,
                                                verbose = verbose, show.message = show.message,
                                                control = control)
							  saveRDS(result, file = paste("temp/DESCEND_result_", str, ".rds", sep = ""))
                              rm(result)
                              gc()
                              print(paste("Results for", str, "saved."))
                              return(list())
                            })

  result.list <- lapply(levels(labels), function(str) {
                             result <- readRDS(paste("temp/DESCEND_result_", str, ".rds", sep = ""))
                             return(result)
                            }) 

  names(result.list) <- levels(labels)

  ## compute Z0 adjusted nonzero fraction
  if (!is.null(Z0)) {
    result.list <- lapply(result.list, 
                          function(ll) {
                            result <- lapply(ll, 
                                             function(lll) {
                                        if (class(lll) != "DESCEND")
                                          return(NA)
                                        temp <- lll@estimates["Z0 effect: beta0", ]
                                        p0.adj <- exp(temp[1])/(1 + exp(temp[1]))
                                        p0.adj.bias <- temp[2] * p0.adj * 
                                          (1 - p0.adj)
                                        p0.adjsd <- temp[3] * p0.adj * (1 - p0.adj)
                                        temp <- c(1- p0.adj, -p0.adj.bias, p0.adjsd,
                                                  sqrt(p0.adj.bias^2 + p0.adjsd^2))

                                        lll@estimates <- rbind(lll@estimates[1:5,,drop = F],
                                                               temp,
                                                  lll@estimates[-(1:5),, drop = F])
                                        rownames(lll@estimates)[6] <- "Z0 Adjusted Nonzero Fraction"
                                        return(lll)
                              
                              })
                            })
  }
  return(list(descend.list.list = result.list, 
              model = list(scaling.consts = scaling.consts,
              Z = Z, 
              Z0 = Z0,
              control = control, family = family, NB.size = NB.size)))
}


#' Perform differential testing of the distribution measurements between two cell populations
#'
#' Permutation test is used to get the p-values in differential testing. To keep the possible correlation structure between genes. We run \code{ceiling(N.genes.null/nrow(count.matrix))} rounds of permutation. For each round, we permute cell label once and deconvolve the true expression distribution for all the selected genes with the permuted cell labels. One may need more rounds of permutation to get smaller p-values. See \code{\link{deTest.more}}. For examples, see \code{\link{descendMultiPop}}

#'
#' @inheritParams descendMultiPop
#' @param descend.multipop.output the returned value of \code{\link{descendMultiPop}}
#' @param de.pair.names a vector of length 2 showing the two cell population names for differential testing. The names should match the values in \code{labels}
#' @param N.genes.null number of permuted genes. The larger the value is, the longer it takes and the smaller the minimum p-value can be (the minimum p-value can be as small as 0.5/N.genes.null)
#' @param alternative the alternative hypotheses for differential testing
#' @param params by default, it is descend.multipop.output$params, but users can provide their own in special cases
#' @param compare.covariates whether also perform differential testing on the coefficients of the covaraites \code{Z} and \code{Z0}. Need to set to FALSE if the number of covariates is not the same in the two cell populations
#'
#' @return a list with elements
#' \item{p.values}{a matrix of p-values calculated from permutation tests. Each row is a gene and each column is a distribution measurement of coefficient (if \code{Z} or \code{Z0} is presented)}
#' \item{scores}{the z-scores matrix in differential testing. Each row is a gene and each column is a distribution measurement of coefficient (if \code{Z} or \code{Z0} is presented)}
#' \item{perm.scores}{the permuted z-scores matrix used as the null distribution of the z-scores in each column. The number of rows equals to \code{N.genes.null}}
#'
#' @rdname deTest
#' @export

deTest <- function(descend.multipop.output,
                   de.pair.names,
                   count.matrix,
                   labels,
                   N.genes.null = 10000,
                   alternative = c("two.sided", "less", "greater"),
                   n.cores = 1, cl = NULL, type = "FORK",
                   verbose = T, show.message = T, params = NULL,
                   compare.covariates = T) {

  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))

  if (is.null(params))
    params <- descend.multipop.output$model

  n.perm.rounds <- ceiling(N.genes.null / nrow(count.matrix))
  N.genes <- nrow(count.matrix)

  idx1 <- which(labels == de.pair.names[1])
  idx2 <- which(labels == de.pair.names[2])

  perm.result1 <- list()
  perm.result2 <- list()

  print(paste("To get", N.genes.null, "permuted genes, total number of permutation rounds is", n.perm.rounds))

  for(i in 1:n.perm.rounds) {
    print(paste("Permutation to compute the null distribution, round", i))
    n <- min(N.genes, N.genes.null - (i - 1) * N.genes)
    perm.idx <- sample(c(idx1, idx2))
    temp.label <- c(rep(de.pair.names[1], length(idx1)),
                    rep(de.pair.names[2], length(idx2)))

    if (n < N.genes)
      gene.idx <- sample(1:N.genes, n)
    else
      gene.idx <- 1:N.genes

    temp <- descendMultiPop(count.matrix[gene.idx, perm.idx, drop = F],
                            labels = temp.label,
                            scaling.consts = params$scaling.consts[perm.idx],
                            Z = params$Z[perm.idx,], 
                            Z0 = params$Z0[perm.idx,],
                            n.cores = n.cores, cl = cl, type = type,
                            do.LRT.test = F,
                            family = params$family, NB.size = params$NB.size,
                            verbose = verbose, show.message = show.message, 
                            control = params$control)

    perm.result1 <- c(perm.result1, temp$descend.list.list[[1]])
    perm.result2 <- c(perm.result2, temp$descend.list.list[[2]])
  }

  perm.ests.list <- lapply(list(perm.result1, perm.result2), getEstimates)
  m <- length(perm.ests.list[[1]])
  perm.scores <- sapply(1:m, function(j) {
                            mat1 <- perm.ests.list[[1]][[j]]
                            mat2 <- perm.ests.list[[2]][[j]]

                            scores <- (mat1[ ,1] - mat2[ ,1])/
                              sqrt(mat1[ ,2]^2 + mat2[ ,2]^2)
                            return(scores)
                            })
  colnames(perm.scores) <- names(perm.ests.list[[1]])

  ori.ests.list <- lapply(descend.multipop.output$descend.list.list[as.factor(de.pair.names)], 
                          getEstimates)
  m1 <- length(ori.ests.list[[1]])
  if (!compare.covariates) {
    m1 <- 5
  }
  ori.scores <- sapply(1:m1, function(j) {
                            mat1 <- ori.ests.list[[1]][[j]]
                            mat2 <- ori.ests.list[[2]][[j]]

                            scores <- (mat1[ ,1] - mat2[ ,1])/
                              sqrt(mat1[ ,2]^2 + mat2[ ,2]^2)
                            return(scores)
                            })
  colnames(ori.scores) <- names(ori.ests.list[[1]][1:m1])

  if (m1 != m) {
    warnings(paste("Number of parameters", m, 
                   "in permutation does not match number of parameters", m1,
                   "in the original result!! Parameters matched by names and unmatched columns are removed. 
                   Be careful when using the result."))
    mm.names <- intersect(colnames(ori.scores), colnames(perm.scores))
    ori.scores <- ori.scores[, mm.names]
    perm.scores <- perm.scores[, mm.names]
    m <- length(mm.names)
  }


  pvals <- sapply(1:m, function(j) {
                  if (alternative == "two.sided") {
                    null.vals <- sort(perm.scores[, j])
                    temp <- sort(ori.scores[, j], index.return = T, na.last = T)
                  } else if (alternative == "less") {
                    null.vals <- sort(-perm.scores[, j])
                    temp <- sort(-ori.scores[, j], index.return = T, na.last = T)
                  } else {
                    null.vals <- sort(perm.scores[, j])
                    temp <- sort(ori.scores[, j], index.return = T, na.last = T)
                  }
                  k2 <- 1
                  ranks <- rep(0, length(temp$x))
                  ranks[is.na(temp$x)] <- NA
                  for(k1 in 1:sum(!is.na(temp$x))) {
                    while (temp$x[k1] > null.vals[k2] && k2 < length(null.vals)) {
                      k2 <- k2 + 1
                    }
                    if (temp$x[k1] < null.vals[k2])
                      ranks[k1] <- k2 - 1
                    else
                      ranks[k1] <- k2
                  }
                  if (alternative == "two.sided")
                    ranks[temp$x < 0 & (!is.na(temp$x))] <- length(null.vals) - 
                    ranks[temp$x < 0 & (!is.na(temp$x))]

                  pp <- (length(null.vals) - ranks + 0.25)/length(null.vals)
                  if (alternative == "two.sided")
                    pp <- pp *2

                  pp[temp$ix] <- pp
                  pp[pp > 1] <- 1

                  return(pp)
                            })

  rownames(pvals) <- rownames(ori.scores)
  colnames(pvals) <- colnames(ori.scores)

  return(list(p.values  = pvals, scores = ori.scores, perm.scores = perm.scores))
}


#' Add more permuted genes for differential testing of the distribution measurements between two cell populations
#'
#' The minimum possible p-value is limited by the number of permuted genes used. This function can add more permuted genes to the previous result to get smaller p-values. For examples, see \code{\link{descendMultiPop}}

#'
#' @inheritParams deTest
#' @param N.more.genes number of more permuted genes. The larger the value is, the longer it takes and the smaller the minimum p-value can be 
#' @param deTest.output returned values of \code{\link{deTest}} or \code{\link{deTest.more}}
#'
#' @return same as \code{\link{deTest}}
#' @rdname deTest.more
#' @export

deTest.more <- function(descend.multipop.output,
                        deTest.output,
                        de.pair.names,
                        count.matrix,
                        labels,
                        N.more.genes = 10000,
                        alternative = c("two.sided", "less", "greater"),
                        n.cores = 1, cl = NULL, type = "FORK",  
                        verbose = T, show.message = T, params = NULL,
                        compare.covariates = T) {

  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))

  params <- descend.multipop.output$model
  count.matrix <- as.matrix(count.matrix)

  n.perm.rounds <- ceiling(N.more.genes / nrow(count.matrix))
  N.genes <- nrow(count.matrix)

  idx1 <- which(labels == de.pair.names[1])
  idx2 <- which(labels == de.pair.names[2])

  perm.result1 <- list()
  perm.result2 <- list()

  print(paste("To get", N.more.genes, 
              "more permuted genes, total number of permutation rounds is", n.perm.rounds))

  for(i in 1:n.perm.rounds) {
    print(paste("Permutation to compute the null distribution, round", i))
    n <- min(N.genes, N.more.genes - (i - 1) * N.genes)
    perm.idx <- sample(c(idx1, idx2))
    temp.label <- c(rep(de.pair.names[1], length(idx1)),
                    rep(de.pair.names[2], length(idx2)))

    if (n < N.genes)
      gene.idx <- sample(1:N.genes, n)
    else
      gene.idx <- 1:N.genes

    temp <- descendMultiPop(count.matrix[gene.idx, perm.idx, drop = F],
                            labels = temp.label,
                            scaling.consts = params$scaling.consts[perm.idx],
                            Z = params$Z[perm.idx,], 
                            Z0 = params$Z0[perm.idx,],
                            n.cores = n.cores, cl = cl, type = type,
                            do.LRT.test = F,
                            family = params$family, NB.size = params$NB.size,
                            verbose = verbose, show.message = show.message,
                            control = params$control)

    perm.result1 <- c(perm.result1, temp$descend.list.list[[1]])
    perm.result2 <- c(perm.result2, temp$descend.list.list[[2]])
  }

  perm.ests.list <- lapply(list(perm.result1, perm.result2), getEstimates)
  m <- length(perm.ests.list[[1]])
  perm.scores <- sapply(1:m, function(j) {
                            mat1 <- perm.ests.list[[1]][[j]]
                            mat2 <- perm.ests.list[[2]][[j]]

                            scores <- (mat1[ ,1] - mat2[ ,1])/
                              sqrt(mat1[ ,2]^2 + mat2[ ,2]^2)
                            return(scores)
                            })
  colnames(perm.scores) <- names(perm.ests.list[[1]])
  perm.scores <- rbind(perm.scores, deTest.output$perm.scores)

  ori.ests.list <- lapply(descend.multipop.output$descend.list.list[as.factor(de.pair.names)], 
                          getEstimates)
  m1 <- length(ori.ests.list[[1]])
  if (!compare.covariates) {
    m1 <- 5
  }
  ori.scores <- sapply(1:m1, function(j) {
                            mat1 <- ori.ests.list[[1]][[j]]
                            mat2 <- ori.ests.list[[2]][[j]]

                            scores <- (mat1[ ,1] - mat2[ ,1])/
                              sqrt(mat1[ ,2]^2 + mat2[ ,2]^2)
                            return(scores)
                            })
  colnames(ori.scores) <- names(ori.ests.list[[1]][1:m1])

  if (m1 != m) {
    warnings(paste("Number of parameters", m, 
                   "in permutation does not match number of parameters", m1,
                   "in the original result!! Parameters matched by names and unmatched columns are removed. 
                   Be careful when using the result."))
    mm.names <- intersect(colnames(ori.scores), colnames(perm.scores))
    ori.scores <- ori.scores[, mm.names]
    perm.scores <- perm.scores[, mm.names]
    m <- length(mm.names)
  }

  pvals <- sapply(1:m, function(j) {
                  if (alternative == "two.sided") {
                    null.vals <- sort(perm.scores[, j])
                    temp <- sort(ori.scores[, j], index.return = T, na.last = T)
                  } else if (alternative == "less") {
                    null.vals <- sort(-perm.scores[, j])
                    temp <- sort(-ori.scores[, j], index.return = T, na.last = T)
                  } else {
                    null.vals <- sort(perm.scores[, j])
                    temp <- sort(ori.scores[, j], index.return = T, na.last = T)
                  }
                  k2 <- 1
                  ranks <- rep(0, length(temp$x))
                  ranks[is.na(temp$x)] <- NA
                  for(k1 in 1:sum(!is.na(temp$x))) {
                    while (temp$x[k1] > null.vals[k2] && k2 < length(null.vals)) {
                      k2 <- k2 + 1
                    }
                    if (temp$x[k1] < null.vals[k2])
                      ranks[k1] <- k2 - 1
                    else
                      ranks[k1] <- k2
                  }
                  if (alternative == "two.sided")
                    ranks[temp$x < 0 & (!is.na(temp$x))] <- length(null.vals) - 
                    ranks[temp$x < 0 & (!is.na(temp$x))]

                  pp <- (length(null.vals) - ranks + 0.5)/length(null.vals)
                  if (alternative == "two.sided")
                    pp <- pp *2

                  pp[temp$ix] <- pp
                  pp[pp > 1] <- 1

                  return(pp)
                            })

  rownames(pvals) <- rownames(ori.scores)
  colnames(pvals) <- colnames(ori.scores)

  return(list(p.values  = pvals, scores = ori.scores, perm.scores = perm.scores))
}

#' Plot the differential testing result of one distribution measurement or coefficient
#'
#' Given a significance level and the measurement name, the significiant genes are visualized as red dots. For examples, see \code{\link{descendMultiPop}}

#'
#' @inheritParams deTest.more
#' @param alpha the nomial significance level
#' @param method default is BH. More details see \code{\link[stats]{p.adjust}}
#' @param measurement.name one of the colnames of \code{deTest.ouput$p.values}
#' @param log control log transformation of x/y axes, see \code{\link[graphics]{plot}}
#'
#' @return a matrix of five columns. Each row is a gene.
#' @rdname plotDeTest
#' @export

plotDeTest <- function(descend.multipop.output,
                       de.pair.names,
                       deTest.output,
                       alpha = 0.05,
                       method = "BH",
                       measurement.name = "Gini", log = "") {

  if (!measurement.name %in% colnames(deTest.output$p.values)) {
    stop("Choose measurement.name as one of the colnames of deTest.output$p.values!")
  }

  output <- descend.multipop.output$descend.list.list[as.factor(de.pair.names)]

  ests <- sapply(output, function(ll)getEstimates(ll)[[measurement.name]][, 1])

  pvals <- deTest.output$p.values[, measurement.name]

  p.adj <- p.adjust(pvals, method = method)
  sig.idx <- p.adj < alpha

  plot(ests, pch = 20, xlab = de.pair.names[1], ylab = de.pair.names[2],
        main = measurement.name, log = log)
  points(ests[sig.idx, ], col = "red", pch = 20)
  abline(a = 0, b= 1, col = "blue")

  stats <- cbind(ests, deTest.output$scores[, measurement.name], pvals, p.adj)
  stats <- stats[order(stats[, 4]), ]

  colnames(stats) <- c(paste("population", de.pair.names), "z-score", "p-value", "Adjusted p-value")
  return(stats)
  
}
