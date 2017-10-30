#' DESCEND applied to two or more cell populations
#'
#' This function is used when two or more cell populations are compared with each other and is a first step for differential testing between any two of the cell populations. The true expression distribution is deconvolved for each cell population separately while \code{Z0} is scaled to have mean 0 (combining all populations) to compute a meaningful \code{Z0} adjusted active fraction.
#'
#' @inheritParams runDescend
#' @param labels a vector of factors or characters, indicating the cell popluation label of each cell. The length of \code{labels} should be the same as the number of columns of \code{count.matrix}
#' @param center.Z0 whether to center Z0 to make \code{Z0} adjusted active fraction more meaningful. Default is TRUE 
#'
#' @return a list with elements
#' \item{descend.list.list}{a list of DESCEND object lists. Each element is a DESCEND object list for one of the cell populations computed from \code{runDescend}.}
#' \item{model}{model parameters, including the actual \code{scaling.consts}, \code{Z}, the rescaled \code{Z0}, \code{control}, \code{family} and \code{NB.size}}
#'
#' @export

descendMultiPop <- function(count.matrix,
                            labels,
                            ercc.matrix = NULL,
                            scaling.consts = NULL,
                            Z = NULL,
                            Z0 = NULL,
                            n.cores = 1, cl = NULL,
                            do.LRT.test = F, 
                            family = c("Poisson", "Negative Binomial"),
                            NB.size = 100,
                            verbose = T, 
                            ercc.trueMol = NULL,
                            center.Z0 = T, 
                            control = list()) {
    
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


  result.list <- lapply(levels(labels), function(str) {
                              print(paste("DESCEND is applied to cell population", 
                                          str, "!"))
                              idx <- which(labels == str)
                              scaling.consts.temp <- scaling.consts[idx]
                          
                              if (!is.null(Z)) 
                                Z.temp <- Z[idx, ]
                              else
                                Z.temp <- NULL

                              if (!is.null(Z0))
                                Z0.temp <- Z0[idx, ]
                              else
                                Z0.temp <- NULL

                              return(runDescend(count.matrix[, idx],
                                                ercc.matrix = NULL,
                                                scaling.consts = scaling.consts.temp,
                                                Z = Z.temp,
                                                Z0 = Z0.temp,
                                                n.cores = n.cores, cl = cl, 
                                                do.LRT.test = do.LRT.test,
                                                family = family, NB.size = NB.size,
                                                verbose = verbose, 
                                                control = control))
                            })
  names(result.list) <- levels(labels)

  ## compute Z0 adjusted active fraction
  if (!is.null(Z0)) {
    result.list <- lapply(result.list, 
                          function(ll) {
                            result <- lapply(ll, 
                                             function(lll) {
                                        temp <- lll@estimates["Z0 effect: beta0", ]
                                        p0.adj <- exp(temp[1])/(1 + exp(temp[1]))
                                        p0.adj.bias <- temp[2] * p0.adj * 
                                          (1 - p0.adj)
                                        p0.adjsd <- temp[3] * p0.adj * (1 - p0.adj)
                                        temp <- c(p0.adj, p0.adj.bias, p0.adjsd)

                                        lll@estimates <- rbind(lll@estimates[1:5, ],
                                                               temp,
                                                  lll@estimates[-(1:5),, drop = F])
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
#' Permutation test is used to get the p-values in differential testing. To keep the possible correlation structure between genes. We run \code{ceiling(N.genes.null/nrow(count.matrix))} rounds of permutation. For each round, we permute cell label once and deconvolve the true expression distribution for all the selected genes with the permuted cell labels. One may need more rounds of permutation to get smaller p-values. See \code{\link{deTest.more}}.
#'
#' @inheritParams descendMultiPop
#' @param descend.multipop.output the returned value of \code{\link{descendMultiPop}}
#' @param de.pair.names a vector of length 2 showing the two cell population names for differential testing. The names should match the values in \code{labels}
#' @param N.genes.null number of permuted genes. The larger the value is, the longer it takes and the smaller the minimum p-value can be (the minimum p-value can be as small as 0.5/N.genes.null)
#' @param alternative the alternative hypotheses for differential testing
#'
#' @return a list with elements
#' \item{p.values}{a matrix of p-values calculated from permutation tests. Each row is a gene and each column is a distribution measurement of coefficient (if \code{Z} or \code{Z0} is presented)}
#' \item{scores}{the z-scores matrix in differential testing. Each row is a gene and each column is a distribution measurement of coefficient (if \code{Z} or \code{Z0} is presented)}
#' \item{perm.scores}{the permuted z-scores matrix used as the null distribution of the z-scores in each column. The number of rows equals to \code{N.genes.null}}
#'
#' @export

deTest <- function(descend.multipop.output,
                   de.pair.names,
                   count.matrix,
                   labels,
                   N.genes.null = 10000,
                   alternative = c("two.sided", "less", "greater"),
                   n.cores = 1, cl = NULL,
                   verbose = T) {

  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))

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

    temp <- descendMultiPop(count.matrix[gene.idx, perm.idx],
                            labels = temp.label,
                            scaling.consts = params$scaling.consts[perm.idx],
                            Z = params$Z[perm.idx,], 
                            Z0 = params$Z0[perm.idx,],
                            n.cores = n.cores, cl = cl,
                            do.LRT.test = F,
                            family = params$family, NB.size = params$NB.size,
                            verbose = verbose, control = params$control)

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
  ori.scores <- sapply(1:m, function(j) {
                            mat1 <- ori.ests.list[[1]][[j]]
                            mat2 <- ori.ests.list[[2]][[j]]

                            scores <- (mat1[ ,1] - mat2[ ,1])/
                              sqrt(mat1[ ,2]^2 + mat2[ ,2]^2)
                            return(scores)
                            })
  colnames(ori.scores) <- names(ori.ests.list[[1]])

  pvals <- sapply(1:m, function(j) {
                  if (alternative == "two.sided") {
                    null.vals <- sort(abs(perm.scores[, j]))
                    temp <- sort(abs(ori.scores[, j]), index.return = T, na.last = T)
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
                  pp <- (length(null.vals) - ranks + 0.5)/length(null.vals)
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
#' The minimum possible p-value is limited by the number of permuted genes used. This function can add more permuted genes to the previous result to get smaller p-values.
#'
#' @inheritParams deTest
#' @param N.more.genes number of more permuted genes. The larger the value is, the longer it takes and the smaller the minimum p-value can be 
#' @param deTest.output returned values of \code{\link{deTest}} or \code{\link{deTest.more}}
#'
#' @return same as \code{\link{deTest}}
#'
#' @export

deTest.more <- function(descend.multipop.output,
                        deTest.output,
                        de.pair.names,
                        count.matrix,
                        labels,
                        N.more.genes = 10000,
                        alternative = c("two.sided", "less", "greater"),
                        n.cores = 1, cl = NULL,
                        verbose = T) {

  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))

  params <- descend.multipop.output$model

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

    temp <- descendMultiPop(count.matrix[gene.idx, perm.idx],
                            labels = temp.label,
                            scaling.consts = params$scaling.consts[perm.idx],
                            Z = params$Z[perm.idx,], 
                            Z0 = params$Z0[perm.idx,],
                            n.cores = n.cores, cl = cl,
                            do.LRT.test = F,
                            family = params$family, NB.size = params$NB.size,
                            verbose = verbose, control = params$control)

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
  ori.scores <- sapply(1:m, function(j) {
                            mat1 <- ori.ests.list[[1]][[j]]
                            mat2 <- ori.ests.list[[2]][[j]]

                            scores <- (mat1[ ,1] - mat2[ ,1])/
                              sqrt(mat1[ ,2]^2 + mat2[ ,2]^2)
                            return(scores)
                            })
  colnames(ori.scores) <- names(ori.ests.list[[1]])

  pvals <- sapply(1:m, function(j) {
                  if (alternative == "two.sided") {
                    null.vals <- sort(abs(perm.scores[, j]))
                    temp <- sort(abs(ori.scores[, j]), index.return = T, na.last = T)
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
                  pp <- (length(null.vals) - ranks + 0.5)/length(null.vals)
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
#' Given a significance level and the measurement name, the significiant genes are visualized as red dots. 
#'
#' @inheritParams deTest.more
#' @param alpha the nomial significance level
#' @param method default is BH. More details see \code{\link[stats]{p.adjust}}
#' @param measurement.name one of the colnames of \code{deTest.ouput$p.values}
#'
#' @return a matrix of five columns. Each row is a gene.
#'
#' @export

plotDeTest <- function(descend.multipop.output,
                       de.pair.names,
                       deTest.output,
                       alpha = 0.05,
                       method = "BH",
                       measurement.name = "Gini") {

  if (!measurement.name %in% colnames(deTest.output$p.values)) {
    stop("Choose measurement.name as one of the colnames of deTest.output$p.values!")
  }

  output <- descend.multipop.output$descend.list.list[as.factor(de.pair.names)]

  ests <- sapply(output, function(ll)getEstimates(ll)[[measurement.name]][, 1])

  pvals <- deTest.output$p.values[, measurement.name]

  p.adj <- p.adjust(pvals, method = method)
  sig.idx <- p.adj < alpha

  plot(ests, pch = 20, xlab = de.pair.names[1], ylab = de.pair.names[2],
        main = measurement.name)
  points(ests[sig.idx, ], col = "red", pch = 20)
  abline(a = 0, b= 1, col = "blue")

  stats <- cbind(ests, deTest.output$scores[, measurement.name], pvals, p.adj)
  stats <- stats[order(stats[, 4]), ]

  colnames(stats) <- c(paste("population", de.pair.names), "z-score", "p-value", "Adjusted p-value")
  return(stats)
  
}
