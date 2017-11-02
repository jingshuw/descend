#' Deconvolve the true gene expression distribution  of a single gene
#'
#' The deconvolution is computed by using the function {\code{\link{deconvG}}}. This function can automatically discretize the underlying distribution and find the proper tuning parameter \code{c0} of the penalty term. Besides, it computes the estimates and standard deviations of five distribution based statistics (active fraction, active intensity, mean, CV and gini coefficient), as well as the estimated coefficients of the covariates on active intensity (Z) and active fraction (Z0), and store them in a DESCEND object.
#'
#' @param y a vector of observed counts across cells for a single gene
#' @param scaling.consts a vector of cell specific scaling constants, either the cell efficiency or the library size
#' @inheritParams deconvG
#' @param plot.density whether plot the density curve of the deconvolved the distribution or not. The zero inflation part has been smoothed into the density curve for visualization. Default is True.
#' @param do.LRT.test whether do LRT test on the coefficients and active fraction or not. Default is True
#' @param verbose verbose the estimation and testing procedures or not. Default is True.
#' @param control settings see {\code{\link{DESCEND.control}}}
#'
#' @return a DESCEND object. See also \code{\link{DESCEND}}
#' @examples
#'
#' X <- rpois(1000, 0.2 * 3)
#' Z <- rnorm(1000)
#' result <- deconvSingle(X, Z = Z, Z0 = Z, scaling.consts = rep(0.2, 1000), do.LRT.test = TRUE)
#' result@estimates
#
#' @import graphics
#' @rdname deconvSingle
#' @export


deconvSingle <- function(y, 
                         scaling.consts = NULL, 
                         Z = NULL,
                         Z0 = NULL,
                         plot.density = T,
                         do.LRT.test = T,
                         family = c("Poisson", "Negative Binomial"),
                         NB.size = 100,
                         verbose = T, control = list()) {
  family <- match.arg(family)
  control <- do.call("DESCEND.control", control)

  if (mean(y == 0) > control$max.sparse[1] || sum(y != 0) < control$max.sparse[2])
    stop("Too sparse data!")

  if (is.null(scaling.consts)) {
    offset <- rep(0, length(y))
  } else {
    if (min(scaling.consts) <= 0)
      stop("Cell specific scaling constants should be positive!")
    offset <- log(scaling.consts)
  }


  temp <- y/exp(offset)
  lam.max <- quantile(temp[temp > 0], probs = control$max.quantile)
  if (control$only.integer) {
    if (lam.max <= control$n.points)
      tau <- log(1:ceiling(lam.max))
    else 
      tau <- c(log(1:5), log(seq(5, lam.max, length.out = control$n.points - 5 + 1)[-1]))
  } else {
    lam.min <- lam.max/(2 * control$n.points - 1)
    tau <- log(seq(lam.min, lam.max, length.out = control$n.points))
  }

  ## rescale Z to guarantee that it does not affect the range of theta, will scale back later
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    Z.means <- log(colMeans(exp(Z)))
    Z <- t(t(Z) - Z.means)
  }

  if (!is.null(Z0)) 
    Z0 <- as.matrix(Z0)

  result <- list()
  value <- Inf
  c0.start <- control$c0.start

  if (!control$zeroInflate)
    p <- control$pDegree
  else {
    p <- control$pDegree + 1
    if (!is.null(Z0))
      p <- p + ncol(Z0)
  }

  for (i in 1:control$nStart) {
    if (verbose)
      print(paste("Compute penalized MLE, the ", i, "th time", sep = ""))
    aStart <- rnorm(p, control$aStart, control$start.sd)
    if (!is.null(Z0))
      bStart <- rnorm(ncol(Z0), control$bStart, control$start.sd)
    if (!is.null(Z))
      gStart <- rnorm(ncol(Z), control$gStart, control$start.sd)

    if (verbose)
      print(paste("c0 is ", c0.start))
    temp <- try(deconvG(tau, y, offset = offset,
                        Z = Z, Z0 = Z0,
                        family = family, NB.size = NB.size,
                        zeroInflate = control$zeroInflate,
                        c0 = c0.start,
                        pDegree = control$pDegree,
                        aStart = aStart, bStart = bStart, gStart = gStart))

    if (verbose)
      print(paste("Relative fake information is:", round(temp$S * 100, 2), "%", sep = ""))

    if (class(temp) == "try-error" || is.na(temp$S))
      return(NA)
   
    iter <- 0
    while (1) {
      if (iter == control$max.c0.iter) {
        temp <- list(value = Inf)
        if (verbose)
          warning("Fail to find a proper c0!")
        break
      }

      if (temp$S < control$rel.info.range[2] & temp$S > control$rel.info.range[1])
        break
      scale.c0 <- control$rel.info.guide/temp$S
      scale.c0 <- min(scale.c0, 100)
      scale.c0 <- max(scale.c0, 0.01)
      c0.start <- max(scale.c0 * c0.start, control$c0.min)
      if (verbose)
        print(paste("c0 is ", c0.start))
      temp <- try(deconvG(tau, y, offset = offset,
                            Z = Z, Z0 = Z0,
                            family = family, NB.size = NB.size,
                            zeroInflate = control$zeroInflate,
                            c0 = c0.start,
                            pDegree = control$pDegree,
                            aStart = aStart, bStart = bStart, gStart = gStart))
      if (verbose)
        print(paste("Relative fake information is:", round(temp$S * 100, 2), "%", sep = ""))
      if (class(temp) == "try-error")
        return(NA)
      iter <- iter + 1
    }
    if (temp$value < value) {
      value <- temp$value
      result <- temp
    }
  }

  if (!is.finite(value))
    return(NA)

  ## collect estimation of statistics
  theta <- result$stats$mat.dist[, 1]
  g <- result$stats$mat.dist[, 2]
  bias.g <- result$stats$mat.dist[, 6]

  mu <- c(sum(theta * g), sum(theta * bias.g),
          sqrt(t(theta) %*% result$cov.g %*% theta))
  sd <- sqrt(sum(theta^2 * g) - mu[1]^2)
  kappa <- sum(theta^2 * g)
  CV2.devg <- theta^2  / mu[1]^2 - 2 * kappa * theta / mu[1]^3
  CV <- sd/mu[1]
  CV.devg <- CV2.devg / (2 * CV)
  CV <- c(CV, sum(CV.devg * bias.g), sqrt(t(CV.devg) %*% result$cov.g %*% CV.devg))
  gamma.devg <- sapply(1:length(g), function(i)sum(abs(theta[i] - theta) * g))
  gamma <- sum(g * gamma.devg) /2
  gini.devg <- gamma.devg / mu[1] - gamma * theta/mu[1]^2
  gini <- c(gamma/mu[1], sum(gini.devg * bias.g), sqrt(t(gini.devg) %*% result$cov.g %*% gini.devg))
  if (control$zeroInflate) {
    p0 <- c(result$stats$mat.dist[1, c(2, 6, 3)])
    if (!is.null(Z)) {
      coefs <- result$stats$mat.coef[1:ncol(Z), , drop = F]
      idx <- (nrow(result$cov) - ncol(Z) + 1):nrow(result$cov)
      cov.coef <- result$cov[idx, idx, drop = F]
      mu.pos.val <- sum(theta[-1] * g[-1]/(1 - g[1]))
      mu.pos.adj <- exp(-sum(coefs[, 1] * Z.means))
      mu.pos.dev <- theta[-1]/(1-g[1]) - mu[1]/(1-g[1])^2
      mu.pos.adj.dev <- -mu.pos.adj * Z.means
      mu <- c(mu[1] * mu.pos.adj, 
              mu.pos.adj * mu[2] + mu[1] * sum(mu.pos.adj.dev * coefs[, 2]),
              sqrt(mu.pos.adj^2 * t(theta) %*% result$cov.g %*% theta+ 
                   2 * mu.pos.adj * mu[1] * t(theta) %*% 
                   result$cov.g.gamma %*% mu.pos.adj.dev + 
                   mu[1]^2 * t(mu.pos.adj.dev) %*% cov.coef %*% mu.pos.adj.dev))
      mu.pos <- c(mu.pos.val * mu.pos.adj, 
                  mu.pos.adj * sum(mu.pos.dev * bias.g[-1]) + mu.pos.val * sum(mu.pos.adj.dev * coefs[, 2]),
                  sqrt(mu.pos.adj^2 * t(mu.pos.dev) %*% result$cov.g[-1, -1] %*% mu.pos.dev+ 
                       2 * mu.pos.adj * mu.pos.val * t(mu.pos.dev) %*% 
                       result$cov.g.gamma[-1, , drop = F] %*% mu.pos.adj.dev + 
                       mu.pos.val^2 * t(mu.pos.adj.dev) %*% cov.coef %*% mu.pos.adj.dev))
    } else {
      mu.pos.val <- sum(theta[-1] * g[-1]/(1 - g[1]))
      mu.pos.dev <- theta[-1]/(1-g[1]) - mu[1]/(1-g[1])^2
      mu.pos <- c(mu.pos.val, 
                  sum(mu.pos.dev * bias.g[-1]), 
                  sqrt(t(mu.pos.dev) %*% result$cov.g[-1, -1] %*% mu.pos.dev))

    }
    sd.pos <- sqrt(sum(theta[-1]^2 * g[-1]/(1 - g[1])) - mu.pos[1]^2)
  } else {
    p0 <- rep(0, 3)
    mu.pos <- mu
  }

  one.minus.p0 <- c(1 - p0[1], -p0[2], p0[3])

  estimates <- rbind(one.minus.p0, mu.pos, mu, CV, 
                      gini)
  colnames(estimates) <- c("est", "bias", "sd")
  rownames(estimates) <- c("Active Fraction", "Active Intensity", "Mean", "CV", "Gini")
  if (!is.null(result$stats$mat.coef))
    estimates <- rbind(estimates, result$stats$mat.coef)

  if (!control$only.integer) {
    if (control$zeroInflate) {
      theta <- theta[-1]
      g1 <- g
      if (control$dense.add.0) {
        g1[2] <- g1[1] + g1[2]
        g1 <- g1[-1]
      } else {
        g1 <- g1[-1]/sum(g1[-1])
      }
    } else
      g1 <- g
    density.points <- cbind(theta = theta, density = g1/(theta[2] - theta[1]))
  } else
    density.points <- cbind(theta, probability = g)

  if (plot.density) {
    if (!control$only.integer)
      plot(density.points, pch = 20, type = "l", ylab = "density", 
           main = "Deconvolved Distribution")
    else
      plot(density.points, pch = 20, ylab = "density", 
           main = "Deconvolved Distribution")
  }

  if (do.LRT.test) {   
    LRT.Z.select <- control$LRT.Z.select
    LRT.Z0.select <- control$LRT.Z0.select
    LRT.Z.values <- control$LRT.Z.values
    condition.list <- list(list(control$zeroInflate, Z, Z0, p, offset))
    test.name <-  c("Full model")
    if (control$zeroInflate) {
      if (!is.null(Z0))
        p1 <- p - 1 - ncol(Z0)
      else
        p1 <- p-1
      condition.list <- c(condition.list, p0 = list(list(!control$zeroInflate, Z, NULL, p1, offset)))
      names(condition.list)[length(condition.list)] <- "Active Fraction = 1"
      test.name <- c(test.name, "Active fraction = 1")
    }
    if (!is.null(Z)) {
      if (length(LRT.Z.select) == 1 || length(LRT.Z.select) == ncol(Z) &&
          (length(LRT.Z.values) == 1 || length(LRT.Z.values) == ncol(Z))) {
        if (length(LRT.Z.select) == 1)
          LRT.Z.select <- rep(LRT.Z.select, ncol(Z))
        if (length(LRT.Z.values) == 1)
          LRT.Z.values <- rep(LRT.Z.values, ncol(Z))
        if (ncol(Z) == 1 && LRT.Z.select) {
          offset1 <- offset + LRT.Z.values * Z
          condition.list <- c(condition.list, Z = list(list(control$zeroInflate, NULL, Z0, p, offset1)))
          names(condition.list)[length(condition.list)] <- paste("Effect of Z: gamma = ",
                                                                 LRT.Z.values, sep = "")
          test.name <- c(test.name, paste("Effect of Z: gamma = ", LRT.Z.values, sep = ""))
        } else
          for(i in 1:ncol(Z)) {
            if (LRT.Z.select[i]) {
              offset1 <- offset + LRT.Z.values[i] * Z[, i]
              condition.list <- c(condition.list, 
                                  list(list(control$zeroInflate, Z[, -i, drop = F], Z0, p, offset1)))
              names(condition.list)[length(condition.list)] <- paste("Effect of Z: gamma", i, 
                                                                      " = ", LRT.Z.values[i], sep = "")
              test.name <- c(test.name, paste("Effect of Z: gamma", i, " = ", LRT.Z.values[i], sep = ""))
            }
          }
      } else 
        warning("length of LRT.Z.select should be either 1 or the same as 
                ncol(Z), LRT test on Z not performed!")
    }
    if (!is.null(Z0) && control$zeroInflate) {
      if (length(LRT.Z0.select) == 1 || length(LRT.Z0.select) == ncol(Z0)) {
        if (length(LRT.Z0.select) == 1)
        LRT.Z0.select <- rep(LRT.Z0.select, ncol(Z0))
        if (ncol(Z0) == 1 && LRT.Z0.select) {
          condition.list <- c(condition.list, Z0 = list(list(control$zeroInflate, Z, NULL, p - 1, offset)))
          names(condition.list)[length(condition.list)] <- "Effect of Z0: beta = 0"
          test.name <- c(test.name, "Effect of Z0: beta = 0")
        } else
          for(i in 1:ncol(Z0)) {
            if (LRT.Z0.select[i]) {
              condition.list <- c(condition.list, 
                                  list(list(control$zeroInflate, Z, Z0[, -i, drop = F], p-1, offset)))
              names(condition.list)[length(condition.list)] <- paste("Effect of Z0: beta", i, "= 0", sep = "")
              test.name <- c(test.name, paste("Effect of Z0: beta", i, " = 0", sep = ""))
            }
          }
      } else 
        warning("length of LRT.Z.select should be either 1 or the same as 
                ncol(Z), LRT test on Z0 not performed!")
    }
 #   print(condition.list)


    if (length(condition.list) > 1) {
      if (verbose)
        print(paste("Compute LRT tests, in total", length(condition.list) - 1, "tests"))
      lk <- lapply(1:length(condition.list), function(k) {
                   if (verbose)
                     print(paste("Compute LRT test, ", test.name[k], sep = ""))
                   ll <- condition.list[[k]]
                   val <- Inf
                   for (i in 1:control$nStart.lrt) {
                     aStart <- rnorm(ll[[4]], control$aStart, control$start.sd)
                     if (!is.null(ll[[3]]))
                       bStart <- rnorm(ncol(ll[[3]]), control$bStart, control$start.sd)
                     if (!is.null(ll[[2]]))
                       gStart <- rnorm(ncol(ll[[2]]), control$gStart, control$start.sd)

                     temp <- try(deconvG(tau, y, offset = ll[[5]],
                                           Z = ll[[2]], Z0 = ll[[3]],
                                           family = family, NB.size = NB.size,
                                           zeroInflate = ll[[1]],
                                           c0 = 0,
                                           pDegree = control$pDegree,
                                           aStart = aStart, bStart = bStart, gStart = gStart,
                                           only.value = T))
                     if (class(temp) == "try-error")
                       temp <- Inf
                     if (temp < val)
                       val <- temp
                   }
                   df <- max(1, p-ll[[4]])
                   if (k == 1)
                     df <- NA
                   if (verbose)
                     print(paste("neg loglikehood is ", as.numeric(val), 
                                 "df is ", df))
                   return(c(as.numeric(val), max(1, p - ll[[4]])))
    })
      if (is.finite(lk[[1]][1])) {
        pval <- sapply(lk[-1], function(vec) {
                       if (is.finite(vec[1]))
                          return(pchisq(2 * (vec[1] - lk[[1]][1]), 
                                        vec[2], lower.tail = F))
                        return(NA)
                      })
        pval <- as.matrix(pval, ncol = 1)
        colnames(pval) <- "P.value"
        rownames(pval) <- names(condition.list)[-1]
      } else
        pval <- NA
    } else
      pval <- NA
  } else {
    pval <- NA
  }

  estimates <- cbind(estimates, DESCEND.sd = sqrt(estimates[, 2]^2 + estimates[, 3]^2))

  result <- new("DESCEND", 
                distribution = result$stats$mat.dist,
                estimates = estimates,
                density.points = density.points)
  if (!is.na(pval))
    result@pval <- pval

  return(result)
}

#' The control parameters of the function \code{\link{runDescend}} and \code{\link{deconvSingle}}
#'
#' @param n.points number of discritized points of the underlying true expression distribution. Default is 50
#' @param nStart number of random starts for the optimization problem (as it is non-convex) to find the global minimum. Default is 2
#' @param nStart.lrt number of random starts for the unpenalized optimization problem for likelihood ratio testing
#' @param max.quantile the maximum quantile of the non-zero observed counts used for finding the range of the underlying true expression distribution. Default is 0.98
#' @param max.sparse a vector of 2 indicating the maximum sparsity allowed for a gene to be computed. The first element is the fraction of zero-counts allowed, the second element is the minimum number of non-zero counts. Both criteria should be satisfied. Default is (0.95, 25). For studying active fraction, one should increase the threshold to get estimates with acceptable accuracy.
#' @param LRT.Z.select a vector of length 1 or the number of columns of \code{Z} indicating for which column of \code{Z} the coefficients are tested against the corresponding value in \code{LRT.Z.values} using LRT. Default is TRUE, meaning that all columns are tested when LRT tests are performed.
#' @param LRT.Z0.select a vector of length 1 or the number of columns of \code{Z0} indicating for which column of \code{Z0} the coefficients are tested against 0 using LRT. Default is TRUE, meaning that all columns are tested when LRT tests are performed.
#' @param LRT.Z.values a vector of length 1 or the number of columns of \code{Z} showing the values that LRT tests on coefficients of \code{Z} are tested against. Default value is 0, meanings that all tests are tested against 0.
#' @param zeroInflate whether to include the zero inflation part to the deconvolved distribution. Default is TRUE
#' @param dense.add.0 whether smooth the density at 0 into \code{density.points} in the output DESCEND object
#' @param only.integer whether set the discrete points to be integers. Default is FALSE
#' @param rel.info.range the relative information range allowed for finding the optimal tuning parameter \code{c0}
#' @param rel.info.guide one parameter inside the \code{rel.info.range} controling for the searching process of \code{c0}
#' @param c0.start the starting value of \code{c0}. Default is 1
#' @param aStart the starting values of the spline coefficients of the deconvolved distribution
#' @param bStart the starting values of the coefficients of Z0
#' @param gStart the starting values of the coefficients of Z
#' @param start.sd standard deviation of the random starting values when nStart > 1
#' @param penalty.Z0 whether add penalty to the coefficients of Z0 or not in optimization. Default is TRUE
#' @param pDegree degree of the spline bases. Default is 5
#' @param max.c0.iter maximum iteration allowed to find the optimal \code{c0}
#' @param c0.min minimum \code{c0} allowed to avoid degeneration
#'
#' @export

DESCEND.control <- function(n.points = 50,
                            nStart = 2,
                            nStart.lrt = 2,
                            max.quantile = 0.95,
                            max.sparse = c(0.99, 20),
                            LRT.Z.select = T,
                            LRT.Z0.select = T,
                            LRT.Z.values = 0,
                            zeroInflate = T,
                            dense.add.0 = T,
                            only.integer = F,
                            rel.info.range = c(0.0005, 0.01),
                            rel.info.guide = 0.005,
                            c0.start = 1,
                            aStart = 1,
                            bStart = 0,
                            gStart = 0, 
                            start.sd = 0.2, 
                            penalty.Z0 = T,
                            pDegree = 5,
                            max.c0.iter = 5,
                            c0.min = 1e-5) {


  n.points <- max(n.points, 10)

  list(n.points = n.points, nStart = nStart, nStart.lrt = nStart.lrt, 
       max.quantile = max.quantile,
       max.sparse = max.sparse,
       rel.info.range = rel.info.range, rel.info.guide = rel.info.guide,
       c0.start = c0.start, 
       aStart = aStart, bStart = bStart, gStart = gStart, start.sd = start.sd,
       penalty.Z0 = penalty.Z0, pDegree = pDegree, max.c0.iter = max.c0.iter,
       c0.min = c0.min,
       LRT.Z.select = LRT.Z.select, 
       LRT.Z0.select = LRT.Z0.select,
       LRT.Z.values = LRT.Z.values,
       zeroInflate = zeroInflate, 
       dense.add.0 = dense.add.0,
       only.integer = only.integer)
}

#' An S4 class object containing the DESCEND result for a single gene
#' 
#' The DESCEND class is a container to store the DESCEND result for a single gene. 
#'
#' @section Slots:
#' \describe{
#' \item{\code{distribution}}{The distribution of the deconvolved distribution with relative statistics}
#' \item{\code{estimates}}{The matrix of distribution measurements and coefficients estimated values, bias, standard deviation and mean square error (named as \code{DESCEND.sd})}
#' \item{\code{pval}}{The p values of the likelihood ratio tests if computed}
#' \item{\code{density.points}}{Smoothed version of the distribution for easier plotting of the distribution}
#' }
#'
#' @import methods
#' @name DESCEND
#' @rdname DESCEND
#' @aliases DESCEND-class
#' @author Jingshu Wang
#' @exportClass DESCEND

require(methods)

methods::setClass("DESCEND", representation(distribution = "matrix", 
                                   estimates = "matrix", 
                                   pval = "matrix", 
                                   density.points = "matrix"))
