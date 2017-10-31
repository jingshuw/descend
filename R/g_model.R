#' Base function of G-modeling deconvolution
#'
#' Base function used by \code{\link{deconvSingle}} to deconvolve the underlying distribution. 
#' We assume X ~ F(T) where F is the noise distribution. We assume that 
#' \deqn{log(T) = offset + \gamma Z + \epsilon}
#' \deqn{P(T = 0) = \beta_0 + \beta_1 Z0}
#' The goal is the recover the distribution of exp(log(T) - offset - gamma Z), which has density g and is discretized at exp(tau) (add 0 when zero inflation is allowed). There can be some warning messages for the optimization process, which can be ignored.
#'
#' @param tau log of the discrete points of the deconvolved distribution
#' @param X a vector of observed counts
#' @param offset a vector with the same length as \code{X}. See details
#' @param family family of the noise distribution, support either "Poisson" or "Negative Binomial" with known tuning parameter
#' @param ignoreZero whether ignore the zero count. If true, then use truncated Poisson / Negative Binomial distribution. Default is False
#' @param zeroInflate whether add zero inflation part to the deconvolved distribution to reflect transcriptional bursting. Default is True. 
#' @param Z covariates for active intensity. Default is NULL. 
#' @param Z0 covariates for active fraction. Used only when zeroInflate is True. Default is NULL.
#' @param c0 the tuning parameter on the L2 penalty term. Default is 1. c0 will be selected automatically in \code{deconvSingle}
#' @param NB.size over-dispersion parameter when the family is Negative Binomial: mu = mu + mu^2/size
#' @param only.value whether not to compute the estimation statistics but only the value of the optimized lieklihood. Used for likelihood ratio test.
#' @param pDegree the degree of the Spline matrix. Default is 5.
#' @param aStart,bStart,gStart initial values of the density parameters, the coefficients of Z0 and coefficients of Z
#' @param ... extra parameters for the \code{\link[stats]{nlm}} function
#'
#' @return a list with elements
#' \item{stats}{a list of two elements. One is the \code{mat.dist}, which is the matrix of the estimated distribution. The other is \code{mat.coef}, which is the matrix of the coefficients of Z and Z0}
#' \item{mle}{the estimated parameters of the the density function}
#' \item{mle.g}{the estimated coefficients of Z}
#' \item{value}{the optimized penalized negative log-likelihood value}
#' \item{S}{the fake information proportion}
#' \item{cov}{the covariance of the parameters}
#' \item{bias}{the bias of the parameters}
#' \item{cov.g}{the covaraince of the estimated density points}
#' \item{cov.g.gamma}{the covariance between the estimated density points and the coefficient of Z}
#' \item{loglik}{the objective function being optimized}
#' \item{statsFunction}{the function computing the relavant statistics}
#'
#' @note This is an extension of the G-modeling package
#' @examples
#'
#' X <- rpois(1000, 0.2 * 3)
#' lam.max <- quantile(X[!X == 0], probs = c(0.98))/0.2
#' tau <- seq(0, lam.max, length.out = 51)[-1]
#' 
#' Z <- rnorm(1000)
#' 
#' result <- deconvG(log(tau), X, zeroInflate = TRUE,
#'                      Z = Z,
#'                      Z0 = Z,
#'                      offset = rep(log(0.2), 1000))
#' 
#' @import stats
#' @rdname deconvG
#' @export

deconvG <- 
function(tau, X,  
         offset, 
         family = c("Poisson", "Negative Binomial"),
         ignoreZero = F, ## whether ignore zero count
         zeroInflate = F, ## whether add zero inflation part to the deconvolved distribution
         Z = NULL,
         Z0 = NULL,
         c0 = 1, 
         NB.size = 100, ## over-dispersion parameter when family is Negative Binomial
         only.value = F,
         pDegree = 5, 
         aStart = 1,
         bStart = 0,
         gStart = 0,
         ...) 
{

#  print(head(X))

  family <- match.arg(family)
  requireNamespace("stats")
  if (!is.null(Z))
    Z <- as.matrix(Z)
  if (!is.null(Z0)) {
    Z0 <- as.matrix(Z0)
    Z0.means <- colMeans(Z0)
    Z0 <- scale(Z0, scale = F)
  }
    
  if (ignoreZero) {
    if (!is.null(Z))
      Z <- Z[X!=0, , drop= F]
    if (!missing(offset))
      offset <- offset[X!=0]
    X <- X[X != 0]
  }

  N <- length(X)
  if (missing(offset))
    offset <- rep(0, N)

  y <- rep(1, N)
  lam <- exp(tau)

  size <- NB.size


##################### Construction step ################################################  

  if (is.null(Z)) {  ## then P is a matrix of N * m
    Q <- scale(splines::ns(lam, pDegree), center = TRUE, 
               scale = FALSE)
    Q <- svd(Q)$u
    Q <- apply(Q, 2, function(w) w/sqrt(sum(w * w))) ## matrix of m * (pDegree + 1)
    if (ignoreZero) {
      if (family == "Poisson")
        P <- sapply(tau, function(theta) 
                    dpois(x = X, 
                          lambda = exp(theta + offset))/(1 - exp(-(exp(theta + offset)))))
      else if (family == "Negative Binomial")
        P <- sapply(tau, function(theta) 
                    dnbinom(x = X, 
                            mu = exp(theta + offset), size = size)/
                    (1 - dnbinom(0, mu = exp(theta + offset), size = size)))
    } else {
      # P is dimension N * m
      if (family == "Poisson")
        P <- sapply(tau, function(theta) dpois(x = X, 
                                               lambda = exp(theta + offset)))
      else if (family == "Negative Binomial")
        P <- sapply(tau, function(theta) 
                    dnbinom(x = X, 
                            mu = exp(theta + offset), size = size))
      if (zeroInflate) {
        #      print("check")
        P <- cbind(as.numeric(X == 0), P)
        lam <- c(0, lam)
        Q <- scale(splines::ns(lam, pDegree), center = TRUE, 
                   scale = FALSE)
        Q <- svd(Q)$u
        Q <- apply(Q, 2, function(w) w/sqrt(sum(w * w))) ## matrix of (m + 1) * (pDegree)

        if (is.null(Z0))
          Q <- cbind(c(1, rep(0, length(tau))), Q)
        else
          Q.fun <- function(Z.vec) {
            add.on <- rbind(c(1, Z.vec), matrix(0, length(tau), length(Z.vec)+1))
            return(cbind(add.on, Q))
          }
      }
    }
  } else { ## when Z is not NULL
    Q <- scale(splines::ns(lam, pDegree), center = TRUE, 
               scale = FALSE)
    Q <- svd(Q)$u
    Q <- apply(Q, 2, function(w) w/sqrt(sum(w * w))) ## matrix of m * (pDegree + 1)
    if (ignoreZero) {
      P.fun <- function(gamma) 
        sapply(tau, function(theta) { 
               if (family == "Poisson")
                 dpois(x = X, 
                       lambda = exp(theta + offset + Z %*% gamma))/ 
               (1 - exp(-(exp(theta + offset + Z %*% gamma))))
               else if (family == "Negative Binomial")
                 dnbinom(x = X, 
                         mu = exp(theta + offset + Z%*% gamma), size = size)/
               (1 - dnbinom(0, mu = exp(theta + offset + Z %*% gamma), size = size))
               })
    } else {
      P.fun <- function(gamma) {
        P <- sapply(tau, function(theta) {
                    if (family == "Poisson")
                      dpois(x = X, 
                            lambda = exp(theta + offset + Z %*% gamma))
                    else if (family == "Negative Binomial")
                      dnbinom(x = X, 
                              mu = exp(theta + offset + Z %*% gamma), size = size)
               })
        if (zeroInflate) 
          P <- cbind(as.numeric(X == 0), P)
        return(P)
      }
      if (zeroInflate) {
        lam <- c(0, lam)
        Q <- scale(splines::ns(lam, pDegree), center = TRUE, 
                   scale = FALSE)
        Q <- svd(Q)$u
        Q <- apply(Q, 2, function(w) w/sqrt(sum(w * w))) ## matrix of (m+1) * (pDegree)

        if (is.null(Z0))
          Q <- cbind(c(1, rep(0, length(tau))), Q)
        else
          Q.fun <- function(Z.vec) {
            add.on <- rbind(c(1, Z.vec), matrix(0, length(tau), length(Z.vec)+1))
            return(cbind(add.on, Q))
          }
      }
    }
  }

######################### Optimization Step ########################################
  # Initialization
  p <- ncol(Q)
  if (!is.null(Z0))
    p <- p + 1 + ncol(Z0)

    
  pGiven <- length(aStart)
  if (pGiven == 1) {
    aStart <- rep(aStart, p)
    if (!is.null(Z0)) {
      if (length(bStart) == 1)
        bStart <- rep(bStart, ncol(Z0))
      aStart[2:(1+ ncol(Z0))] <- bStart
    }
  }
  else {
    if (pGiven != p) 
      stop(sprintf("Wrong length (%d) for initial parameter, expecting length (%d)", 
                   pGiven, p))
  }

  if (!is.null(Z)) {
    if (length(gStart) == 1)
      gStart <- rep(gStart, ncol(Z))
  }

  if (!is.null(Z0)) 
    g.fun <- function(Z.vec, a) {
      Qi <- Q.fun(Z.vec)
      gi <- exp(Qi %*% a)
      gi <- as.vector(gi /sum(gi))
      return(gi)
    }


  ## optimization objective function
  if (is.null(Z)) {
    loglik <- function(a) {
      if (is.null(Z0)) {
        g <- exp(Q %*% a)
        g <- as.vector(g/sum(g))
        f <- as.vector(P %*% g)
        ## objective value
        value <- -sum(y * log(f)) + c0 * sum(a^2)^0.5
        Pt <- P/f
        W <- g * (t(Pt) - 1)
        ## derivative
        qw <- t(Q) %*% W
      } else {
        m.temp <- 1 + ncol(Z0)
        g.temp <- Q %*% a[-(1:m.temp)]
        G.temp1 <- cbind(rep(1, nrow(Z0)), Z0) %*% a[1:m.temp]
        G <- matrix(rep(g.temp, nrow(Z0)), nrow = length(g.temp))
        G[1, ] <- G[1, ] + G.temp1
        G <- exp(G)
        G <- t(G) / colSums(G)
        f <- rowSums(P * G)
        value <- -sum(y * log(f)) + c0 * sum(a^2)^0.5
        Pt <- P/f
        W <- t(G * (Pt - 1)) # dimension is as G: N * (m + 1)

        qw.temp2 <- t(Q) %*% W
        qw.temp1 <- cbind(rep(1, nrow(Z0)), Z0) * W[1, ]
        qw <- rbind(t(qw.temp1), qw.temp2)

      }
      aa <- sqrt(sum(a^2))
      sDot <- c0 * a/aa

      attr(value, "gradient") <- -(qw %*% y) + sDot
      value
    }

    result <- stats::nlm(f = loglik, p = aStart, gradtol = 1e-5,
                         ...)
    mle <- result$estimate
    value <- loglik(mle)
    mle.a <- result$estimate
    mle.g <- NULL

  } else {
    loglik <- function(a.gamma) {
      a <- a.gamma[1:p]
      gamma <- a.gamma[(p+1):length(a.gamma)]
      P <- P.fun(gamma)
      if (is.null(Z0)) {
        g <- exp(Q %*% a)
        g <- as.vector(g/sum(g))
        f <- as.vector(P %*% g)
        ## objective value
        value <- -sum(y * log(f)) + c0 * sum(a^2)^0.5
        Pt <- P/f
        W <- g * (t(Pt) - 1)
        ## derivative
        qw <- t(Q) %*% W
      } else {
        m.temp <- 1 + ncol(Z0)
        g.temp <- Q %*% a[-(1:m.temp)]
        G.temp1 <- cbind(rep(1, nrow(Z0)), Z0) %*% a[1:m.temp]
        G <- matrix(rep(g.temp, nrow(Z0)), nrow = length(g.temp))
        G[1, ] <- G[1, ] + G.temp1
        G <- exp(G)
        G <- t(G) / colSums(G)

        f <- rowSums(P * G)
        value <- -sum(y * log(f)) + c0 * sum(a^2)^0.5
        Pt <- P/f
        W <- t(G * (Pt - 1)) # dimension is as G: N * (m + 1)

        qw.temp2 <- t(Q) %*% W
        qw.temp1 <- cbind(rep(1, nrow(Z0)), Z0) * W[1, ]
        qw <- rbind(t(qw.temp1), qw.temp2)
      }

      aa <- sqrt(sum(a^2))
      sDot <- c0 * a/aa
      dot.l.a <- -rowSums(qw) + sDot # length p

      # derivative for gamma
      eta <- rep(1, N) %*% t(tau) + as.vector(Z %*% gamma) + offset
      if (family == "Poisson" || family == "Negative Binomial")
        mu <- exp(eta)
      if (zeroInflate){
        P.tilde <- P[, -1]
        if (is.null(Z0))
          g.tilde <- g[-1]
        else
          G.tilde <- G[, -1]
      } else {
        P.tilde <- P
        g.tilde <- g
      }
    #  print(head(X, 3))
      if (family == "Poisson")
        A <- (X - mu) * P.tilde # dimension is N * m
      if (family == "Negative Binomial")
        A <- (X - mu) / (1 + mu/size) * P.tilde
      if (is.null(Z0))
        agf <- A %*% g.tilde / f
      else
        agf <- rowSums(A * G.tilde)/f
      dot.l.gamma <- -t(Z) %*% agf

      attr(value, "gradient") <- c(dot.l.a, dot.l.gamma)

      value
    }
    result <- stats::nlm(f = loglik, p = c(aStart, gStart), 
                         gradtol = 1e-5, 
                        ...)

    mle <- result$estimate
    mle.a <- mle[1:p]
    mle.g <- mle[-(1:p)]
    value <- loglik(mle)
  }


################################# inference functions #####################################

  ## used when Z is NULL
  statsFunction <- function(a) {
    if (is.null(Z0)) {
      g <- exp(Q %*% a)
      g <- as.vector(g/sum(g))
      f <- as.vector(P %*% g)
      ## objective value

      value <- -sum(y * log(f)) + c0 * sum(a^2)^0.5
      Pt <- P/f
      W <- g * (t(Pt) - 1)
      ## derivative
      qw <- t(Q) %*% W
      qG <- matrix(rep(t(Q) %*% g, ncol(W)), ncol = ncol(W))
      Wplus <- rowSums(W)
      qWq <- t(Q) %*% (Wplus * Q)
    } else {
      m.temp <- 1 + ncol(Z0)
      g.temp <- Q %*% a[-(1:m.temp)]
      G.temp1 <- cbind(rep(1, nrow(Z0)), Z0) %*% a[1:m.temp]
      G <- matrix(rep(g.temp, nrow(Z0)), nrow = length(g.temp))
      G[1, ] <- G[1, ] + G.temp1
      G <- exp(G)
      G <- t(G) / colSums(G)

      f <- rowSums(P * G)

      value <- sum(y * log(f)) + c0 * sum(a^2)^0.5
      Pt <- P/f
      W <- t(G * (Pt - 1)) # dimension is as G: N * (m + 1)
      qw <- sapply(1:nrow(Z0), function(i) {
                   Qi <- Q.fun(Z0[i, ])
                   return(t(Qi) %*% W[, i])
                         })
      g <- colMeans(G)
      qG <- sapply(1:nrow(Z0), function(i) {
                         Qi <- Q.fun(Z0[i, ])
                         return(t(Qi) %*% G[i, ])
                         })
      qWq <- matrix(0, nrow(qw), nrow(qw))
      for (i in 1:nrow(Z0)) {
        Qi <- Q.fun(Z0[i, ])
        qWq <- qWq + t(Qi) %*% (W[, i] * Qi)
      }
    }


    yHat <- if (length(y) == 1 && y == 1 || length(y) == length(X)) 
      y
    else sum(y) * f
    # estimate W
    I1 <- qw %*% (yHat * t(qw))

    aa <- sqrt(sum(a^2))
    sDot <- c0 * a/aa
    sDotDot <- (c0/aa) * (diag(length(a)) - outer(a, a)/aa^2)
    R <- sum(diag(sDotDot))/sum(diag(I1))
    temp <- qw %*% t(qG)
    I2 <- solve(I1 + temp + t(temp) - qWq + sDotDot)

    bias <- as.vector(-I2 %*% sDot)
    Cov <- I2 %*% (I1 %*% t(I2))

    if (is.null(Z0))
      mat.coef <- NULL
    else {
      idx <- 2:(1+ncol(Z0))
      mat.coef <- cbind(mle[idx], 
                         bias[idx],
                         sqrt(diag(Cov[idx, idx, drop = F])))
      Q0 <- cbind(rbind(c(1, -Z0.means), matrix(0, nrow(Q) - 1, 1 + ncol(Z0))), Q)
      temp <- as.vector(exp(Q0 %*% a))
      beta0 <- log(temp[1]/sum(temp[-1]))
      temp.dev <- temp * Q0
      beta0.dev <- c(1/temp[1], -rep(1/sum(temp[-1]), length(temp)-1))
      beta0.dev <- t(temp.dev) %*% beta0.dev
      beta0 <- c(beta0, t(beta0.dev) %*% bias, 
                 sqrt(t(beta0.dev) %*% Cov %*% beta0.dev))
      mat.coef <- rbind(beta0, mat.coef)
      rownames(mat.coef) <- paste("Z0 effect: beta", 0:ncol(Z0), sep = "")
      colnames(mat.coef) <- c("est", "bias", "sd")
    }


    if (is.null(Z0)) {
      Dq <- (diag(g) - outer(g, g)) %*% Q
    } else {
      Dq <- matrix(0, ncol(G), p)
      sapply(1:nrow(G), function(i) {
             Qi <- Q.fun(Z0[i, ])
             temp <- (diag(G[i, ]) - outer(G[i, ], G[i, ])) %*% Qi
             Dq <<- Dq + temp
                         })
      Dq <- Dq/nrow(G)
    }
    bias.g <- Dq %*% bias
    Cov.g <- Dq %*% Cov %*% t(Dq)
    se.g <- diag(Cov.g)^0.5
    D <- diag(length(lam))
    D[lower.tri(D)] <- 1
    accuG <- cumsum(g)
    Cov.G <- D %*% (Cov.g %*% t(D))
    se.G <- diag(Cov.G)^0.5
    mat <- cbind(lam, g, se.g, accuG, se.G, bias.g)
    colnames(mat) = c("theta", "g", "SE.g", "G", "SE.G", 
                      "Bias.g")
    list(S = R, cov = Cov, bias = bias, cov.g = Cov.g, cov.g.gamma= NULL,
         mat = list(mat.dist = mat, mat.coef = mat.coef))
  }

  

  ## Used when Z is not NULL
  statsFunctionZ <- function(a, gamma) {
    P <- P.fun(gamma)
    if (is.null(Z0)) {
      g <- exp(Q %*% a)
      g <- as.vector(g/sum(g))
      f <- as.vector(P %*% g)
      ## objective value
      value <- -sum(y * log(f)) + c0 * sum(a^2)^0.5
      Pt <- P/f
      W <- g * (t(Pt) - 1)
      ## derivative
      qw <- t(Q) %*% W
      qG <- matrix(rep(t(Q) %*% g, ncol(W)), ncol = ncol(W))
      Wplus <- rowSums(W)
      qWq <- t(Q) %*% (Wplus * Q)
    } else {
      m.temp <- 1 + ncol(Z0)
      g.temp <- Q %*% a[-(1:m.temp)]
      G.temp1 <- cbind(rep(1, nrow(Z0)), Z0) %*% a[1:m.temp]
      G <- matrix(rep(g.temp, nrow(Z0)), nrow = length(g.temp))
      G[1, ] <- G[1, ] + G.temp1
      G <- exp(G)
      G <- t(G) / colSums(G)

      f <- rowSums(P * G)
      value <- sum(y * log(f)) + c0 * sum(a^2)^0.5
      Pt <- P/f
      W <- t(G * (Pt - 1)) # dimension is as G: N * (m + 1)
      qw <- sapply(1:nrow(Z0), function(i) {
                   Qi <- Q.fun(Z0[i, ])
                   return(t(Qi) %*% W[, i])
                      })
      g <- colMeans(G)
      qG <- sapply(1:nrow(Z0), function(i) {
                         Qi <- Q.fun(Z0[i, ])
                         return(t(Qi) %*% G[i, ])
                         })
      qWq <- matrix(0, nrow(qw), nrow(qw))
      for (i in 1:nrow(Z0)) {
        Qi <- Q.fun(Z0[i, ])
        qWq <- qWq + t(Qi) %*% (W[, i] * Qi)
      }
    }

    ## compute Igg
    eta <- rep(1, N) %*% t(tau) + as.vector(Z %*% gamma) + offset
    if (family == "Poisson" || family == "Negative Binomial") {
      mu <- exp(eta)
    } 
    if (zeroInflate){
      P.tilde <- P[, -1]
      if (is.null(Z0)) {
        g.tilde <- g[-1]
        Q.tilde <- Q[-1, ]
      }
      else {
        G.tilde <- G[, -1]
        g.tilde <- G.tilde[1, ]
        Q.tilde <- cbind(matrix(0, nrow(Q)-1, 1 + ncol(Z0)), Q[-1, ])
      }
    } else {
      P.tilde <- P
      g.tilde <- g
      Q.tilde <- Q
    }
    if (family == "Poisson") {
      A <- (X - mu) * P.tilde # dimension is N * m
      B <- ((X - mu)^2 - mu) * P.tilde
    }
    if (family == "Negative Binomial") {
      A <- (X - mu) / (1 + mu/size) * P.tilde
      B <- ((X - mu)^2 / (1 + mu/size)^2 - mu/(1 + mu/size) - 
            (X - mu) * mu / (1 + mu/size)^2/size) * P.tilde
    }


    if (is.null(Z0)) {
      agf <- as.vector(A %*% g.tilde / f)
      bgf <- as.vector(B %*% g.tilde / f)
    }
    else {
      agf <- rowSums(A * G.tilde)/f
      bgf <- rowSums(B * G.tilde)/f
    }


    dot.li.gamma <- t(agf * Z)

    deriv.comb <- rbind(qw, dot.li.gamma)
    aa <- sqrt(sum(a^2))
    sDot <- c0 * a/aa
    sDotDot <- (c0/aa) * (diag(length(a)) - outer(a, a)/aa^2)

    I1full <- deriv.comb %*% t(deriv.comb)
    Iaa <- I1full[1:length(a), 1:length(a)]
    Igg <- t(Z) %*% (as.vector(agf^2 - bgf) * Z)
    if (is.null(Z0))
      Iag <- t(Q.tilde) %*% (g.tilde * (t(P.tilde/f) %*% (agf * Z) - t(A/f) %*% Z))
    else
      Iag <- t(Q.tilde) %*% (t(G.tilde * P.tilde/f) %*% (agf * Z) - t(G.tilde * A/f) %*% Z)
    I2full <- rbind(cbind(Iaa, Iag), cbind(t(Iag), Igg))
    I2full[1:length(a), 1:length(a)] <- I2full[1:length(a), 1:length(a)] + 
                                        sDotDot

    Ifull.inv <- solve(I2full)

      ## ratio of artificial to genuine information
    R <- sum(diag(sDotDot))/sum(diag(I1full[1:length(a), 1:length(a)])) 


    bias <- as.vector(-Ifull.inv %*% c(sDot, rep(0, ncol(Z))))
    Cov <- Ifull.inv %*% (I1full %*% t(Ifull.inv))
    if (is.null(Z0)) {
      mat.coef <- cbind(mle.g, bias[-(1:length(a))],
                        sqrt(diag(Cov[-(1:length(a)), -(1:length(a)), drop = F])))
      rownames(mat.coef) <- c(paste("Z:gamma", 1:ncol(Z), sep = ""))
      colnames(mat.coef) <- c("est", "bias", "sd")
    }
    else {
      idx1 <- 2:(1+ncol(Z0))
      idx <- c((length(a)+1):length(bias), idx1)
      mat.coef <- cbind(mle[idx], 
                         bias[idx],
                         sqrt(diag(Cov[idx, idx, drop = F])))

      Q0 <- cbind(rbind(c(1, -Z0.means), matrix(0, nrow(Q) - 1, 1 + ncol(Z0))), Q)
      temp <- as.vector(exp(Q0 %*% a))
      beta0 <- log(temp[1]/sum(temp[-1]))
      temp.dev <- temp * Q0
      beta0.dev <- c(1/temp[1], -rep(1/sum(temp[-1]), length(temp)-1))
      beta0.dev <- t(temp.dev) %*% beta0.dev
      beta0 <- c(beta0, t(beta0.dev) %*% bias[1:length(a)], 
                 sqrt(t(beta0.dev) %*% Cov[1:length(a), 1:length(a), drop = F] %*% beta0.dev))
      mat.coef <- rbind(mat.coef[1:ncol(Z), ],
                        beta0,
                        mat.coef[-(1:ncol(Z)), ])

      rownames(mat.coef) <- c(paste("Z effect: gamma", 1:ncol(Z), sep = ""),
                              paste("Z0 effect: beta", 0:ncol(Z0), sep = ""))
      colnames(mat.coef) <- c("est", "bias", "sd")
    }

    if (is.null(Z0)) {
      Dq <- (diag(g) - outer(g, g)) %*% Q
    } else {
      Dq <- matrix(0, ncol(G), p)
      sapply(1:nrow(G), function(i) {
             Qi <- Q.fun(Z0[i, ])
             temp <- (diag(G[i, ]) - outer(G[i, ], G[i, ])) %*% Qi
             Dq <<- Dq + temp
                      })
      Dq <- Dq/nrow(G)
    }

    bias.g <- Dq %*% bias[1:p]
    Cov.g <- Dq %*% Cov[1:p, 1:p] %*% t(Dq)
    Cov.g.gamma <- Dq %*% Cov[1:p, -(1:p), drop = F]
    se.g <- diag(Cov.g)^0.5
    D <- diag(length(lam))
    D[lower.tri(D)] <- 1
    accuG <- cumsum(g)
    Cov.G <- D %*% (Cov.g %*% t(D))
    se.G <- diag(Cov.G)^0.5
    mat <- cbind(lam, g, se.g, accuG, se.G, bias.g)
    colnames(mat) = c("theta", "g", "SE.g", "G", "SE.G", 
                      "Bias.g")
    list(S = R, cov = Cov, bias = bias, cov.g = Cov.g, cov.g.gamma = Cov.g.gamma,
         mat = list(mat.dist = mat, mat.coef = mat.coef))
  }


  if (only.value)
    return(value)

  if (is.null(Z))
    stats <- statsFunction(mle.a)
  else
    stats <- statsFunctionZ(mle.a, mle.g)
  list(stats = stats$mat, 
       mle = mle.a, 
       mle.g = mle.g,
       value = value,
       S = stats$S, 
       cov = stats$cov, 
       bias = stats$bias,
       cov.g = stats$cov.g, 
       cov.g.gamma = stats$cov.g.gamma,
       loglik = loglik, 
       statsFunction = statsFunction)
}



