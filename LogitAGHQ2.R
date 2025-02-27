# EM on Classification Model
# Code by Dr Hadrien Charvat
# Modified by Dr Maxime Bonjour & Damien Georges

# Based on : http://tinyheero.github.io/2016/01/03/gmm-em.html

# Version 1 : 20/04/2020
# Version 2 : 20/12/2023

Correct <- function(Res) {
  ## If evaluates at NA, return maximum value
  if (is.na(Res)) {
    Res <- .Machine$double.xmax
  }
  ## If evaluates at +/-infinity, return maximum value
  else if (abs(Res) == Inf) {
    Res <- sign(Res) * .Machine$double.xmax
  }
  return(Res)
}

Clust.Lik.C <- function(w, beta, data, num, denom, pos, cst) {
  beta[pos] <- beta[pos] + w
  Temp <- exp(data %*% beta)
  Denom <- 1 / (1 + Temp)
  lcond <- (num * log(Temp) + denom * log(Denom))
  res <- exp(sum(lcond) + length(num) * cst)
  return(res)
}

Clust.LogLik.C <- function(w, beta, data, num, denom, pos, cst) {
  beta[pos] <- beta[pos] + w
  lTemp <- (data %*% beta)
  lDenom <- -log(1 + exp(lTemp))
  lres <- sum(num * lTemp + denom * lDenom) + length(num) * cst
  return(lres)
}


Clust.LogLik.M <- function(beta, data, num, denom, pos, nbN, Sigma) {
  CstInf <- -1000
  CstSup <- 1000
  Cst <- 0

  ## Number of random effects
  n.rand <- dim(Sigma)[1]

  ## GH Quadrature nodes and weights
  GH <- gauss.quad(nbN, kind = "hermite")
  x.H <- GH$nodes
  eval(parse(text = paste("x.H <- t(as.matrix(expand.grid(", paste(rep("x.H", n.rand), collapse = ","), ")))")))
  x.H.2 <- rep(1, n.rand) %*% (x.H * x.H)
  eval(parse(text = paste("w.H <- t(as.matrix(expand.grid(", paste(rep("GH$weights", n.rand), collapse = ","), ")))")))
  log.rho.H <- log(w.H)
  log.rho.H <- apply(log.rho.H, 2, sum)
  log.rho.H[log.rho.H == -Inf] <- -.Machine$double.xmax
  nbNT <- length(log.rho.H)

  InvCov <- try(solve(Sigma), silent = TRUE)
  if (inherits(InvCov, "try-error")) {
    Res <- .Machine$double.xmax
  } else {
    ## (i) Define the integrand of the integral we want to calculate
    MyCst <- -log(prod(diag(chol(Sigma)))) - (n.rand / 2) * log(2 * pi)
    Integrand <- function(w) {
      exp(MyCst) * Clust.Lik.C(w, beta, data, num, denom, pos, Cst) * exp(-0.5 * w %*% InvCov %*% w)
    }

    ## (ii) Take the logarithm of the integrand
    ## (with a minus sign because we are going
    ## to give it to a minimisation algorithm...)
    Log.Int <- function(w) {
      Res <- -(MyCst + Clust.LogLik.C(w, beta, data, num, denom, pos, Cst)
        - 0.5 * w %*% InvCov %*% w)
      return(Correct(Res))
    }

    OptInt <- nlm(Log.Int, rep(0, n.rand))

    ## (iii) Find the mode of LogInt (i.e., the
    ## value for which the derivative becomes 0)
    MuHat <- OptInt$estimate

    ## (iv) Find the value of the second derivative
    ## of LogInt in order to define the standard error
    ## of the normal distribution used to approximate the integrand
    Hessian <- hessian(Log.Int, MuHat, method.args = list(r = 4))

    SigmaHat <- try(solve(Hessian), silent = TRUE)
    if (inherits(SigmaHat, "try-error")) {
      Res <- 1 / 0
    } else {
      SqrtSigma <- try(chol(SigmaHat), silent = TRUE)
      if (inherits(SqrtSigma, "try-error")) {
        Res <- 1 / 0
      } else {
        Nodes <- MuHat + sqrt(2) * (SqrtSigma %*% x.H)
        LCstMult <- log(prod(diag(SqrtSigma))) + n.rand / 2 * log(2)

        ## Approximate the integral by AGHQ
        Res <- 0
        for (k in 1:nbNT) {
          Res <- Res + exp(LCstMult + log.rho.H[k] + x.H.2[k] - Log.Int(Nodes[, k]))
        }
      }
    }

    Try <- 1
    while ((is.nan(Res) | (!is.nan(Res) & (Res == Inf | Res == 0))) & (Try < 25)) {
      if ((Log.Int(MuHat) == -.Machine$double.xmax) | (!is.nan(Res) & (Res == 1 / 0))) {
        CstSup <- Cst
      } else if ((Log.Int(MuHat) == .Machine$double.xmax) | (!is.nan(Res) & (Res == 0))) {
        CstInf <- Cst
      }
      Cst <- 0.5 * (CstInf + CstSup)

      OptInt <- nlm(Log.Int, rep(0, n.rand))
      MuHat <- OptInt$estimate
      Hessian <- hessian(Log.Int, MuHat, method.args = list(r = 4))
      SigmaHat <- try(solve(Hessian), silent = TRUE)
      if (inherits(SigmaHat, "try-error")) {
        Res <- 1 / 0
      } else {
        SqrtSigma <- try(chol(SigmaHat), silent = TRUE)
        if (inherits(SqrtSigma, "try-error")) {
          Res <- 1 / 0
        } else {
          Nodes <- MuHat + sqrt(2) * (SqrtSigma %*% x.H)
          LCstMult <- log(prod(diag(SqrtSigma))) + n.rand / 2 * log(2)
          Res <- 0
          for (k in 1:nbNT) {
            Res <- Res + exp(LCstMult + log.rho.H[k] + x.H.2[k] - Log.Int(Nodes[, k]))
          }
        }
      }
      Try <- Try + 1
    }
    ## Return (minus) the logarithm of the
    ## cluster-specific marginal likelihood
    Res <- -log(Res) + length(num) * Cst
  }

  return(as.vector(Correct(Res)))
}
