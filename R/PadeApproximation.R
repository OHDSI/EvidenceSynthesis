# Copyright 2021 Observational Health Data Sciences and Informatics
#
# This file is part of EvidenceSynthesis
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


pade <- function(x, beta, a0, a1, a2, b1, b2) {
  delta <- x - beta
  numerator <- cbind(1, delta, delta^2) %*% c(a0, a1, a2)
  denominator <- cbind(1, delta, delta^2) %*% c(1, b1, b2)
  return(numerator / denominator)
}

padeGradient<- function(x, beta, a0, a1, a2, b1, b2) {
  delta <- x - beta
  Pm <- cbind(1, delta, delta^2) %*% c(a0, a1, a2)
  Qn <- cbind(1, delta, delta^2) %*% c(1, b1, b2)
  GradPade_num <- Qn*(a1 + 2*a2*delta) - Pm*(b1+2*b2*delta)
  GradPade <- GradPade_num / (Qn^2)
  return(as.vector(GradPade))
}

fitPadeUsingSumSquares <- function(beta, ll, weighByLikelihood = TRUE) {

  sumSquares <- function(p, maxAnchor = TRUE, idx = NULL) {
    constraints <- getPadeConstraints(p[1], p[2], p[3], p[4], p[5], p[6])
    approxLl <- padeConstrained(beta, p[1], p[2], p[3], p[4], p[5], p[6], constraints$minBeta, constraints$maxBeta, constraints$minD1, constraints$maxD1)
    if (maxAnchor) {
      approxLl <- approxLl - max(approxLl)
    } else {
      approxLl <- approxLl - approxLl[idx]
    }
    if (weighByLikelihood) {
      result <- sum(weights * (approxLl - ll)^2)
    } else {
      result <- sum((approxLl - ll)^2)
    }
    return(result)
  }

  beta <- beta[!is.nan(ll)]
  ll <- ll[!is.nan(ll)]

  # Scale to standard (translate in log space so max at 0):
  ll <- ll - max(ll)
  weights <- exp(ll)
  # Weights shouldn't be too small:
  weights[weights < 0.001] <- 0.001

  mode <- beta[ll == 0][1]
  if (mode == min(beta) || mode == max(beta)) {
    mode <- 0
  }

  if (min(ll) > -1e-06) {
    return(data.frame(beta = 0,
                      a0 = 0,
                      a1 = 0,
                      a2 = 0,
                      b1 = 0,
                      b2 = 0))
  }

  fit <- tryCatch({
    suppressWarnings(nlm(sumSquares, c(mode, 0, 0 ,0 ,0, 0), maxAnchor = TRUE))
  }, error = function(e) {
    list(estimate = rep(0, 6), minimum = Inf)
  })
  result <- data.frame(beta = fit$estimate[1], a0 = fit$estimate[2], a1 = fit$estimate[3], a2 = fit$estimate[4], b1 = fit$estimate[5], b2 = fit$estimate[6])
  minimum <- fit$minimum

  fit <- tryCatch({
    suppressWarnings(optim(c(mode, 0, 0 ,0 ,0, 0), sumSquares, maxAnchor = TRUE))
  }, error = function(e) {
    list(par = rep(0, 6), value = Inf)
  })
  if (fit$value < minimum) {
    result <- data.frame(beta = fit$par[1], a0 = fit$par[2], a1 = fit$par[3], a2 = fit$par[4], b1 = fit$par[5], b2 = fit$par[6])
    minimum <- fit$value
  }

  # Scale to standard (translate in log space so intersects at (0,0)):
  idx <- which(abs(beta) == min(abs(beta)))
  ll <- ll - ll[idx]
  fit <- tryCatch({
    suppressWarnings(nlm(sumSquares, c(mode, 0, 0 ,0 ,0, 0), maxAnchor = FALSE, idx = idx))
  }, error = function(e) {
    list(estimate = rep(0, 6), minimum = Inf)
  })
  if (fit$minimum < minimum) {
    result <- data.frame(beta = fit$estimate[1], a0 = fit$estimate[2], a1 = fit$estimate[3], a2 = fit$estimate[4], b1 = fit$estimate[5], b2 = fit$estimate[6])
    minimum <- fit$minimum
  }
  fit <- tryCatch({
    suppressWarnings(optim(c(mode, 0, 0 ,0 ,0, 0), sumSquares, maxAnchor = FALSE, idx = idx))
  }, error = function(e) {
    list(par = rep(0, 6), value = Inf)
  })
  if (fit$value < minimum) {
    result <- data.frame(beta = fit$par[1], a0 = fit$par[2], a1 = fit$par[3], a2 = fit$par[4], b1 = fit$par[5], b2 = fit$par[6])
    minimum <- fit$value
  }
  threshold <- 0.05
  if (minimum/length(beta) > threshold) {
    warn(paste("Mean squared error greater than ", threshold, ". Probably bad fit."))
  }
  return(result)
}


getPadeConstraints <- function(beta, a0, a1, a2, b1, b2) {
  # Solve where second derivative of Pade is 0:
  aa <- -2*b2*(a2*b1-a1*b2)
  bb <- 6*b2*(a0*b2-a2)
  cc <- 2*(a2*b1-a1*b2) - 4*b2*(a1-a0*b1)- 2*(a2-a0*b2)*b1
  dd <- 2*(a2-a0*b2) - 2*b1*(a1-a0*b1)
  pp <- c(aa,bb,cc,dd)
  rootD2 <- tryCatch(cubic(pp), error = function(e) c(-Inf, Inf))
  if (is.complex(rootD2)) {
    rootD2 <- c(-Inf, Inf)
  }

  # Set min and max beta between points where second derivate is 0:
  minBeta <- max(rootD2[rootD2 < beta])
  maxBeta <- min(rootD2[rootD2 > beta])

  # Compute first derivatives at min and max beta:
  minD1 <- max(0, padeGradient(minBeta, beta, a0, a1, a2, b1, b2))
  maxD1 <- min(0, padeGradient(maxBeta, beta, a0, a1, a2, b1, b2))

  return(data.frame(minBeta = minBeta,
                    maxBeta = maxBeta,
                    minD1 = minD1,
                    maxD1 = maxD1))
}

# Copied from RConics:
cubic <- function(p) {
  if (abs(p[1]) < .Machine$double.eps^0.95) {
    stop("Coefficient of highest power must not be zero!\n")
  }
  if (!all(is.numeric(p)) || length(p) != 4) {
    stop("p is not a numeric or/and has not 4 elements!\n")
  }
  if (any(is.complex(p))) {
    stop("the coefficient must be real")
  }
  a <- numeric(3)
  for (i in 2:4) {
    a[i - 1] = p[i]/p[1]
  }
  Q <- (a[1]^2 - 3 * a[2])/9
  R <- (2 * a[1]^3 - 9 * a[1] * a[2] + 27 * a[3])/54
  x <- numeric(3)
  if (R^2 < Q^3) {
    theta <- acos(R/sqrt(Q^3))
    x[1] <- -2 * sqrt(Q) * cos(theta/3) - a[1]/3
    x[2] <- -2 * sqrt(Q) * cos((theta + 2 * pi)/3) - a[1]/3
    x[3] <- -2 * sqrt(Q) * cos((theta - 2 * pi)/3) - a[1]/3
  }
  else {
    A = -sign(R) * (abs(R) + sqrt(R^2 - Q^3))^(1/3)
    if (isTRUE(all.equal(0, A))) {
      B <- 0
    }
    else {
      B <- Q/A
    }
    x[1] <- (A + B) - a[1]/3
    x[2] <- -0.5 * (A + B) - a[1]/3 + sqrt(3) * complex(real = 0,
                                                        imaginary = 1) * (A - B)/2
    x[3] <- -0.5 * (A + B) - a[1]/3 - sqrt(3) * complex(real = 0,
                                                        imaginary = 1) * (A - B)/2
  }
  return(x)
}

padeConstrained <- function(x, beta, a0, a1, a2, b1, b2, minBeta = -Inf, maxBeta = Inf, minD1 = 0, maxD1 = 0) {
  ll <- rep(0, length(x))

  idx <- which(x < minBeta)
  if (length(idx) > 0) {
    ll[idx] <- c(pade(minBeta, beta, a0, a1, a2, b1, b2)) + (x[idx] - minBeta) * minD1
  }
  idx <- which(x > maxBeta)
  if (length(idx) > 0) {
    ll[idx] <- c(pade(maxBeta, beta, a0, a1, a2, b1, b2)) + (x[idx] - maxBeta) * maxD1
  }
  idx <- which(x >= minBeta & x <= maxBeta)
  if (length(idx) > 0) {
    ll[idx] <- pade(x[idx], beta, a0, a1, a2, b1, b2)
  }
  return(ll)
}

#' Approximate Cox regression using Pade approximation
#'
#' @param survivalData    A data frame with the person-level survival data. See details.
#' @param bounds          The bounds on the effect size used to fit the approximation.
#'
#' @details
#' The `survivalData` data frame should have the following columns:
#' - **x**: the binary exposure variable.
#' - **time**: the time to either the outcome or censoring.
#' - **y**: indicator whether at the end of `time` the outcome occurred (1) or censoring (0).
#' - **stratumId**: an integer indicating the stratum (e.g. propensity-score strata).
#'
#' @return
#' A data frame with the parameters of the Pade approximation.
#'
#' @export
approximateLikelihoodUsingPade <- function(survivalData,
                                           bounds = c(log(0.1), log(10))){
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(survivalData, add = errorMessages)
  checkmate::assertNames(colnames(survivalData), must.include = c("x", "y", "time", "stratumId"), add = errorMessages)
  checkmate::assertNumeric(bounds, len = 2, add = errorMessages)
  checkmate::assertTRUE(bounds[1] < bounds[2], add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  ensure_installed("pda")

  cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                            data = survivalData,
                                            modelType = "cox")
  cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
  if (cyclopsFit$return_flag=="ILLCONDITIONED" || any(is.na(confint(cyclopsFit, 1)))) {
    message("Computing Pade approximation using sum of squares method")
    x <- seq(bounds[1], bounds[2], length.out = 100)
    ll <- Cyclops::getCyclopsProfileLogLikelihood(cyclopsFit, 1, x)$value
    return(fitPadeUsingSumSquares(x, ll))
  } else {
    message("Computing Pade approximation using derivatives")
    beta <- coef(cyclopsFit)
    confint(cyclopsFit, 1)
    llDeriv <- computeLlDerivatives(survivalData, beta)

    c0 <- llDeriv[1]
    c1 <- llDeriv[2]
    c2 <- llDeriv[3]/2
    c3 <- llDeriv[4]/6
    c4 <- llDeriv[5]/24

    b2 <- (c3^2-c2*c4)/(c2^2-c1*c3)
    b1 <- -c3/c2-c1/c2*b2
    a0 <- c0
    a1 <- c1+c0*b1
    a2 <- c2+c1*b1+c0*b2
    return(data.frame(beta = beta,
                      a0 = a0,
                      a1 = a1,
                      a2 = a2,
                      b1 = b1,
                      b2 = b2))
  }
}

# rcpp_aggregate <- function(x, indices, simplify = TRUE, cumulative = FALSE, reversely = FALSE) {
#   .Call('_pda_rcpp_aggregate', PACKAGE = 'pda', x, indices, simplify, cumulative, reversely)
# }

computeLlSingleDerivatives <- function(bbar,ipdata) {
  if(sum(ipdata$status==TRUE)==0){
    logL=0;logL_D1=0;logL_D2=0;logL_D3=0;logL_D4=0
  }else{

    T_all<-sort(unique(ipdata$time[ipdata$status==TRUE]))
    nt<-length(T_all)
    t_max <- max(ipdata$time)+1
    #ipdata0 <- ipdata
    ipdata <- ipdata[!ipdata$time<min(T_all),]

    #generate dataframe in format expected by ODAC
    pfdata <- cbind(T_all, 0, rep(0))
    colnames(pfdata)<-colnames(ipdata)
    pfdata <- rbind(ipdata, pfdata)
    pfdata <- pfdata[order(pfdata$time),]
    pfdata$interval <- cut(pfdata$time, breaks = c(T_all, t_max), labels = 1:nt, right=FALSE)
    pfdata <- pfdata[order(pfdata$interval),]
    pfdata$interval[is.na(pfdata$interval)]<-nt
    X <- c(pfdata$x)
    # summary stats: U, W, Z
    eXb <- c(exp(X*bbar))
    UWZ <- eXb * cbind(1, X, X^2, X^3, X^4)
    # rcpp_aggregate() is a function written in rcpp for calculating column-wise (reverse) cumsum
    # credit to Dr Wenjie Wang
    UWZ <- pda:::rcpp_aggregate(x = UWZ, indices = pfdata$interval, cumulative = T, reversely = T)

    # since fake X=0, cumulative W and Z will be the same,
    # but exp(Xb)=1, so need to remove cumulated ones from each time pts
    U <- UWZ[,1] - c(nt:1)
    W <- UWZ[,2]
    Z <- UWZ[,3]
    Z3 <- UWZ[,4]
    Z4 <- UWZ[,5]
    #d <- c(table(ipdata[ipdata$status==TRUE,c(1,2)]))
    ipdataE <- ipdata[ipdata$status==TRUE,]
    d<-c(aggregate(status~time,ipdataE, FUN=sum)$status)
    #X <- as.matrix(ipdata[ipdata$status==TRUE, -c(1,2)])
    X<-c(aggregate(x~time,ipdataE, FUN=mean)$x) #breslow's method for ties
    eXb <- c(exp(X*bbar))
    logL_D1 <- sum(d*X) - sum(d * W / U,na.rm = TRUE)
    logL_D2 <- sum(d * (W^2 - U*Z) / U^2)
    logL_D3 <- -sum(d * (Z3/U - 3*Z*W / (U^2) + 2*W^3/(U^3)))
    logL_D4 <- sum(d* (-Z4/U+4*Z3*W/(U^2)+3*(Z/U)^2-12*Z*W^2/(U^3)+6*W^4/(U^4)))
    logL <- sum(d*log(eXb/U))
  }
  derivatives <- list(logL=logL,logL_D1=logL_D1,logL_D2=logL_D2,logL_D3=logL_D3,logL_D4=logL_D4)
  return(derivatives)
}

computeLlDerivatives <- function(InputData, InputBeta){  # get pade approximation for each site
  InputData <- cbind.data.frame(InputData$time,InputData$y,InputData$x,InputData$stratumId)
  colnames(InputData) <- c("time","status","x","stratumId")
  nStrata <- max(InputData$stratumId)
  logL.deriv <- rep(0)
  for (StrataID in 1:nStrata){
    ipdata <- InputData[InputData$stratumId==StrataID,c("time","status","x")]
    logL.deriv<-logL.deriv+unlist(computeLlSingleDerivatives(InputBeta,ipdata))
  }
  return(logL.deriv)
}


# Borrowed from devtools:
# https://github.com/hadley/devtools/blob/ba7a5a4abd8258c52cb156e7b26bb4bf47a79f0b/R/utils.r#L44
is_installed <- function(pkg, version = 0) {
  installed_version <-
    tryCatch(
      utils::packageVersion(pkg),
      error = function(e) {
        NA
      }
    )
  !is.na(installed_version) && installed_version >= version
}

# Borrowed and adapted from devtools:
# https://github.com/hadley/devtools/blob/ba7a5a4abd8258c52cb156e7b26bb4bf47a79f0b/R/utils.r#L74
ensure_installed <- function(pkg) {
  if (!is_installed(pkg)) {
    msg <-
      paste0(sQuote(pkg), " must be installed for this functionality.")
    if (interactive()) {
      message(msg, "\nWould you like to install it?")
      if (menu(c("Yes", "No")) == 1) {
        install.packages(pkg)
      } else {
        stop(msg, call. = FALSE)
      }
    } else {
      stop(msg, call. = FALSE)
    }
  }
}

