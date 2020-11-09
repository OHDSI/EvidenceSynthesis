# Copyright 2020 Observational Health Data Sciences and Informatics
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

#' Approximate a likelihood function
#'
#' @description
#' Approximate the likelihood function using a parametric (normal, skew-normal, or custom parametric),
#' or grid approximation. The approximation does not reveal person-level information, and can
#' therefore be shared among data sites. When counts are low, a normal approximation might not be
#' appropriate.
#'
#' @param cyclopsFit      A model fitted using the [Cyclops::fitCyclopsModel()] function.
#' @param parameter       The parameter in the `cyclopsFit` object to profile.
#' @param approximation   The type of approximation. Valid options are `'normal'`, `'skew normal'`,
#'                        `'custom'`, or `'grid'`.
#' @param bounds          The bounds on the effect size used to fit the approximation.
#'
#' @seealso
#' [computeConfidenceInterval], [computeFixedEffectMetaAnalysis], [computeBayesianMetaAnalysis]
#'
#' @return
#' A vector of parameters of the likelihood approximation.
#'
#' @examples
#' # Simulate some data for this example:
#' populations <- simulatePopulations()
#'
#' cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'                                           data = populations[[1]],
#'                                           modelType = "cox")
#' cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#' approximation <- approximateLikelihood(cyclopsFit, "x")
#' approximation
#'
#' # (Estimates in this example will vary due to the random simulation)
#'
#' @export
approximateLikelihood <- function(cyclopsFit,
                                  parameter = 1,
                                  approximation = "custom",
                                  bounds = c(log(0.1), log(10))) {
  if (!approximation %in% c("normal", "skew normal", "custom", "grid"))
    stop("'approximation' argument should be 'normal', 'skew normal', 'custom', or 'grid'.")
  if (!is(cyclopsFit, "cyclopsFit"))
    stop("'cyclopsFit' argument should be of type 'cyclopsFit'")
  
  if (approximation == "grid") {
    x <- log(seq(exp(bounds[1]), exp(bounds[2]), by = 0.01))
    result <- getLikelihoodProfile(cyclopsFit, "x", x)
    names(result) <- x
    return(result)
  } else if (approximation == "normal") {
    if (cyclopsFit$return_flag != "SUCCESS") {
      return(data.frame(rr = NA, ci95Lb = NA, ci95Ub = NA, logRr = NA, seLogRr = NA))
    } else {
      mode <- coef(cyclopsFit)
      ci95 <- confint(cyclopsFit, 1, level = 0.95)
      return(data.frame(rr = exp(mode),
                        ci95Lb = exp(ci95[2]),
                        ci95Ub = exp(ci95[3]),
                        logRr = mode,
                        seLogRr = (ci95[3] - ci95[2])/(2 * qnorm(0.975))))
    }
  } else {
    x <- seq(bounds[1], bounds[2], length.out = 100)
    ll <- getLikelihoodProfile(cyclopsFit, parameter, x)
    if (approximation == "skew normal") {
      fun <- skewNormal
    } else {
      fun <- customFunction
    }
    result <- fitLogLikelihoodFunction(x, ll, fun = fun)
    if (approximation == "skew normal") {
      names(result)[names(result) == "gamma"] <- "alpha"
    }
    return(result)
  }
}

#' A custom function to approximate a log likelihood function
#'
#' @details
#' A custom parametric function designed to approximate the shape of the Cox log likelihood function.
#' When `gamma = 0` this function is the normal distribution.
#'
#' @param x       The log(hazard ratio) for which to approximate the log likelihood.
#' @param mu      The position parameter.
#' @param sigma   The scale parameter.
#' @param gamma   The skew parameter.
#'
#' @examples
#' customFunction(x = 0:3, mu = 0, sigma = 1, gamma = 0)
#'
#' @return
#' The approximate log likelihood for the given x.
#'
#' @export
customFunction <- function(x, mu, sigma, gamma) {
  return(((exp(gamma * (x - mu)))) * ((-(x - mu)^2)/(2 * sigma^2)))
  
  # Derivative: -(exp(gamma * (x - mu)) * (gamma * (x - mu) + 2) * (x - mu))/(2 * sigma^2)
}

#' The skew normal function to approximate a log likelihood function
#'
#' @details
#' The skew normal function. When `alpha = 0` this function is the normal distribution.
#'
#' @param x       The log(hazard ratio) for which to approximate the log likelihood.
#' @param mu      The position parameter.
#' @param sigma   The scale parameter.
#' @param alpha   The skew parameter.
#'
#' @examples
#' skewNormal(x = 0:3, mu = 0, sigma = 1, alpha = 0)
#'
#' @return
#' The approximate log likehood for the given x.
#'
#' @references
#' Azzalini, A. (2013). The Skew-Normal and Related Families. Institute of Mathematical Statistics
#' Monographs. Cambridge University Press.
#'
#' @export
skewNormal <- function(x, mu, sigma, alpha) {
  if (is.infinite(sigma)) {
    return(rep(0, length(x)))
  }
  return(log(2) + dnorm(x, mu, sigma, log = TRUE) + pnorm(alpha * (x - mu), 0, sigma, log.p = TRUE))
}

fitLogLikelihoodFunction <- function(beta, ll, weighByLikelihood = TRUE, fun = customFunction) {
  
  sumSquares <- function(p, maxAnchor = TRUE, idx = NULL) {
    approxLl <- fun(beta, p[1], p[2], p[3])
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
    return(data.frame(mu = 0, sigma = Inf, gamma = 0))
  }
  
  fit <- tryCatch({
    suppressWarnings(nlm(sumSquares, c(mode, 1, 0), maxAnchor = TRUE))
  }, error = function(e) {
    list(estimate = c(0, 0, 0), minimum = Inf)
  })
  result <- data.frame(mu = fit$estimate[1], sigma = fit$estimate[2], gamma = fit$estimate[3])
  minimum <- fit$minimum
  
  fit <- tryCatch({
    suppressWarnings(optim(c(mode, 1, 0), sumSquares, maxAnchor = TRUE))
  }, error = function(e) {
    list(par = c(0, 0, 0), value = Inf)
  })
  if (fit$value < minimum) {
    result <- data.frame(mu = fit$par[1], sigma = fit$par[2], gamma = fit$par[3])
    minimum <- fit$value
  }
  
  # Scale to standard (translate in log space so intersects at (0,0)):
  idx <- which(abs(beta) == min(abs(beta)))
  ll <- ll - ll[idx]
  fit <- tryCatch({
    suppressWarnings(nlm(sumSquares, c(mode, 1, 0), maxAnchor = FALSE, idx = idx))
  }, error = function(e) {
    list(estimate = c(0, 0, 0), minimum = Inf)
  })
  if (fit$minimum < minimum) {
    result <- data.frame(mu = fit$estimate[1], sigma = fit$estimate[2], gamma = fit$estimate[3])
    minimum <- fit$minimum
  }
  fit <- tryCatch({
    suppressWarnings(optim(c(mode, 1, 0), sumSquares, maxAnchor = FALSE, idx = idx))
  }, error = function(e) {
    list(par = c(0, 0, 0), value = Inf)
  })
  if (fit$value < minimum) {
    result <- data.frame(mu = fit$par[1], sigma = fit$par[2], gamma = fit$par[3])
    minimum <- fit$value
  }
  threshold <- 0.05
  if (minimum/length(beta) > threshold) {
    warn(paste("Mean squared error greater than ", threshold, ". Probably bad fit."))
  }
  result$sigma <- abs(result$sigma)
  return(result)
}

getLikelihoodProfile <- function(cyclopsFit, parameter, x) {
  ll <- Cyclops::getCyclopsProfileLogLikelihood(cyclopsFit, parameter, x)$value
  return(ll)
}

#' Compute the point estimate and confidence interval given a likelihood function approximation
#'
#' @details 
#' Compute the point estimate and confidence interval given a likelihood function approximation.
#'
#' @param approximation   An approximation of the likelihood function as fitted using the
#'                        [approximateLikelihood()] function.
#' @param alpha           The alpha (expected type I error).
#'
#' @return
#' A data frame containing the point estimate, and upper and lower bound of the confidence interval.
#' 
#' @examples
#' # Simulate some data for this example:
#' populations <- simulatePopulations()
#'
#' cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'                                           data = populations[[1]],
#'                                           modelType = "cox")
#' cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#' approximation <- approximateLikelihood(cyclopsFit, "x")
#' computeConfidenceInterval(approximation)
#' 
#' @export
computeConfidenceInterval <- function(approximation, alpha = 0.05) {
  # Determine type based on data structure:
  if ("logRr" %in% colnames(approximation)) {
    inform("Detected data following normal distribution")
    estimate <- data.frame(rr = exp(approximation$logRr),
                           lb = exp(approximation$logRr + qnorm(alpha/2) * approximation$seLogRr),
                           ub = exp(approximation$logRr + qnorm(1 - alpha/2) * approximation$seLogRr),
                           logRr = approximation$logRr,
                           seLogRr = approximation$seLogRr)
    return(estimate)
  } else if ("gamma" %in% colnames(approximation)) {
    inform("Detected data following custom parameric distribution")
    estimate <- computeEstimateFromApproximation(approximationFuntion = customFunction,
                                                 a = alpha,
                                                 mu = approximation$mu,
                                                 sigma = approximation$sigma,
                                                 gamma = approximation$gamma)
    return(estimate)
  } else if ("alpha" %in% colnames(approximation)) {
    inform("Detected data following skew normal distribution")
    estimate <- computeEstimateFromApproximation(approximationFuntion = skewNormal,
                                                 a = alpha,
                                                 mu = approximation$mu,
                                                 sigma = approximation$sigma,
                                                 alpha = approximation$alpha)
    return(estimate)
  } else {
    inform("Detected data following grid distribution")
    x <- as.numeric(colnames(approximation))
    if (any(is.na(x))) {
      abort("Expecting grid data, but not all column names are numeric")
    }
    estimate <- computeEstimateFromGrid(approximation, alpha = alpha)
    return(estimate)
  }
}