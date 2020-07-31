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
#' Approximate the likelihood function using a parametric (normal, skew-normal, or custom parametric), or grid 
#' approximation. The approximation does not reveal person-level information, and can therefore be shared amongst
#' data sites. When counts are low, a normal approximation might not be appropriate. 
#'
#' @param cyclopsFit    A model fitted using the [Cyclops::fitCyclopsModel()] function.
#' @param parameter     The parameter in the `cyclopsFit` object to profile.
#' @param approximation The type of approximation. Valid options are `'normal'`, `'skew normal'`, `'custom'`, or `'grid'`.
#' @param bounds        The bounds on the effect size used to fit the approximation. 
#'
#' @return
#' A vector of parameters of the likelihood approximation.
#' 
#' @examples
#' populations <- simulatePopulations()
#' cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), 
#'                                           data = populations[[1]], 
#'                                           modelType = "cox")
#' cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#' approximation <-  approximateLikelihood(cyclopsFit, "x")
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
    result <- Cyclops::getCyclopsProfileLogLikelihood(cyclopsFit, "x", x)$value
    names(result) <- x
    return(result)
  } else if (approximation == "normal") {
    if (cyclopsFit$return_flag != "SUCCESS") {
      return(data.frame(rr = NA,
                        ci95Lb = NA,
                        ci95Ub = NA,
                        logRr = NA,
                        seLogRr = NA))
    } else {
      mode <- coef(cyclopsFit)
      ci95 <- confint(cyclopsFit, 1, level = .95)
      return(data.frame(rr = exp(mode),
                        ci95Lb = exp(ci95[2]),
                        ci95Ub = exp(ci95[3]),
                        logRr = mode,
                        seLogRr = (ci95[3] - ci95[2])/(2 * qnorm(0.975))))
    }
  } else {
    x <- seq(bounds[1], bounds[2], length.out = 100)
    ll <- Cyclops::getCyclopsProfileLogLikelihood(cyclopsFit, parameter, x)$value
    if (approximation == "skew normal") {
      fun = skewNormal
    } else {
      fun = customFunction
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
#' @param x       The log(hazard ratio) for which to approximate the log likelihood.
#' @param mu      The position parameter.
#' @param sigma   The scale parameter.
#' @param gamma   The skew parameter.
#'
#' @return The approximate log likelihood for the given x.
#' 
#' @export
customFunction <- function(x, mu, sigma, gamma) {
  return(((exp(gamma*(x - mu)))) * ((-(x - mu)^2)/(2*sigma^2)))
  
  # Derivative:
  #  -(exp(gamma * (x - mu)) * (gamma * (x - mu) + 2) * (x - mu))/(2 * sigma^2) 
}

#' The skew normal function to approximate a log likelihood function
#'
#' @param x       The log(hazard ratio) for which to approximate the log likelihood.
#' @param mu      The position parameter.
#' @param sigma   The scale parameter.
#' @param alpha   The skew parameter.
#'
#' @return The approximate log likehood for the given x.
#' 
#' @export
skewNormal <- function(x, mu, sigma, alpha) {
  return(log(2) + dnorm(x, mu, sigma, log = TRUE) + pnorm(alpha*(x - mu), 0, sigma, log.p = TRUE))
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
  
  # Scale to standard (translate in log space so max at 0):
  ll <- ll - max(ll)
  weights <- exp(ll)
  # Weights shouldn't be too small:
  weights[weights < 1e-5] <- 1e-5

  mode <- beta[ll == 0][1]
  if (mode == min(beta) || mode == max(beta)) {
    mode = 0
  }
  
  fit <- tryCatch({
    suppressWarnings(nlm(sumSquares, c(mode,1,0), maxAnchor = TRUE))
  }, error = function(e){
    list(estimate = c(0,0,0), minimum = Inf)
  })
  result <- data.frame(mu = fit$estimate[1],
                       sigma  = fit$estimate[2],
                       gamma = fit$estimate[3])
  minimum <- fit$minimum
  
  fit <- tryCatch({
    suppressWarnings(optim(c(mode,1,0), sumSquares, maxAnchor = TRUE))
  }, error = function(e){
    list(par = c(0,0,0), value = Inf)
  })
  if (fit$value < minimum) {
    result <- data.frame(mu = fit$par[1],
                         sigma  = fit$par[2],
                         gamma = fit$par[3])
    minimum <- fit$value
  }
  
  # Scale to standard (translate in log space so intersects at (0,0)):
  idx <- which(abs(beta) == min(abs(beta)))
  ll <- ll - ll[idx]
  fit <- tryCatch({
    suppressWarnings(nlm(sumSquares, c(mode,1,0), maxAnchor = FALSE, idx = idx))
  }, error = function(e){
    list(estimate = c(0,0,0), minimum = Inf)
  })
  if (fit$minimum < minimum) {
    result <- data.frame(mu = fit$estimate[1],
                         sigma  = fit$estimate[2],
                         gamma = fit$estimate[3])
    minimum <- fit$minimum
  }
  fit <- tryCatch({
    suppressWarnings(optim(c(mode,1,0), sumSquares, maxAnchor = FALSE, idx = idx))
  }, error = function(e){
    list(par = c(0,0,0), value = Inf)
  })
  if (fit$value < minimum) {
    result <- data.frame(mu = fit$par[1],
                         sigma  = fit$par[2],
                         gamma = fit$par[3])
    minimum <- fit$value
  }
  threshold <- 0.05
  if (minimum / length(beta) > threshold) {
    warning("Mean squared error greater than ", threshold, ". Probably bad fit.")
  }
  result$sigma <- abs(result$sigma)
  return(result)
}
