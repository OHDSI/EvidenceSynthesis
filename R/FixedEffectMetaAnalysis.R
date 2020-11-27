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

#' Compute a fixed-effect meta-analysis
#'
#' @description
#' Compute a fixed-effect meta-analysis using a choice of various likelihood approximations.
#'
#' @param data    A data frame containing either normal, skew-normal, custom parametric, or grid
#'                likelihood data. One row per database.
#' @param alpha   The alpha (expected type I error) used for the confidence intervals.
#'
#' @seealso
#' [approximateLikelihood], [computeBayesianMetaAnalysis]
#'
#' @return
#' The meta-analytic estimate, expressed as the point estimate hazard ratio (rr), its 95 percent
#' confidence interval (lb, ub), as well as the log of the point estimate (logRr), and the standard
#' error (seLogRr).
#'
#' @examples
#' # Simulate some data for this example:
#' populations <- simulatePopulations()
#'
#' # Fit a Cox regression at each data site, and approximate likelihood function:
#' fitModelInDatabase <- function(population) {
#'   cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'                                             data = population,
#'                                             modelType = "cox")
#'   cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#'   approximation <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "custom")
#'   return(approximation)
#' }
#' approximations <- lapply(populations, fitModelInDatabase)
#' approximations <- do.call("rbind", approximations)
#'
#' # At study coordinating center, perform meta-analysis using per-site approximations:
#' computeFixedEffectMetaAnalysis(approximations)
#'
#' # (Estimates in this example will vary due to the random simulation)
#'
#' @export
computeFixedEffectMetaAnalysis <- function(data, alpha = 0.05) {
  # Determine type based on data structure:
  if ("logRr" %in% colnames(data)) {
    inform("Detected data following normal distribution")
    data <- cleanData(data, c("logRr", "seLogRr"), minValues = c(-100, 1e-05))
    m <- meta::metagen(TE = data$logRr,
                       seTE = data$seLogRr,
                       studlab = rep("", nrow(data)),
                       byvar = NULL,
                       sm = "RR",
                       level.comb = 1 - alpha)
    ffx <- summary(m)$fixed
    estimate <- data.frame(rr = exp(ffx$TE),
                           lb = exp(ffx$lower),
                           ub = exp(ffx$upper),
                           logRr = ffx$TE,
                           seLogRr = ffx$seTE)
    return(estimate)
  } else if ("gamma" %in% colnames(data)) {
    inform("Detected data following custom parameric distribution")
    data <- cleanData(data, c("mu", "sigma", "gamma"), minValues = c(-100, 1e-05, -100))
    estimate <- computeEstimateFromApproximation(approximationFuntion = combineLogLikelihoodFunctions, 
                                                 a = alpha, 
                                                 fits = data,
                                                 fun = customFunction)
    return(estimate)
  } else if ("alpha" %in% colnames(data)) {
    inform("Detected data following skew normal distribution")
    data <- cleanData(data,
                      c("mu", "sigma", "alpha"),
                      minValues = c(-100, 1e-05, -10000),
                      maxValues = c(100, 10000, 10000))
    estimate <- computeEstimateFromApproximation(approximationFuntion = combineLogLikelihoodFunctions, 
                                                 a = alpha, 
                                                 fits = data,
                                                 fun = skewNormal)
    return(estimate)
  } else if (is.list(data) && !is.data.frame(data)) {
    inform("Detected (pooled) patient-level data")
    population <- poolPopulations(data)
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                              data = population,
                                              modelType = "cox")
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
    mode <- coef(cyclopsFit)
    ci95 <- confint(cyclopsFit, 1, level = 0.95)
    estimate <- data.frame(rr = exp(mode),
                           lb = exp(ci95[2]),
                           ub = exp(ci95[3]),
                           logRr = mode,
                           seLogRr = (ci95[3] - ci95[2])/(2 * qnorm(0.975)))
    return(estimate)
  } else {
    inform("Detected data following grid distribution")
    x <- as.numeric(colnames(data))
    if (any(is.na(x))) {
      abort("Expecting grid data, but not all column names are numeric")
    }
    grid <- apply(data, 2, sum)
    estimate <- computeEstimateFromGrid(grid, alpha = alpha)
    return(estimate)
  }
}

combineLogLikelihoodFunctions <- function(x, fits, fun = customFunction) {
  ll <- apply(fits, 1, function(fit) fun(x, fit[1], fit[2], fit[3]))
  ll <- sum(ll)
  return(ll)
}

computeEstimateFromApproximation <- function(approximationFuntion, a = 0.05, ...) {
  fit <- suppressWarnings(optim(0, function(x) -approximationFuntion(x, ...)))
  logRr <- fit$par
  threshold <- -fit$value - qchisq(1 - a, df = 1)/2

  precision <- 1e-07

  # Binary search for upper bound
  L <- logRr
  H <- 10
  ub <- Inf
  while (H >= L) {
    M <- L + (H - L)/2
    llM <- approximationFuntion(M, ...)
    metric <- threshold - llM
    if (metric > precision) {
      H <- M
    } else if (-metric > precision) {
      L <- M
    } else {
      ub <- M
      break
    }
    if (M == logRr) {
      warn("Error finding upper bound")
      break
    } else if (M == 10) {
      warn("Confidence interval upper bound out of range")
      break
    }
  }

  # Binary search for lower bound
  L <- -10
  H <- logRr
  lb <- -Inf
  while (H >= L) {
    M <- L + (H - L)/2
    llM <- approximationFuntion(M, ...)
    metric <- threshold - llM
    if (metric > precision) {
      L <- M
    } else if (-metric > precision) {
      H <- M
    } else {
      lb <- M
      break
    }
    if (M == logRr) {
      warn("Error finding lower bound")
      break
    } else if (M == -10) {
      warn("Confidence interval lower bound out of range")
      break
    }
  }
  result <- data.frame(rr = exp(logRr),
                       lb = exp(lb),
                       ub = exp(ub),
                       logRr = logRr,
                       seLogRr = (ub - lb)/(2 * qnorm(0.975)))
  return(result)
}

poolPopulations <- function(populations) {
  highestId <- 0
  for (i in 1:length(populations)) {
    populations[[i]]$stratumId <- populations[[i]]$stratumId + highestId
    highestId <- max(populations[[i]]$stratumId) + 1
  }
  pooledPop <- do.call("rbind", populations)
  return(pooledPop)
}

computeEstimateFromGrid <- function(grid, alpha = 0.05) {
  maxIdx <- which(grid == max(grid))[1]
  logRr <- as.numeric(names(grid)[maxIdx])
  threshold <- unlist(grid[maxIdx]) - qchisq(1 - alpha, df = 1)/2
  lbIdx <- min(which(grid[1:maxIdx] > threshold))
  if (lbIdx == 1) {
    warn("Lower bound of confidence interval out of range")
  }
  lb <- as.numeric(names(grid)[lbIdx])
  ubIdx <- maxIdx + max(which(grid[(maxIdx + 1):length(grid)] > threshold))
  if (lbIdx == length(grid)) {
    warn("Upper bound of confidence interval out of range")
  }
  ub <- as.numeric(names(grid)[ubIdx])
  result <- data.frame(rr = exp(logRr),
                       lb = exp(lb),
                       ub = exp(ub),
                       logRr = logRr,
                       seLogRr = (ub - lb)/(2 * qnorm(0.975)))
  return(result)
}
