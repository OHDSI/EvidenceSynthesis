# Copyright 2023 Observational Health Data Sciences and Informatics
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
#'     data = population,
#'     modelType = "cox"
#'   )
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
  type <- detectApproximationType(data)
  data <- cleanApproximations(data)
  if (type == "normal") {
    m <- meta::metagen(
      TE = data$logRr,
      seTE = data$seLogRr,
      studlab = rep("", nrow(data)),
      byvar = NULL,
      sm = "RR",
      level.comb = 1 - alpha
    )
    ffx <- summary(m)$fixed
    estimate <- data.frame(
      rr = exp(ffx$TE),
      lb = exp(ffx$lower),
      ub = exp(ffx$upper),
      logRr = ffx$TE,
      seLogRr = ffx$seTE
    )
    return(estimate)
  } else if (type == "custom") {
    estimate <- computeEstimateFromApproximation(
      approximationFuntion = combineLogLikelihoodFunctions,
      a = alpha,
      fits = data,
      fun = customFunction
    )
    return(estimate)
  } else if (type == "skew normal") {
    estimate <- computeEstimateFromApproximation(
      approximationFuntion = combineLogLikelihoodFunctions,
      a = alpha,
      fits = data,
      fun = skewNormal
    )
    return(estimate)
  } else if (type == "pooled") {
    population <- poolPopulations(data)
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
      data = population,
      modelType = "cox"
    )
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
    mode <- coef(cyclopsFit)
    ci95 <- confint(cyclopsFit, 1, level = 0.95)
    estimate <- data.frame(
      rr = exp(mode),
      lb = exp(ci95[2]),
      ub = exp(ci95[3]),
      logRr = mode,
      seLogRr = (ci95[3] - ci95[2]) / (2 * qnorm(0.975))
    )
    return(estimate)
  } else if (type == "adaptive grid") {
    estimate <- computeFixedEffectAdaptiveGrid(data, alpha)
    return(estimate)
  } else if (type == "grid") {
    data <- cleanData(data,
      colnames(data),
      minValues = rep(-1e6, ncol(data)),
      maxValues = rep(0, ncol(data)),
      grid = TRUE
    )
    grid <- apply(data, 2, sum)
    estimate <- computeEstimateFromGrid(grid, alpha = alpha)
    return(estimate)
  } else {
    abort(sprintf("Approximation type '%s' not supported by this function", type))
  }
}

estimate <- computeFixedEffectAdaptiveGrid <- function(data, alpha) {
  gridPoints <- sort(unique(do.call(c, lapply(data, function(x) x$point))))
  values <- rep(0, length(gridPoints))
  for (i in seq_along(data)) {
    cleanedData <- as.data.frame(data[[i]])
    cleanedData$value <- cleanedData$value - max(cleanedData$value)
    cleanedData <- cleanData(cleanedData,
      c("point", "value"),
      minValues = c(-100, -1e6),
      maxValues = c(100, 0)
    )
    cleanedData <- cleanedData[order(cleanedData$point), ]

    # Compute likelihood at unioned grid points, using linear interpolation where needed:
    cursor <- 1
    for (j in seq_along(gridPoints)) {
      if (cursor == 1 && gridPoints[j] <= cleanedData$point[cursor]) {
        value <- cleanedData$value[cursor]
      } else if (cursor == nrow(cleanedData) && gridPoints[j] >= cleanedData$point[cursor]) {
        value <- cleanedData$value[cursor]
      } else {
        if (gridPoints[j] > cleanedData$point[cursor]) {
          cursor <- cursor + 1
        }
        if (gridPoints[j] == cleanedData$point[cursor]) {
          value <- cleanedData$value[cursor]
        } else {
          x1 <- cleanedData$point[cursor - 1]
          x2 <- cleanedData$point[cursor]
          y1 <- cleanedData$value[cursor - 1]
          y2 <- cleanedData$value[cursor]
          value <- y1 + ((y2 - y1) * (gridPoints[j] - x1) / (x2 - x1))
        }
      }
      values[j] <- values[j] + value
    }
  }
  names(values) <- gridPoints
  estimate <- computeEstimateFromGrid(grid = values, alpha = alpha)
  return(estimate)
}

combineLogLikelihoodFunctions <- function(x, fits, fun = customFunction) {
  ll <- apply(fits, 1, function(fit) fun(x, fit[1], fit[2], fit[3]))
  ll <- sum(ll)
  return(ll)
}

computeEstimateFromApproximation <- function(approximationFuntion, a = 0.05, ...) {
  fit <- suppressWarnings(optim(0, function(x) -approximationFuntion(x, ...)))
  logRr <- fit$par
  threshold <- -fit$value - qchisq(1 - a, df = 1) / 2

  precision <- 1e-07

  # Binary search for upper bound
  L <- logRr
  H <- 10
  ub <- Inf
  while (H >= L) {
    M <- L + (H - L) / 2
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
    M <- L + (H - L) / 2
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
  result <- data.frame(
    rr = exp(logRr),
    lb = exp(lb),
    ub = exp(ub),
    logRr = logRr,
    seLogRr = (ub - lb) / (2 * qnorm(0.975))
  )
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
  if (maxIdx == 1 || maxIdx == length(grid)) {
    warn("Point estimate out of range")
    result <- data.frame(
      rr = NA,
      lb = NA,
      ub = NA,
      logRr = NA,
      seLogRr = NA
    )
    return(result)
  }

  logRr <- as.numeric(names(grid)[maxIdx])
  threshold <- unlist(grid[maxIdx]) - qchisq(1 - alpha, df = 1) / 2
  lbIdx <- min(which(grid[1:maxIdx] > threshold))
  # TODO: use interpolation to get more precise CI
  if (lbIdx == 1) {
    warn("Lower bound of confidence interval out of range")
  }
  lb <- as.numeric(names(grid)[lbIdx])
  ubIdx <- maxIdx + max(which(grid[(maxIdx + 1):length(grid)] > threshold))
  if (lbIdx == length(grid)) {
    warn("Upper bound of confidence interval out of range")
  }
  ub <- as.numeric(names(grid)[ubIdx])
  result <- data.frame(
    rr = exp(logRr),
    lb = exp(lb),
    ub = exp(ub),
    logRr = logRr,
    seLogRr = (ub - lb) / (2 * qnorm(0.975))
  )
  return(result)
}
