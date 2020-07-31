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
#' Compute a  fixed-effect meta-analysis.
#'
#' @param data        A data frame containing either normal, skew-normal, custom parametric, or grid likelihood data. One row per database.
#' @param alpha       The alpha (expected type I error) used for the confidence intervals.
#' 
#' @examples 
#' populations <- simulatePopulations()
#' 
#' fitModelInDatabase <- function(population) {
#'   cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), 
#'                                             data = population, modelType = "cox")
#'   cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#'   approximation <-  approximateLikelihood(cyclopsFit, "x")
#'   return(approximation)
#' }
#' data <- lapply(populations, fitModelInDatabase)
#' data <- do.call("rbind", data)
#' computeFixedEffectMetaAnalysis(data)
#'
#' @export
computeFixedEffectMetaAnalysis <- function(data, alpha = 0.05) {
  # Determine type based on data structure:
  if ("logRr" %in% colnames(data)) {
    writeLines("Detected data following normal distribution")
    m <- meta::metagen(data$logRr, data$seLogRr, level.comb = 1 - alpha)
    ffx <- summary(m)$fixed
    estimate <- data.frame(rr = ffx$TE,
                           lb = ffx$lower,
                           ub = ffx$upper,
                           p = ffx$p,
                           logRr = log(ffx$TE),
                           seLogRr = ffx$seTE)
    return(estimate)
  } else if ("gamma" %in% colnames(data)) {
    writeLines("Detected data following custom parameric distribution")
    estimate <- computeEstimateFromCombiLl(data, alpha = alpha, fun = customFunction)
    return(estimate)
  } else if ("alpha" %in% colnames(data)) {
    writeLines("Detected data following skew normal distribution")
    estimate <- computeEstimateFromCombiLl(data, alpha = alpha, fun = skewNormal)
  } else if (is.list(data) && !is.data.frame(data)) {
    writeLines("Detected (pooled) patient-level data")
    population <- poolPopulations(data)
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
    mode <- coef(cyclopsFit)
    ci95 <- confint(cyclopsFit, 1, level = .95)
    estimate <- data.frame(rr = exp(mode),
                      lb = exp(ci95[2]),
                      ub = exp(ci95[3]),
                      logRr = mode,
                      seLogRr = (ci95[3] - ci95[2])/(2 * qnorm(0.975)))
    return(estimate)
  } else {
    writeLines("Detected data following grid distribution")
    x <- as.numeric(colnames(data))
    if (any(is.na(x))) {
      stop("Expecting grid data, but not all column names are numeric") 
    }
    estimate <- computeEstimateFromCombiGrids(data, alpha = alpha)
    return(estimate)
  }
}

combineLogLikelihoodFunctions <- function(x, fits, fun = customFunction) {
  ll <- apply(fits, 1, function(fit) fun(x, fit[1], fit[2], fit[3]))
  ll <- sum(ll)
  return(ll)
}

computeEstimateFromCombiLl <- function(fits, alpha = 0.05, fun = customFunction) {
  if (all(c("mu", "sigma", "gamma") %in% colnames(fits)))  {
    fits <- fits[, c("mu", "sigma", "gamma")]
  } else if (all(c("mu", "sigma", "alpha") %in% colnames(fits)))  {
    fits <- fits[, c("mu", "sigma", "alpha")]
  } else {
    stop("Expecting columns 'mu', 'sigma', and 'gamma' or 'alpha', but found columns '", paste(colnames(fits), collapse = "', '"), "'")
  }
  fit <- suppressWarnings(optim(0, function(x) -combineLogLikelihoodFunctions(x, fits, fun)))
  logRr <- fit$par           
  threshold <- -fit$value - qchisq(1 - alpha, df = 1) / 2
  
  precision = 1e-07
  
  # Binary search for upper bound
  L <- logRr
  H <- 10
  ub <- Inf
  while (H >= L) {
    M <- L + (H - L)/2
    llM <- combineLogLikelihoodFunctions(M, fits, fun)
    metric <- threshold - llM
    # writeLines(paste('M =', M, 'Metric = ',metric))
    if (metric > precision) {
      H <- M
    } else if (-metric > precision) {
      L <- M
    } else {
      ub <- M
      # print(M)
      break
    }
    if (M == logRr) {
      warning("Error finding upper bound")
      break
    } else if (M == 10) {
      warning("Confidence interval upper bound out of range")
      break
    }
  }
  
  # Binary search for lower bound  
  L <- -10
  H <- logRr
  lb <- -Inf
  while (H >= L) {
    M <- L + (H - L)/2
    llM <- combineLogLikelihoodFunctions(M, fits, fun)
    metric <- threshold - llM
    # writeLines(paste('M =', M, 'Metric = ',metric))
    if (metric > precision) {
      L <- M
    } else if (-metric > precision) {
      H <- M
    } else {
      lb <- M
      # print(M)
      break
    }
    if (M == logRr) {
      warning("Error finding lower bound")
      break
    } else if (M == -10) {
      warning("Confidence interval lower bound out of range")
      break
    }
  }
  result <- data.frame(rr = exp(logRr),
                       lb = exp(lb),
                       ub = exp(ub),
                       logRr = logRr,
                       seLogRr = (ub - lb)/(2*qnorm(0.975)))
  return(result)
}

poolPopulations <- function(populations) {
  highestId <- 0
  for (i in 1:length(populations)) {
    # Making sure stratum IDs are unique:
    populations[[i]]$stratumId <- populations[[i]]$stratumId + highestId
    highestId <- max(populations[[i]]$stratumId)
  }
  pooledPop <- do.call("rbind", populations) 
  return(pooledPop)
}
  
computeEstimateFromCombiGrids <- function(grids, alpha = 0.05) {
  grid <- apply(grids, 2, sum)
  maxIdx <- which(grid == max(grid))[1]
  logRr <- as.numeric(names(grid)[maxIdx])
  threshold <- grid[maxIdx] - qchisq(1 - alpha, df = 1) / 2
  lbIdx <- min(which(grid[1:maxIdx] > threshold))
  if (lbIdx == 1) {
    warning("Lower bound of confidence interval out of range")
  }
  lb <- as.numeric(names(grid)[lbIdx])
  ubIdx <- maxIdx + max(which(grid[(maxIdx + 1):length(grid)] > threshold))
  if (lbIdx == length(grid)) {
    warning("Upper bound of confidence interval out of range")
  }
  ub <- as.numeric(names(grid)[ubIdx])
  result <- data.frame(rr = exp(logRr),
                       lb = exp(lb),
                       ub = exp(ub),
                       logRr = logRr,
                       seLogRr = (ub - lb)/(2*qnorm(0.975)))
  return(result)
}
