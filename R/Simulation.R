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

#' Create simulation settings
#'
#' @description
#' Create an object specifying a simulation. Currently only Cox proportional hazard models are
#' supported.
#'
#' @param nSites                Number of database sites to simulate.
#' @param n                     Number of subjects per site. Either a single number, or a vector of
#'                              length nSites.
#' @param treatedFraction       Fraction of subjects that is treated. Either a single number, or a
#'                              vector of length nSites.
#' @param nStrata               Number of strata per site. Either a single number, or a vector of
#'                              length nSites.
#' @param minBackgroundHazard   Minimum background hazard. Either a single number, or a vector of
#'                              length nSites.
#' @param maxBackgroundHazard   Maximum background hazard. Either a single number, or a vector of
#'                              length nSites.
#' @param hazardRatio           Hazard ratio.
#' @param randomEffectSd        Standard deviation of the log(hazardRatio). Fixed effect if equal to 0.
#' @param siteEffects           Fixed site effects (if assuming varying site-specific effects). Same effects if 0.
#'
#' @seealso
#' [simulatePopulations]
#'
#' @return
#' An object of type `simulationSettings`, to be used in the [simulatePopulations()] function.
#'
#' @examples
#' settings <- createSimulationSettings(nSites = 1, hazardRatio = 2)
#' populations <- simulatePopulations(settings)
#'
#' # Fit a Cox regression for the simulated data site:
#' cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'   data = populations[[1]],
#'   modelType = "cox"
#' )
#' cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#' coef(cyclopsFit)
#'
#' # (Estimates in this example will vary due to the random simulation)
#'
#' @export
createSimulationSettings <- function(nSites = 5,
                                     n = 10000,
                                     treatedFraction = 0.2,
                                     nStrata = 10,
                                     minBackgroundHazard = 2e-07,
                                     maxBackgroundHazard = 2e-05,
                                     hazardRatio = 2,
                                     randomEffectSd = 0,
                                     siteEffects = 0) {
  expand <- function(x) {
    if (length(x) == 1) {
      return(rep(x, nSites))
    } else {
      return(x)
    }
  }

  settings <- list(
    nSites = nSites,
    n = expand(n),
    treatedFraction = expand(treatedFraction),
    nStrata = expand(nStrata),
    minBackgroundHazard = expand(minBackgroundHazard),
    maxBackgroundHazard = expand(maxBackgroundHazard),
    hazardRatio = hazardRatio,
    randomEffectSd = randomEffectSd,
    siteEffects = expand(siteEffects)
  )
  class(settings) <- "simulationSettings"
  return(settings)
}

#' Simulate survival data for multiple databases
#'
#' @param settings   Either an object of type `simulationSettings`, created by the
#'                   [createSimulationSettings()] function or an object of type `sccsSimulationSettings`
#'                   as created by the [createSccsSimulationSettings()] function.
#'
#' @return
#' A object of class `simulation`, which is a list of population data frames. Depending on the type of
#' simulation, the data frames have different columns: Cox simulations will have the columns
#' `rowId`, `stratumId`, `x`, `time`, and `y`. SCCS simulations will have the columns `stratumId`, `a`,
#' `x1...xN`, `time`, and `y`.
#'
#' @examples
#' settings <- createSimulationSettings(nSites = 1, hazardRatio = 2)
#' populations <- simulatePopulations(settings)
#'
#' # Fit a Cox regression for the simulated data site:
#' cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'   data = populations[[1]],
#'   modelType = "cox"
#' )
#' cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#' coef(cyclopsFit)
#'
#' # (Estimates in this example will vary due to the random simulation)
#'
#' @export
simulatePopulations <- function(settings = createSimulationSettings()) {
  UseMethod("simulatePopulations", settings)
}

#' @export
simulatePopulations.simulationSettings <- function(settings = createSimulationSettings()) {
  stopifnot(class(settings) == "simulationSettings")
  thetas <- rnorm(
    n = settings$nSites,
    mean = log(settings$hazardRatio),
    sd = settings$randomEffectSd
  )
  thetas = thetas + settings$siteEffects
  hazardRatios <- exp(thetas)

  simulateSite <- function(i) {
    population <- data.frame(
      rowId = 1:settings$n[i],
      stratumId = round(runif(settings$n[i],
        min = 1,
        max = settings$nStrata[i]
      )),
      y = 0,
      x = as.numeric(runif(settings$n[i]) < settings$treatedFraction[i])
    )
    strataBackgroundHazard <- runif(settings$nStrata[i],
      min = settings$minBackgroundHazard[i],
      max = settings$maxBackgroundHazard[i]
    )
    population$hazard <- strataBackgroundHazard[population$stratumId]
    oldTotalHazard <- sum(population$hazard)
    population$hazard[population$x == 1] <- population$hazard[population$x == 1] * hazardRatios[i]

    # Normalize so higher hazard ratios don't come with more statistical power:
    newTotalHazard <- sum(population$hazard)
    population$hazard <- population$hazard * oldTotalHazard / newTotalHazard

    population$timeToOutcome <- 1 + round(rexp(n = settings$n[i], population$hazard))
    population$timeToCensor <- 1 + round(rexp(n = settings$n[i], 0.01))
    population$time <- population$timeToOutcome
    idx <- population$timeToCensor < population$timeToOutcome
    population$time[idx] <- population$timeToCensor[idx]
    population$y <- as.integer(!idx)
    return(population[, c("rowId", "stratumId", "x", "time", "y")])
  }
  simulation <- lapply(1:settings$nSites, simulateSite)
  attr(simulation, "simulationSettings") <- settings
  attr(simulation, "thetas") <- thetas
  class(simulation) <- "simulation"
  return(simulation)
}

#' @export
print.simulation <- function(x, ...) {
  writeLines("Simulation object")
  writeLines("")
  writeLines(paste("Number of sites: ", length(x)))
}

#' @export
summary.simulation <- function(object, ...) {
  summarizeSite <- function(site) {
    return(data.frame(
      subjects = nrow(site),
      treated = sum(site$x),
      outcomesTreated = sum(site$y[site$x ==
        1]), outcomesUntreated = sum(site$y[site$x == 0])
    ))
  }
  siteSummaries <- lapply(object, summarizeSite)
  siteSummaries <- do.call("rbind", siteSummaries)
  siteSummaries$theta <- attr(object, "thetas")
  class(siteSummaries) <- "summary.simulation"
  return(siteSummaries)
}

#' @export
print.summary.simulation <- function(x, ...) {
  class(x) <- "data.frame"
  rownames(x) <- 1:nrow(x)
  printCoefmat(x)
}


#' Simulate survival data across a federated data network, with negative control outcomes as well.
#'
#' @description
#' A function to simulate patient-level survival data for a hypothetical exposure, with simulated bias
#' and data-source-specific random effects. Patient-level data for negative control outcomes are simulated
#' as well to reflect systematic error.
#'
#' @param meanExposureEffect Average exposure effect; has to be on the log-scale
#' @param meanBias Average bias for the bias distribution
#' @param biasStd Standard deviation for the bias distribution
#' @param meanSiteEffect Average of the data-source-specific effects (typically should be zero)
#' @param siteEffectStd Standard deviation for data-source-specific effects
#' @param mNegativeControls Number of negative control outcomes
#' @param nSites Number of data sources
#' @param sitePop Population size per data source
#' @param seed Random seed
#' @param ... Arguments that will be passed to other functions
#'
#' @seealso [computeHierarchicalMetaAnalysis]
#'
#' @export
simulateMetaAnalysisWithNegativeControls <- function(meanExposureEffect = log(2),
                                                     meanBias = 0.5,
                                                     biasStd = 0.2,
                                                     meanSiteEffect = 0,
                                                     siteEffectStd = 0.1,
                                                     mNegativeControls = 10,
                                                     nSites = 10,
                                                     sitePop = 2000,
                                                     seed = 42,
                                                     ...){

  set.seed(seed)
  metaPopulations = list()

  # per outcome biases
  biases = rnorm(mNegativeControls, mean = meanBias, sd = biasStd)

  # per site variations
  siteEffects = rnorm(nSites, mean = meanSiteEffect, sd = siteEffectStd)

  # negative controls data
  for(i in 1:mNegativeControls){
    simulationSettings = createSimulationSettings(nSites = nSites,
                                                  n = sitePop,
                                                  hazardRatio = exp(biases[i]),
                                                  siteEffects = siteEffects,
                                                  ...)
    metaPopulations[[i]] = simulatePopulations(simulationSettings)
  }

  # main exposure data
  simulationSettings = createSimulationSettings(nSites = nSites,
                                                n = sitePop,
                                                hazardRatio = exp(meanExposureEffect + rnorm(1, meanBias, biasStd)),
                                                siteEffects = siteEffects,
                                                ...)
  metaPopulations[[mNegativeControls + 1]] = simulatePopulations(simulationSettings)


  hyperParameters = list(meanExposureEffect = meanExposureEffect,
                         meanBias = meanBias,
                         biasStd = biasStd,
                         meanSiteEffect = meanSiteEffect,
                         siteEffectStd = siteEffectStd,
                         mNegativeControls = mNegativeControls,
                         nSites = nSites,
                         sitePop = sitePop)

  attr(metaPopulations, "hyperParameters") = hyperParameters

  return(metaPopulations)
}

#' Create likelihood approximations from individual-trajectory data
#'
#' @param populations   Individual-level population data
#' @param approximation Type of approximation method
#'
#' @export
createApproximations <- function(populations, approximation) {
  UseMethod("createApproximations", populations)
}

#' @export
createApproximations.simulation <- function(populations, approximation) {

  fitModelInDatabase <- function(population, approximation) {
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                              data = population,
                                              modelType = "cox"
    )
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData,
                                           fixedCoefficients = c(!approximation %in% c("normal", "grid with gradients"))
    )
    res <- approximateLikelihood(cyclopsFit, "x",
                                 approximation = approximation)
    return(res)
  }
  data <- lapply(populations, fitModelInDatabase, approximation = approximation)
  if (!approximation %in% c("adaptive grid", "grid with gradients")) {
    data <- do.call("rbind", data)
  }
  return(data)
}

#' @export
createApproximations.sccsSimulation <- function(populations, approximation) {

  fitModelInDatabase <- function(population, approximation) {
    if (nrow(population) == 0) {
      return(NULL)
    }
    xColnames <- colnames(population)[grep("x[0-9]+", colnames(population))]
    formula <- as.formula(paste("y ~ a +", paste0(xColnames, collapse = " + "), " + strata(stratumId) + offset(log(time))"))
    cyclopsData <- Cyclops::createCyclopsData(formula, data = population, modelType = "cpr")
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData, fixedCoefficients = rep(!approximation %in% c("normal", "grid with gradients"), 1 + length(xColnames)))
    res <- EvidenceSynthesis::approximateLikelihood(cyclopsFit, "a", approximation = approximation)
    return(res)
  }
  data <- lapply(populations, fitModelInDatabase, approximation = approximation)
  if (approximation %in% c("adaptive grid", "grid with gradients")) {
    data <- data[!sapply(data, is.null)]
  } else {
    data <- do.call("rbind", data)
  }
  return(data)
}

#' Create SCCS simulation settings
#'
#' @description
#' Create an object specifying a simulation for the Self-Controlled Case Series (SCCS).
#'
#' @param nSites                Number of database sites to simulate.
#' @param n                     Number of subjects per site. Either a single number, or a vector of
#'                              length nSites.
#' @param atRiskTimeFraction    Fraction of patient time when at risk (exposed). Either a single number, or a
#'                              vector of length nSites.
#' @param timePartitions        Number of time partitions for seasonal covariates. Either a single number,
#'                              or a vector of length nSites.
#' @param timeCovariates        Number of covariates to represent seasonality. Either a single number,
#'                              or a vector of length nSites.
#' @param timeEffectSize        Strength of the seasonality effect. Either a single number,
#'                              or a vector of length nSites.
#' @param minBackgroundRate     Minimum background outcome rate. Either a single number, or a vector of
#'                              length nSites.
#' @param maxBackgroundRate     Maximum background outcome rate. Either a single number, or a vector of
#'                              length nSites.
#' @param rateRatio             The incidence rate ratio.
#' @param randomEffectSd        Standard deviation of the log(hazardRatio). Fixed effect if equal to 0.
#'
#' @seealso
#' [simulatePopulations]
#'
#' @return
#' An object of type `simulationSccsSettings`, to be used in the [simulatePopulations()] function.
#'
#' @examples
#' settings <- createSccsSimulationSettings(nSites = 1, rateRatio = 2)
#' populations <- simulatePopulations(settings)
#'
#' # Fit a SCCS regression for the simulated data site:
#' cyclopsData <- Cyclops::createCyclopsData(
#'   y ~ a + x1 + x2 + x3 + x4 + x5 + strata(stratumId) + offset(log(time)),
#'   data = populations[[1]],
#'   modelType = "cpr"
#' )
#' cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#' coef(cyclopsFit)
#'
#' # (Estimates in this example will vary due to the random simulation)
#'
#' @export
createSccsSimulationSettings <- function(nSites = 5,
                                         n = 10000,
                                         atRiskTimeFraction = 0.1,
                                         timePartitions = 24,
                                         timeCovariates = 5,
                                         timeEffectSize = log(2),
                                         minBackgroundRate = 0.001,
                                         maxBackgroundRate = 0.01,
                                         rateRatio = 2,
                                         randomEffectSd = 0) {
  expand <- function(x) {
    if (length(x) == 1) {
      return(rep(x, nSites))
    } else {
      return(x)
    }
  }
  settings <- list(
    nSites = nSites,
    n = expand(n),
    atRiskTimeFraction = expand(atRiskTimeFraction),
    timePartitions = expand(timePartitions),
    timeCovariates = expand(timeCovariates),
    timeEffectSize = expand(timeEffectSize),
    minBackgroundRate = expand(minBackgroundRate),
    maxBackgroundRate = expand(maxBackgroundRate),
    rateRatio = rateRatio,
    randomEffectSd = randomEffectSd
  )
  class(settings) <- "sccsSimulationSettings"
  return(settings)
}

#' @export
simulatePopulations.sccsSimulationSettings <- function(settings = createSimulationSettings()) {
  thetas <- rnorm(
    n = settings$nSites,
    mean = log(settings$rateRatio),
    sd = settings$randomEffectSd
  )
  rateRatios <- exp(thetas)

  simulateSite <- function(i) {
    personBackgroundRate <- runif(settings$n[i],
                                  min = settings$minBackgroundRate[i],
                                  max = settings$maxBackgroundRate[i]
    )
    if (settings $timeCovariates[i] == 0) {
      designMatrix <- NULL
    } else {
      seasonKnots <- 0.5 + seq(0, 12, length.out = settings$timeCovariates[i] + 1)
      designMatrix <- cyclicSplineDesign(1:12, seasonKnots)
      while(nrow(designMatrix) < settings$timePartitions[i]) {
        designMatrix <- rbind(designMatrix, designMatrix)
      }
      designMatrix <- designMatrix[1:settings$timePartitions[i], ]
    }

    # Simulate time periods with covariates
    population <- data.frame(
      stratumId = rep(seq_len(settings$n[i]), each = settings$timePartitions[i]),
      time = 1/settings$timePartitions[i],
      a = 0,
      deltaA = 0,
      backgroundRate = rep(personBackgroundRate, each = settings$timePartitions[i])
    )
    if (settings $timeCovariates[i] != 0) {
      x <- do.call(rbind, replicate(settings$n[i], designMatrix, simplify = FALSE))
      population <- cbind(population, x)
      colnames(population)[-(1:5)] <- paste0("x", seq_len(settings$timeCovariates[i]))
    }

    # Simulate time at risk and merge with time periods
    tarStarts <- runif(settings$n[i], min = 0, max = 1 - settings$atRiskTimeFraction[i])
    tarEnds <- tarStarts + settings$atRiskTimeFraction[i]
    tarStartIdx <- (tarStarts + seq_len(settings$n[i]) - 1) * settings$timePartitions[i] + 1
    tarEndIdx <- (tarEnds + seq_len(settings$n[i]) - 1) * settings$timePartitions[i] + 1
    prePartitionTime <- (tarStartIdx - floor(tarStartIdx)) / settings$timePartitions[i]
    postPartitionTime <- (ceiling(tarEndIdx) - tarEndIdx) / settings$timePartitions[i]
    tarStartIdx <- floor(tarStartIdx)
    tarEndIdx <- floor(tarEndIdx)
    population$deltaA[tarStartIdx] <- 1
    population$deltaA[tarEndIdx] <- population$deltaA[tarEndIdx] - 1
    population$a <- cumsum(population$deltaA)
    population$deltaA <- NULL
    prePartitions <- population[tarStartIdx, ]
    postPartitions <- population[tarEndIdx, ]
    population$time[tarStartIdx] <- population$time[tarStartIdx] - prePartitionTime
    population$a[tarStartIdx] <- 1
    prePartitions$time <- prePartitionTime
    prePartitions$a <- 0
    population$time[tarEndIdx] <- population$time[tarEndIdx] - postPartitionTime
    population$a[tarEndIdx] <- 1
    postPartitions$time <- postPartitionTime
    postPartitions$a <- 0
    population <- rbind(population, prePartitions, postPartitions)

    # Simulate outcomes
    timeEffect <- runif(settings$timeCovariates[i], max = settings$timeEffectSize[i])
    population$rate <- population$backgroundRate * exp(population$a * thetas[i] + rowSums(population[, ncol(population)-(seq_len(settings$timeCovariates[i])) + 1] * timeEffect))

    # Normalize so higher hazard ratios don't come with more statistical power:
    population$totalRate <- population$rate * population$time
    totalRates <- aggregate(totalRate ~ stratumId, data = population, sum)
    population$totalRate <- NULL
    population <- merge(population, totalRates, by = "stratumId")
    population$rate <- population$rate * population$backgroundRate / population$totalRate

    population$y <- rpois(nrow(population), population$rate * population$time)
    population$backgroundRate <- NULL
    population$rate <- NULL

    # Remove subjects with no events
    stratumIds <- unique(population$stratumId[population$y > 0])
    population <- population[population$stratumId %in% stratumIds, ]
    population$rowId <- seq_len(nrow(population))

    return(population)
  }
  simulation <- lapply(1:settings$nSites, simulateSite)
  attr(simulation, "sccsSimulationSettings") <- settings
  attr(simulation, "thetas") <- thetas
  class(simulation) <- "sccsSimulation"
  return(simulation)
}

#' @export
print.sccsSimulation <- function(x, ...) {
  writeLines("SCCS simulation object")
  writeLines("")
  writeLines(paste("Number of sites: ", length(x)))
}

#' @export
summary.sccsSimulation<- function(object, ...) {
  summarizeSite <- function(site) {
    return(data.frame(
      cases = length(unique(site$stratumId)),
      exposedCases = length(unique(site$stratumId[site$a == 1 & site$y > 0]))
    ))
  }
  siteSummaries <- lapply(object, summarizeSite)
  siteSummaries <- do.call("rbind", siteSummaries)
  siteSummaries$theta <- attr(object, "thetas")
  class(siteSummaries) <- "summary.sccsSimulation"
  return(siteSummaries)
}

#' @export
print.summary.sccsSimulation <- function(x, ...) {
  class(x) <- "data.frame"
  rownames(x) <- 1:nrow(x)
  printCoefmat(x)
}

cyclicSplineDesign <- function(x, knots, ord = 4) {
  nk <- length(knots)
  if (ord < 2) {
    stop("order too low")
  }
  if (nk < ord) {
    stop("too few knots")
  }
  knots <- sort(knots)
  k1 <- knots[1]
  if (min(x) < k1 || max(x) > knots[nk]) {
    stop("x out of range")
  }
  xc <- knots[nk - ord + 1]
  knots <- c(k1 - (knots[nk] - knots[(nk - ord + 1):(nk - 1)]), knots)
  ind <- x > xc
  X1 <- splines::splineDesign(knots, x, ord, outer.ok = TRUE)
  x[ind] <- x[ind] - max(knots) + k1
  if (sum(ind)) {
    X2 <- splines::splineDesign(knots, x[ind], ord, outer.ok = TRUE)
    X1[ind, ] <- X1[ind, ] + X2
  }
  X1
}
