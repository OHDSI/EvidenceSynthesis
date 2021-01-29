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
#'                                           data = populations[[1]],
#'                                           modelType = "cox")
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
                                     randomEffectSd = 0) {
  expand <- function(x) {
    if (length(x) == 1) {
      return(rep(x, nSites))
    } else {
      return(x)
    }
  }
  
  settings <- list(nSites = nSites,
                   n = expand(n),
                   treatedFraction = expand(treatedFraction),
                   nStrata = expand(nStrata),
                   minBackgroundHazard = expand(minBackgroundHazard),
                   maxBackgroundHazard = expand(maxBackgroundHazard),
                   hazardRatio = hazardRatio,
                   randomEffectSd = randomEffectSd)
  class(settings) <- "simulationSettings"
  return(settings)
}

#' Simulate survival data for multiple databases
#'
#' @param settings   An object of type `simulationSettings`, created by the
#'                   [createSimulationSettings()] function.
#'
#' @return
#' A object of class `simulation`, which is a list of populations, each a data frame with columns
#' `rowId`, `stratumId`, `x`, `time`, and `y`.
#'
#' @examples
#' settings <- createSimulationSettings(nSites = 1, hazardRatio = 2)
#' populations <- simulatePopulations(settings)
#'
#' # Fit a Cox regression for the simulated data site:
#' cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'                                           data = populations[[1]],
#'                                           modelType = "cox")
#' cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#' coef(cyclopsFit)
#'
#' # (Estimates in this example will vary due to the random simulation)
#'
#' @export
simulatePopulations <- function(settings = createSimulationSettings()) {
  stopifnot(class(settings) == "simulationSettings")
  thetas <- rnorm(n = settings$nSites,
                  mean = log(settings$hazardRatio),
                  sd = settings$randomEffectSd)
  hazardRatios <- exp(thetas)
  
  simulateSite <- function(i) {
    population <- data.frame(rowId = 1:settings$n[i],
                             stratumId = round(runif(settings$n[i],
                                                     min = 1,
                                                     max = settings$nStrata[i])),
                             y = 0,
                             x = as.numeric(runif(settings$n[i]) < settings$treatedFraction[i]))
    strataBackgroundHazard <- runif(settings$nStrata[i],
                                    min = settings$minBackgroundHazard[i],
                                    max = settings$maxBackgroundHazard[i])
    population$hazard <- strataBackgroundHazard[population$stratumId]
    oldTotalHazard <- sum(population$hazard)
    population$hazard[population$x == 1] <- population$hazard[population$x == 1] * hazardRatios[i]
    
    # Normalize so higher hazard ratios don't come with more statistical power:
    newTotalHazard <- sum(population$hazard)
    population$hazard <- population$hazard * oldTotalHazard/newTotalHazard
    
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
    return(data.frame(subjects = nrow(site),
                      treated = sum(site$x),
                      outcomesTreated = sum(site$y[site$x ==
                                                     1]), outcomesUntreated = sum(site$y[site$x == 0])))
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
