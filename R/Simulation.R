# Copyright 2019 Observational Health Data Sciences and Informatics
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
#' An object specifing a simulation. Currently only Cox proportional hazard models are supported.
#'
#' @param nSites          Number of database sites to simulate.
#' @param n               Number of subjects per site. Either a single number, or a vector of length nSites.
#' @param treatedFraction Fraction of subjects that is treated. Either a single number, or a vector of length nSites.
#' @param nStrata         Number of strata per site. Either a single number, or a vector of length nSites.
#' @param minBackgroundHazard  Minimum background hazard. Either a single number, or a vector of length nSites.
#' @param maxBackgroundHazard  Maximum background hazard. Either a single number, or a vector of length nSites.
#' @param hazardRatio     Hazard ratio. 
#' @param randomEffectSd  Standard deviation of the log(hazardRatio). Fixed effect if equal to 0.
#'
#' @return
#' An object of type \code{simulationSettings}.
#' 
#' @export
createSimulationSettings <- function(nSites = 5,
                                     n = 10000,
                                     treatedFraction = 0.2,
                                     nStrata = 10,
                                     minBackgroundHazard = 0.0000002,
                                     maxBackgroundHazard = 0.00002,
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

#' Run simulation
#'
#' @param settings   An object of type \code{simulationSettings}, created by the 
#'                   \code{\link{createSimulationSettings}} function.
#'
#' @return
#' A object of class \code{simulation}..
#'
#' @examples
#' settings <- createSimulationSettings()
#' simulation <- runSimulation(settings)
#' 
#' @export
runSimulation <- function(settings) {
  stopifnot(class(settings) == "simulationSettings")
  hazardRatios <- exp(rnorm(n = settings$nSites,
                            mean = log(settings$hazardRatio),
                            sd =  settings$randomEffectSd))

  simulateSite <- function(i) {
    population <- data.frame(rowId = 1:settings$n[i],
                             stratumId = round(runif(settings$n[i] ,min = 1, max = settings$nStrata[i])),
                             y = 0,
                             x = as.numeric(runif(settings$n[i]) < settings$treatedFraction[i]))
    strataBackgroundHazard <- runif(settings$nStrata[i], 
                                    min = settings$minBackgroundHazard[i], 
                                    max = settings$maxBackgroundHazard[i])
    population$hazard <-  strataBackgroundHazard[population$stratumId]
    population$hazard[population$x == 1] <- population$hazard[population$x == 1] * hazardRatios[i]
    population$timeToOutcome <- 1 + round(rexp(n = settings$n[i], population$hazard))
    population$timeToCensor <- 1 + round(runif(n = settings$n[i], min = 0, max = 499))
    population$time <- population$timeToOutcome
    population$time[population$timeToCensor < population$timeToOutcome] <- population$timeToCensor[population$timeToCensor < population$timeToOutcome]
    population$y <- as.integer(population$timeToCensor > population$timeToOutcome)
    return(population[, c("rowId", "stratumId", "x", "time", "y")])
  }
  simulation <- lapply(1:settings$nSites, simulateSite)
  attr(simulation, "simulationSettings") <- settings
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
                        outcomesTreated = sum(site$y[site$x == 1]),
                        outcomesUntreated = sum(site$y[site$x == 0])))
  }
  siteSummaries <- lapply(object, summarizeSite)
  siteSummaries <- do.call("rbind", siteSummaries)
  class(siteSummaries) <- "summary.simulation"
  return(siteSummaries)
}

#' @export
print.summary.simulation <- function(x, ...) {
  class(x) <- "data.frame"
  rownames(x) <- 1:nrow(x)
  printCoefmat(x)
}

#' Estimate effect by pooling data
#' 
#' @description 
#' Estimate the hazard ratio across site by pooling the data.
#'
#' @param simulation An object of type \code{simulation} as created by the 
#'                   \code{\link{runSimulation}} function.
#' @param useCyclops Should Cyclops be used? Cyclops doesn't assume a normal distribution for the 
#'                   likelihood function, and uses likelihood profiling instead.
#'
#' @return
#' A data frame with one row containing the point estimate and confidence interval.
#' 
#' @export
estimateByPooling <- function(simulation, useCyclops = TRUE) {
  stopifnot(class(simulation) == "simulation")
  
  maxStrataIds <- sapply(simulation, function(x) max(x$stratumId))
  startStratumId <- c(0, cumsum(maxStrataIds)[1:(length(maxStrataIds - 1))])
  addStartStratumId <- function(i) {
    simulation[[i]]$stratumId <- simulation[[i]]$stratumId + startStratumId[i]
    return(simulation[[i]])
  }
  pooled <- lapply(1:length(simulation), addStartStratumId) 
  pooled <- do.call("rbind", pooled)
  estimate <- fitIndividualModel(population = pooled, useCyclops = useCyclops)
  return(estimate)
}

#' Estimate effect by using standard meta-analysis
#' 
#' @description 
#' Estimate the effect size by first estimating the effect size at each site, and then combine using
#' standard meta-analysis techniques.
#'
#' @param simulation         An object of type \code{simulation} as created by the 
#'                           \code{\link{runSimulation}} function.
#' @param assumeRandomEffect Should meta-analysis for random effects be used?
#' @param useCyclops Should Cyclops be used? Cyclops doesn't assume a normal distribution for the 
#'                   likelihood function, and uses likelihood profiling instead.
#'
#' @return
#' A data frame with one row containing the point estimate and confidence interval.
#' 
#' @export
estimateUsingStandardMetaAnalysis <- function(simulation, 
                                              assumeRandomEffect = TRUE,
                                              useCyclops = FALSE) {
  stopifnot(class(simulation) == "simulation")

  estimates <- lapply(simulation, fitIndividualModel, useCyclops = useCyclops)
  estimates <- do.call("rbind", estimates)
  estimates <- estimates[estimates$seLogRr < 100, ]
  meta <- meta::metagen(TE = estimates$logRr, 
                        seTE = estimates$seLogRr, 
                        studlab = rep("", nrow(estimates)), 
                        byvar = NULL,
                        sm = "RR")
  s <- summary(meta)
  rnd <- s$random
  return(data.frame(rr = exp(rnd$TE),
                    ci95Lb = exp(rnd$lower),
                    ci95Ub = exp(rnd$upper),
                    logRr = rnd$TE,
                    seLogRr = rnd$seTE))
}

fitIndividualModel <- function(population, useCyclops = FALSE) {
  if (useCyclops) {
    # Using Cyclops (doesn't really assume normality)
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), 
                                              data = population, 
                                              modelType = "cox")
    fit <- Cyclops::fitCyclopsModel(cyclopsData)
    mode <- coef(fit)
    ci95 <- confint(fit, 1, level = .95)
    return(data.frame(rr = exp(mode),
                      ci95Lb = exp(ci95[2]),
                      ci95Ub = exp(ci95[3]),
                      logRr = mode,
                      seLogRr = (ci95[3] - ci95[2])/(2 * qnorm(0.975))))
  } else {
    # Using vanilla Cox regression:
    fit <- survival::coxph(Surv(time, y) ~ x + strata(stratumId), population)
    ci95 <- confint(fit)
    return(data.frame(rr = exp(fit$coefficients[1]),
                      ci95Lb = exp(ci95[1]),
                      ci95Ub = exp(ci95[2]),
                      logRr = fit$coefficients[1],
                      seLogRr = sqrt(fit$var[1,1])))
  }
}

#' Evaluate a meta-analysis approach
#'
#' @param simulationSettings   An object of type \code{simulationSettings}, created by the 
#'                             \code{\link{createSimulationSettings}} function.
#' @param iterations           Number of times the simulation should be run.
#' @param metaAnalysisFunction The function performing the meta-analysis.
#' @param seed                 The random seed to be used. Set before the first iteration.
#' @param ...                  Additional arguments that will be passed on to the meta-
#'                             analysis function.
#' 
#' @details 
#' The meta-analysis function should accept a \code{simulation} object, and produce an estimate data frame.
#'
#' @return
#' Something
#' 
#' @export
evaluateMetaAnalysisApproach <- function(simulationSettings, 
                                         iterations = 100, 
                                         metaAnalysisFunction,
                                         seed = 0,
                                         ...) {
  set.seed(seed)
  doIteration <- function(i) {
     simulation <- runSimulation(simulationSettings)
     estimate <- metaAnalysisFunction(simulation, ...)
     return(estimate)
  }
  estimates <- plyr::llply(1:iterations, doIteration, .progress = "text")
  estimates <- do.call("rbind", estimates)
  result <- data.frame(coverage = mean(estimates$ci95Lb <= simulationSettings$hazardRatio &
                                         estimates$ci95Ub >= simulationSettings$hazardRatio))
  return(result)
}

plotLikelihood <- function(simulation) {
  fitCyclopsModels <- function(population) {
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
    fit <- Cyclops::fitCyclopsModel(cyclopsData)
  }
  models <- lapply(simulation, fitCyclopsModels)
  
  # getLimits <- function(fit) {
  #   mode <- coef(fit)
  #   ci99 <- confint(fit, 1, level = .99)
  #   if (is.na(ci99[2])) {
  #     ci99[2] <- mode - (ci99[3] - mode) * 2
  #   }
  #   if (is.na(ci99[3])) {
  #     ci99[3] <- mode + (mode - ci99[2]) * 2
  #   }
  #   return(data.frame(lower = ci99[2], upper = ci99[3]))
  # }
  # limits <- lapply(models, getLimits)
  # limits <- do.call("rbind", limits)
  # x <- seq(from = min(limits$lower), to = max(limits$upper), length.out = 100)

  getCis <- function(fit, level = 0.95) {
    mode <- coef(fit)
    ci95 <- confint(fit, 1, level = level)
    if (is.na(ci95[2]) | is.na(ci95[3])) {
      return(data.frame(lower = NA, mode = NA, upper = NA))
    } else {
      return(data.frame(lower = ci95[2], mode = mode, upper = ci95[3]))
    }
  }
  limits <- lapply(models, getCis, level = 0.99)
  limits <- do.call("rbind", limits)
  x <- seq(from = min(limits$lower, na.rm = TRUE), to = max(limits$upper, na.rm = TRUE), length.out = 100)
  
  
  cis <- lapply(models, getCis)
  cis <- do.call("rbind", cis)
  cis$site <- paste("Site", 1:nrow(cis))
  
  computeY <- function(i) {
    y <- rep(0, 100)
    for (j in 1:100) {
      # Set starting coefficient, then tell Cyclops that it is fixed:
      temp <- Cyclops::fitCyclopsModel(models[[i]]$cyclopsData, startingCoefficients = x[j], fixedCoefficients = 1)
      y[j] <- temp$log_likelihood
    }
    return(data.frame(x = x, y = y - min(y), site = paste("Site", i)))
  }
  vizData <- lapply(1:length(models), computeY)
  vizData <- do.call("rbind", vizData)
  vizData$y <- vizData$y - min(vizData$y)
  vizData$site <- factor(vizData$site, levels = paste("Site", 1:length(models)))
  cis$site <- factor(cis$site, levels = paste("Site", 1:length(models)))
    
  plot <- ggplot2::ggplot(vizData, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_area(color = rgb(0, 0, 0.8), fill = rgb(0, 0, 0.8), alpha = 0.5) +
    ggplot2::geom_vline(size = 1, ggplot2::aes(xintercept = mode), data = cis) +
    ggplot2::geom_vline(size = 1, linetype = "dotted", ggplot2::aes(xintercept = lower), data = cis) +
    ggplot2::geom_vline(size = 1, linetype = "dotted", ggplot2::aes(xintercept = upper), data = cis) +
    ggplot2::xlab("Beta") +
    ggplot2::ylab("Log likelihood") +
    ggplot2::facet_grid(site~., scales = "free_y") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank())
  return(plot)
}
