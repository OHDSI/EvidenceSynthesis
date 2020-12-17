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

cleanData <- function(data,
                      columns,
                      minValues = rep(-100, length(columns)),
                      maxValues = rep(100, length(columns))) {
  for (i in 1:length(columns)) {
    column <- columns[i]
    if (any(is.infinite(data[, column]))) {
      warn(paste("Estimate(s) with infinite",
                 column,
                 "detected. Removing before computing meta-analysis."))
      data <- data[!is.infinite(data[, column]), ]
    }
    if (any(is.na(data[, column]))) {
      warn(paste("Estimate(s) with NA",
                 column,
                 "detected. Removing before computing meta-analysis."))
      data <- data[!is.na(data[, column]), ]
    }
    if (any(data[, column] > maxValues[i])) {
      warn(sprintf("Estimate(s) with extremely high %s (>%s) detected. Removing before computing meta-analysis.",
                   column,
                   maxValues[i]))
      data <- data[data[, column] <= maxValues[i], ]
    }
    if (any(data[, column] < minValues[i])) {
      warn(sprintf("Estimate(s) with extremely low %s (<%s) detected. Removing before computing meta-analysis.",
                   column,
                   minValues[i]))
      data <- data[data[, column] >= minValues[i], ]
    }
  }
  if (nrow(data) == 0) {
    warn("No estimates left after removing estimates with NA, infinite or extreme values")
  }
  return(data)
}

createNaEstimate <- function(type) {
  estimate <- data.frame(mu = NA,
                         mu95Lb = NA,
                         mu95Ub = NA,
                         muSe = NA,
                         tau = NA,
                         tau95Lb = NA,
                         tau95Ub = NA)
  attr(estimate, "traces") <- matrix(ncol = 2)
  attr(estimate, "type") <- type
  return(estimate)
}

isRmdCheck <- function() {
  return(Sys.getenv("_R_CHECK_PACKAGE_NAME_", "") != "")
}

isUnitTest <- function() {
  return(tolower(Sys.getenv("TESTTHAT", "")) == "true")
   
}

#' Compute a Bayesian random-effects meta-analysis
#'
#' @description
#' Compute a Bayesian meta-analysis using the Markov chain Monte Carlo (MCMC) engine BEAST.
#' A normal and half-normal prior are used for the mu and tau parameters, respectively, with standard
#' deviations as defined by the `priorSd` argument.
#'
#' @param data                 A data frame containing either normal, skew-normal, custom parametric,
#'                             or grid likelihood data, with one row per database.
#' @param chainLength          Number of MCMC iterations.
#' @param burnIn               Number of MCMC iterations to consider as burn in.
#' @param subSampleFrequency   Subsample frequency for the MCMC.
#' @param priorSd              A two-dimensional vector with the standard deviation of the prior for mu
#'                             and tau, respectively.
#' @param alpha                The alpha (expected type I error) used for the credible intervals.
#' @param seed                 The seed for the random number generator.
#'
#' @seealso
#' [approximateLikelihood], [computeFixedEffectMetaAnalysis]
#'
#' @return
#' A data frame with the point estimates and 95% credible intervals for the mu and tau parameters (the
#' mean and standard deviation of the distribution from which the per-site effect sizes are drawn).
#' Attributes of the data frame contain the MCMC trace and the detected approximation type.
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
#' estimate <- computeBayesianMetaAnalysis(approximations)
#' estimate
#'
#' # (Estimates in this example will vary due to the random simulation)
#' 
#' @export
computeBayesianMetaAnalysis <- function(data,
                                        chainLength = 1100000,
                                        burnIn = 1e+05,
                                        subSampleFrequency = 100,
                                        priorSd = c(2, 0.5),
                                        alpha = 0.05,
                                        seed = 1) {
  if (!supportsJava8()) {
    inform("Java 8 or higher is required, but older version was found. Cannot compute estimate.")
    estimate <- data.frame(mu = as.numeric(NA),
                           mu95Lb = as.numeric(NA),
                           mu95Ub = as.numeric(NA),
                           muSe = as.numeric(NA),
                           tau = as.numeric(NA),
                           tau95Lb = as.numeric(NA),
                           tau95Ub = as.numeric(NA),
                           logRr = as.numeric(NA),
                           seLogRr = as.numeric(NA),
                           row.names = NULL)
    traces <- matrix(runif(700), ncol = 7, nrow = 100)
    attr(estimate, "traces") <- traces
    attr(estimate, "type") <- "Unknown"
    attr(estimate, "ess") <- NA
    return(estimate)
  }
  if (isRmdCheck() && !isUnitTest()) {
    inform(paste("Function is executed as an example in R check:",
                 "Reducing chainLength and burnIn to reduce compute time.",
                 "Result may be unreliable"))
    chainLength <- 110000
    burnIn <- 10000
  }
  
  # Determine type based on data structure:
  if ("logRr" %in% colnames(data)) {
    inform("Detected data following normal distribution")
    type <- "normal"
    data <- cleanData(data, c("logRr", "seLogRr"), minValues = c(-100, 1e-05))
    if (nrow(data) == 0)
      return(createNaEstimate(type))
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.NormalDataModel")
    for (i in 1:nrow(data)) {
      dataModel$addLikelihoodParameters(as.numeric(c(data$logRr[i], data$seLogRr[i])),
                                        as.numeric(c(NA, NA)))
    }
    dataModel$finish()
  } else if ("gamma" %in% colnames(data)) {
    inform("Detected data following custom parameric distribution")
    type <- "custom"
    data <- cleanData(data, c("mu", "sigma", "gamma"), minValues = c(-100, 1e-05, -100))
    if (nrow(data) == 0)
      return(createNaEstimate(type))
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.ParametricDataModel")
    for (i in 1:nrow(data)) {
      dataModel$addLikelihoodParameters(as.numeric(c(data$mu[i], data$sigma[i], data$gamma[i])),
                                        as.numeric(c(NA, NA)))
    }
    dataModel$finish()
  } else if ("alpha" %in% colnames(data)) {
    inform("Detected data following skew normal distribution")
    type <- "skew normal"
    data <- cleanData(data,
                      c("mu", "sigma", "alpha"),
                      minValues = c(-100, 1e-05, -10000),
                      maxValues = c(100, 10000, 10000))
    if (nrow(data) == 0)
      return(createNaEstimate(type))
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.SkewNormalDataModel")
    for (i in 1:nrow(data)) {
      dataModel$addLikelihoodParameters(as.numeric(c(data$mu[i], data$sigma[i], data$alpha[i])),
                                        as.numeric(c(NA, NA)))
    }
    dataModel$finish()
  } else if (is.list(data) && !is.data.frame(data)) {
    inform("Detected (pooled) patient-level data")
    type <- "pooled"
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.CoxDataModel")
    for (i in 1:length(data)) {
      dataModel$addLikelihoodData(as.integer(data[[i]]$stratumId),
                                  as.integer(data[[i]]$y),
                                  as.numeric(data[[i]]$time),
                                  as.numeric(data[[i]]$x))
    }
    dataModel$finish()
  } else {
    inform("Detected data following grid distribution")
    type <- "grid"
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.ExtendingEmpiricalDataModel")
    x <- as.numeric(colnames(data))
    if (any(is.na(x))) {
      stop("Expecting grid data, but not all column names are numeric")
    }
    data <- as.matrix(data)
    for (i in 1:nrow(data)) {
      dataModel$addLikelihoodParameters(x, data[i, ])
    }
    dataModel$finish()
  }
  
  inform("Performing MCMC. This may take a while")
  
  prior <- rJava::.jnew("org.ohdsi.metaAnalysis.HalfNormalOnStdDevPrior", 0, as.numeric(priorSd[2]))
  
  metaAnalysis <- rJava::.jnew("org.ohdsi.metaAnalysis.Runner",
                               rJava::.jcast(rJava::.jnew("org.ohdsi.metaAnalysis.MetaAnalysis",
                                                          rJava::.jcast(dataModel, "org.ohdsi.metaAnalysis.DataModel"),
                                                          rJava::.jcast(prior, "org.ohdsi.metaAnalysis.ScalePrior"),
                                                          as.numeric(priorSd[1])), "org.ohdsi.metaAnalysis.Analysis"),
                               as.integer(chainLength),
                               as.integer(burnIn),
                               as.integer(subSampleFrequency),
                               as.numeric(seed))
  metaAnalysis$setConsoleWidth(getOption("width"))
  metaAnalysis$run()
  parameterNames <- metaAnalysis$getParameterNames()
  trace <- metaAnalysis$getTrace(as.integer(3))
  traces <- matrix(ncol = length(parameterNames) - 2, nrow = length(trace))
  traces[, 1] <- trace
  for (i in 4:length(parameterNames)) {
    trace <- metaAnalysis$getTrace(as.integer(i))
    traces[, i - 2] <- trace
  }
  hdiMu <- HDInterval::hdi(traces[, 1], credMass = 1 - alpha)
  hdiTau <- HDInterval::hdi(traces[, 2], credMass = 1 - alpha)
  mu <- mean(traces[, 1])
  estimate <- data.frame(mu = mu,
                         mu95Lb = hdiMu[1],
                         mu95Ub = hdiMu[2],
                         muSe = sqrt(mean((traces[, 1] - mu)^2)),
                         tau = median(traces[, 2]),
                         tau95Lb = hdiTau[1],
                         tau95Ub = hdiTau[2],
                         logRr = mu,
                         seLogRr = sqrt(mean((traces[, 1] - mu)^2)),
                         row.names = NULL)
  attr(estimate, "traces") <- traces
  attr(estimate, "type") <- type
  attr(estimate, "ess") <- coda::effectiveSize(traces)
  return(estimate)
}
