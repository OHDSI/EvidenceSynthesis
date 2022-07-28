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

#' Approximate simple Bayesian posterior
#'
#' @description
#' Approximate a Bayesian posterior from a `Cyclops` likelihood profile and normal prior
#' using the Markov chain Monte Carlo engine BEAST.
#'
#' @param likelihoodProfile    Named vector containing grid likelihood data from `Cyclops`.
#' @param chainLength          Number of MCMC iterations.
#' @param burnIn               Number of MCMC iterations to consider as burn in.
#' @param subSampleFrequency   Subsample frequency for the MCMC.
#' @param priorMean            Prior mean for the regression parameter
#' @param priorSd              Prior standard deviation for the regression parameter
#' @param startingValue        Initial state for regression parameter
#' @param seed                 Seed for the random number generator.
#'
#'
#' @return
#' A data frame with the point estimates and 95% credible intervals for the regression parameter.
#' Attributes of the data frame contain the MCMC trace for diagnostics.
#'
#' @examples
#' # Simulate some data for this example:
#' population <- simulatePopulations(createSimulationSettings(nSites = 1))[[1]]
#'
#' # Fit a Cox regression at each data site, and approximate likelihood function:
#' cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'                                           data = population,
#'                                           modelType = "cox")
#' cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#' likelihoodProfile <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "grid")
#'
#' # Run MCMC
#' mcmcTraces <- approximateSimplePosterior(likelihoodProfile = likelihoodProfile,
#'                                          priorMean = 0, priorSd = 100)
#'
#' # Report posterior expectation
#' mean(mcmcTraces$theta)
#'
#' # (Estimates in this example will vary due to the random simulation)
#'
#' @export
approximateSimplePosterior <- function(likelihoodProfile,
                                       chainLength = 1100000,
                                       burnIn = 1e+05,
                                       subSampleFrequency = 100,
                                       priorMean = 0,
                                       priorSd = 0.5,
                                       startingValue = 0.0,
                                       seed = 1) {
  if (!supportsJava8()) {
    inform("Java 8 or higher is required, but older version was found. Cannot compute estimate.")
    return(NULL) # TODO Return something more informative
  }
  if (isRmdCheck() && !isUnitTest()) {
    inform(paste("Function is executed as an example in R check:",
                 "Reducing chainLength and burnIn to reduce compute time.",
                 "Result may be unreliable"))
    chainLength <- 110000
    burnIn <- 10000
  }

  inform("Detected data following grid distribution")
  type <- "grid"
  dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.ExtendingEmpiricalDataModel")
  x <- as.numeric(names(likelihoodProfile))
  if (any(is.na(x))) {
    stop("Expecting grid data, but not all column names are numeric")
  }
  dataModel$addLikelihoodParameters(x, likelihoodProfile)
  dataModel$finish()

  inform("Performing MCMC. This may take a while")
  analysis <- rJava::.jnew("org.ohdsi.simpleDesign.ProfileNormalAnalysis",
                           rJava::.jcast(dataModel, "org.ohdsi.metaAnalysis.EmpiricalDataModel"),
                           priorMean, priorSd, startingValue)
  runner <- rJava::.jnew("org.ohdsi.mcmc/Runner",
                         rJava::.jcast(analysis, "org.ohdsi.mcmc.Analysis"),
                         as.integer(chainLength),
                         as.integer(burnIn),
                         as.integer(subSampleFrequency),
                         as.numeric(seed))

  runner$setConsoleWidth(getOption("width"))
  runner$run()

  runner$processSamples() # TODO Comment out at some point

  names <- runner$getParameterNames()
  traces <- as.data.frame(sapply(1:length(names), function(i) { runner$getTrace(as.integer(i)) }))
  colnames(traces) <- names

  return(traces)
}





#' Approximate Bayesian posterior for hierarchical Normal model
#'
#' @description
#' Approximate a Bayesian posterior from a set of`Cyclops` likelihood profiles
#' under a hierarchical normal model
#' using the Markov chain Monte Carlo engine BEAST.
#'
#' @param likelihoodProfiles     List of grid likelihoods profiled with `Cyclops`.
#' @param chainLength            Number of MCMC iterations.
#' @param burnIn                 Number of MCMC iterations to consider as burn in.
#' @param subSampleFrequency     Subsample frequency for the MCMC.
#' @param effectPriorMean        Prior mean for global parameter
#' @param effectPriorSd          Prior standard deviation for the global parameter
#' @param nu0                    Prior "sample size" for precision
#'                               (with precision ~ gamma(nu0/2, nu0*sigma0/2))
#' @param sigma0                 Prior "variance" for precision
#'                               (with precision ~ gamma(nu0/2, nu0*sigma0/2))
#' @param effectStartingValue    Initial value for global & local parameter
#' @param precisionStartingValue Initial value for the precision
#' @param seed                   Seed for the random number generator.
#'
#'
#' @return
#' A data frame with the point estimates and 95% credible intervals for the
#' the global and local parameter, as well as the global precision.
#' Attributes of the data frame contain the MCMC trace for diagnostics.
#'
#' @examples
#' #TBD
#'
#' @export
approximateHierarchicalNormalPosterior <- function(likelihoodProfiles,
                                                   chainLength = 1100000,
                                                   burnIn = 1e+05,
                                                   subSampleFrequency = 100,
                                                   effectPriorMean = 0,
                                                   effectPriorSd = 0.5,
                                                   nu0 = 1,
                                                   sigma0 = 1,
                                                   effectStartingValue = 0.0,
                                                   precisionStartingValue = 1,
                                                   seed = 1) {
  if (!supportsJava8()) {
    inform("Java 8 or higher is required, but older version was found. Cannot compute estimate.")
    return(NULL) # TODO Return something more informative
  }
  if (isRmdCheck() && !isUnitTest()) {
    inform(paste("Function is executed as an example in R check:",
                 "Reducing chainLength and burnIn to reduce compute time.",
                 "Result may be unreliable"))
    chainLength <- 110000
    burnIn <- 10000
  }

  # specify data model for each likelihood profile in the list
  # and compile them into an R list
  inform("Detected data following grid distribution")
  type <- "grid"

  dataModelList <- list()
  for(likelihood in likelihoodProfiles){
    this.dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.ExtendingEmpiricalDataModel")
    x = as.numeric(names(likelihood))
    if (any(is.na(x))) {
      stop("Expecting grid data, but not all column names are numeric")
    }
    this.dataModel$addLikelihoodParameters(x, likelihood)
    this.dataModel$finish()

    # cast data model to class "org.ohdsi.metaAnalysis.EmpiricalDataModel" here instead of later on
    this.dataModel = rJava::.jcast(this.dataModel, "org.ohdsi.metaAnalysis.EmpiricalDataModel")

    dataModelList = c(dataModelList, this.dataModel)
  }

  # translate nu0 and sigma0 to shape and rate parameters for gamma prior
  precisionShape = nu0/2
  precisionRate = nu0 * sigma0/2 # assuming we use rate (not scale) for gamma distribution in BEAST

  # specify model and run the MCMC
  inform("Performing MCMC. This may take a while")
  analysis <- rJava::.jnew("org.ohdsi.simpleDesign.ProfileHierarchicalNormalAnalysis", #<- use Hierarchical here
                           rJava::.jarray(dataModelList),  #<- not sure this is correct...
                           effectPriorMean, effectPriorSd, effectStartingValue, #<- prior and init info for beta
                           precisionShape, precisionRate, precisionStartingValue) #<- prior and init info for precision
  runner <- rJava::.jnew("org.ohdsi.mcmc/Runner",
                         rJava::.jcast(analysis, "org.ohdsi.mcmc.Analysis"),
                         as.integer(chainLength),
                         as.integer(burnIn),
                         as.integer(subSampleFrequency),
                         as.numeric(seed))

  runner$setConsoleWidth(getOption("width"))
  runner$run()

  runner$processSamples() # TODO Comment out at some point

  names <- runner$getParameterNames()
  traces <- as.data.frame(sapply(1:length(names), function(i) { runner$getTrace(as.integer(i)) }))
  colnames(traces) <- names

  return(traces)
}
