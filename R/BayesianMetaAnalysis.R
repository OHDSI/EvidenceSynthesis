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
createNaEstimate <- function(type) {
  estimate <- data.frame(
    mu = NA,
    mu95Lb = NA,
    mu95Ub = NA,
    muSe = NA,
    tau = NA,
    tau95Lb = NA,
    tau95Ub = NA
  )
  attr(estimate, "traces") <- matrix(ncol = 2)
  attr(estimate, "type") <- type
  return(estimate)
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
#' @param robust               Whether or not to use a t-distribution model; default: FALSE.
#' @param df                   Degrees of freedom for the t-model, only used if robust is TRUE.
#' @param seed                 The seed for the random number generator.
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
                                        robust = FALSE,
                                        df = 4,
                                        seed = 1) {
  if (!supportsJava8()) {
    inform("Java 8 or higher is required, but older version was found. Cannot compute estimate.")
    estimate <- data.frame(
      mu = as.numeric(NA),
      mu95Lb = as.numeric(NA),
      mu95Ub = as.numeric(NA),
      muSe = as.numeric(NA),
      tau = as.numeric(NA),
      tau95Lb = as.numeric(NA),
      tau95Ub = as.numeric(NA),
      logRr = as.numeric(NA),
      seLogRr = as.numeric(NA),
      row.names = NULL
    )
    traces <- matrix(runif(700), ncol = 7, nrow = 100)
    attr(estimate, "traces") <- traces
    attr(estimate, "type") <- "Unknown"
    attr(estimate, "ess") <- NA
    return(estimate)
  }
  if (isRmdCheck() && !isUnitTest()) {
    inform(paste(
      "Function is executed as an example in R check:",
      "Reducing chainLength and burnIn to reduce compute time.",
      "Result may be unreliable"
    ))
    chainLength <- 110000
    burnIn <- 10000
    Sys.sleep(1) # To avoid CRAN message about CPU time being more than 2.5. times elapsed time
  }

  # refactored: using utils function to create a `dataModel` object
  dataModel = constructDataModel(data)


  inform("Performing MCMC. This may take a while")
  prior <- rJava::.jnew("org.ohdsi.metaAnalysis.HalfNormalOnStdDevPrior", 0, as.numeric(priorSd[2]))
  if (robust) {
    metaAnalysis <- rJava::.jnew(
      "org.ohdsi.mcmc.Runner",
      rJava::.jcast(
        rJava::.jnew(
          "org.ohdsi.metaAnalysis.RobustMetaAnalysis",
          rJava::.jcast(dataModel, "org.ohdsi.metaAnalysis.DataModel"),
          rJava::.jcast(prior, "org.ohdsi.metaAnalysis.ScalePrior"),
          as.numeric(priorSd[1]),
          as.numeric(df)
        ),
        "org.ohdsi.mcmc.Analysis"
      ),
      as.integer(chainLength),
      as.integer(burnIn),
      as.integer(subSampleFrequency),
      as.numeric(seed)
    )
  } else {
    metaAnalysis <- rJava::.jnew(
      "org.ohdsi.mcmc.Runner",
      rJava::.jcast(rJava::.jnew(
        "org.ohdsi.metaAnalysis.MetaAnalysis",
        rJava::.jcast(dataModel, "org.ohdsi.metaAnalysis.DataModel"),
        rJava::.jcast(prior, "org.ohdsi.metaAnalysis.ScalePrior"),
        as.numeric(priorSd[1])
      ), "org.ohdsi.mcmc.Analysis"),
      as.integer(chainLength),
      as.integer(burnIn),
      as.integer(subSampleFrequency),
      as.numeric(seed)
    )
  }

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
  estimate <- data.frame(
    mu = mu,
    mu95Lb = hdiMu[1],
    mu95Ub = hdiMu[2],
    muSe = sqrt(mean((traces[, 1] - mu)^2)),
    tau = median(traces[, 2]),
    tau95Lb = hdiTau[1],
    tau95Ub = hdiTau[2],
    logRr = mu,
    seLogRr = sqrt(mean((traces[, 1] - mu)^2)),
    row.names = NULL
  )
  attr(estimate, "traces") <- traces
  attr(estimate, "type") <- type
  attr(estimate, "ess") <- coda::effectiveSize(traces)
  return(estimate)
}
