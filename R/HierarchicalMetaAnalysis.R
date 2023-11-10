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


# utility function to summarize MCMC samples -- mean, HDI, std, etc.
summarizeChain <- function(chain, alpha = 0.05){
  avg = mean(chain, na.rm = TRUE)
  hdi = HDInterval::hdi(chain, credMass = 1 - alpha)
  se = sqrt(mean((chain - avg)^2))

  return(c(avg, hdi, se))
}


#' Compute a Bayesian random-effects hierarchical meta-analysis
#'
#' @description
#' Compute a Bayesian hierarchical meta-analysis (two-level model) to learn the global effect
#' with bias correction via negative control outcomes analysis.
#' Bayesian inference is performed using the Markov chain Monte Carlo (MCMC) engine BEAST.
#' Nnormal priors are used for the global effect, outcome-specific effects, and data-source-specific effects;
#' a half normal prior is used for the standard deviation; a gamma prior is used for the precision parameters.
#'
#' @param data                 A data frame containing either normal, skew-normal, custom parametric,
#'                             or grid likelihood data, with one row per database.
#' @param chainLength          Number of MCMC iterations.
#' @param burnIn               Number of MCMC iterations to consider as burn in.
#' @param subSampleFrequency   Subsample frequency for the MCMC.
#' @param alpha                The alpha (expected type I error) used for the credible intervals.
#' @param effectPriorStd       Standard deviation for the average outcome / data-source effect.
#' @param globalExposureEffectPrior  Prior mean and standard deviation for the global main effect.
#' @param effectPrecisionPrior  Shape and scale for the gamma prior of the precision term in the
#'                             random effects model (normal) for individual outcome / data-source effects.
#' @param errorPrecisionPrior  Shape and scale for the gamma prior of the precision term in the
#'                             normal model for random errors.
#' @param errorPrecisionStartValue Initial value for the error distribution's precision term.
#' @param includeSourceEffect   Whether or not to consider the data-source-specific random effects. Default is TRUE.
#' @param includeExposureEffect Whether or not to estimate the main effect of interest. Default is TRUE.
#' @param seed                 Seed for the random number generator.
#' @seealso
#' [approximateLikelihood], [computeBayesianMetaAnalysis]
#'
#' @return
#' A data frame with the point estimates, 95% credible intervals and sample standard errors
#' for the (de-biased) global main effect, the average outcome effect, the average data source effect,
#' and precision of random errors.
#' Attributes of the data frame contain the MCMC trace and the detected approximation type.
#'
#' @examples
#' data("hmaLikelihoodList")
#' estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = hmaLikelihoodList,
#' chainLength = 110000,
#' burnIn = 1e+04,
#' subSampleFrequency = 10,
#' seed = 666)
#'
#' @export
computeHierarchicalMetaAnalysis <- function(data,
                                            chainLength = 1100000,
                                            burnIn = 1e+05,
                                            subSampleFrequency = 100,
                                            alpha = 0.05,
                                            effectPriorStd = 1.0,
                                            globalExposureEffectPrior = c(0.0, 2.0),
                                            effectPrecisionPrior = c(1.0, 1.0),
                                            errorPrecisionPrior = c(1.0, 1.0),
                                            errorPrecisionStartValue = 1.0,
                                            includeSourceEffect = TRUE,
                                            includeExposureEffect = TRUE,
                                            seed = 1){
  # checks...
  if (!supportsJava8()) {
    inform("Java 8 or higher is required, but older version was found. Cannot compute estimate.")
    estimate <- data.frame(
      mean = NA,
      LB = NA,
      UB = NA,
      se = NA,
      parameter = NA
    ) # returning an empty dataframe if Java version is not okay
    traces <- NULL
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

  # make sure the data is a list of stuff
  if (typeof(data) != "list") {
    stop("data must a list of likelihood functions or data models!")
  }

  # build data models
  ## build a reference list of string labels and integer labels...
  labelReferences = buildLabelReferences(data)

  ## get the type of data in each element of the list
  type = detectApproximationType(data[[1]])

  ## construct list of dataModel objects
  dataModelList <- rJava::.jnew("java.util.ArrayList")
  for(i in 1:length(data)){
    thisDataModel = constructDataModel(data[[i]])
    #thisDataModel = rJava::.jcast(thisDataModel, "org.ohdsi.metaAnalysis.DataModel")
    dataModelList$add(rJava::.jcast(thisDataModel,
                                    "org.ohdsi.metaAnalysis.DataModel"))
  }

  # convert to Java
  ## configuration
  hmaConfiguration = rJava::.jnew("org.ohdsi.metaAnalysis.HierarchicalMetaAnalysis$HierarchicalMetaAnalysisConfiguration")

  hmaConfiguration$hierarchicalLocationHyperStdDev = as.numeric(effectPriorStd)
  hmaConfiguration$gammaHyperShape = as.numeric(effectPrecisionPrior[1])
  hmaConfiguration$gammaHyperScale = as.numeric(effectPrecisionPrior[2])
  hmaConfiguration$exposureHyperLocation = as.numeric(globalExposureEffectPrior[1])
  hmaConfiguration$exposureHyperStdDev = as.numeric(globalExposureEffectPrior[2])
  hmaConfiguration$tauShape = as.numeric(errorPrecisionPrior[1])
  hmaConfiguration$tauScale = as.numeric(errorPrecisionPrior[2])
  hmaConfiguration$startingTau = as.numeric(errorPrecisionStartValue)
  hmaConfiguration$includeSecondary = as.logical(includeSourceEffect)
  hmaConfiguration$includeExposure = as.logical(includeExposureEffect)
  hmaConfiguration$seed = rJava::.jlong(seed) # this is surplus though...

  # construct the analysis
  hierarchicalMetaAnalysis <- rJava::.jnew(
    "org.ohdsi.mcmc.Runner",
    rJava::.jcast(rJava::.jnew(
      "org.ohdsi.metaAnalysis.HierarchicalMetaAnalysis",
      rJava::.jcast(dataModelList, "java.util.List"),
      hmaConfiguration
    ), "org.ohdsi.mcmc.Analysis"),
    as.integer(chainLength),
    as.integer(burnIn),
    as.integer(subSampleFrequency),
    as.numeric(seed)
  )

  # run analysis
  hierarchicalMetaAnalysis$setConsoleWidth(getOption("width"))
  hierarchicalMetaAnalysis$run()

  # extract results and chains
  parameterNames <- hierarchicalMetaAnalysis$getParameterNames()

  # cat(sprintf("configuration: include exposure? %s \n\n",
  #             hmaConfiguration$includeExposure))

  ## get the MCMC chains one by one
  trace <- hierarchicalMetaAnalysis$getTrace(as.integer(3))
  traces <- matrix(ncol = length(parameterNames) - 2, nrow = length(trace))
  traces[, 1] <- trace
  for (i in 4:length(parameterNames)) {
    trace <- hierarchicalMetaAnalysis$getTrace(as.integer(i))
    traces[, i - 2] <- trace
  }
  ## get summary on key parameters
  parameterNames = parameterNames[-c(1:2)]

  # cat(sprintf("All parameter names: %s \n\n",
  #             paste(parameterNames, collapse = ",")))

  mainParameters = c("tau", "outcome.mean", "outcome.scale")
  if(includeSourceEffect){
    mainParameters = c(mainParameters, c("source.mean", "source.scale"))
  }
  if(includeExposureEffect){
    mainParameters = c(mainParameters, "exposure")
  }
  mainParameterIndices = which(parameterNames %in% mainParameters)

  estimates = apply(traces[,mainParameterIndices], 2, summarizeChain, alpha)
  estimates = data.frame(t(estimates), row.names = NULL)
  names(estimates) = c("mean", "LB", "UB", "se")
  estimates$parameter = mainParameters

  ## add attributes
  type = detectApproximationType(data[[1]])
  ## add parameter names to the trace matrix
  attr(estimates, "traces") <- traces
  colnames(traces) = parameterNames
  attr(estimates, "type") <- type
  attr(estimates, "ess") <- coda::effectiveSize(traces)

  return(estimates)
}



