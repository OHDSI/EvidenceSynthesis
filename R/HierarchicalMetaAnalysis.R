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


#' Utility function to summarize MCMC samples (posterior mean, median, HDI, std, etc.)
#'
#' @param chain A vector of posterior samples from MCMC.
#' @param alpha Alpha level for the credible interval.
#'
#' @export
summarizeChain <- function(chain, alpha = 0.05){
  avg = mean(chain, na.rm = TRUE)
  med = median(chain, na.rm = TRUE)
  hdi = HDInterval::hdi(chain, credMass = 1 - alpha)
  se = sqrt(mean((chain - avg)^2))

  res = c(avg, med, hdi, se)
  names(res) = c("mean", "median", "LB", "UB", "se")

  return(res)
}

#' Generate settings for the Bayesian random-effects hierarchical meta-analysis model
#'
#' @description
#' This function generates a settings list for fitting a Bayesian hierarchical meta-analysis model.
#' See `computeHierarchicalMetaAnalysis()` for more details.
#'
#' @param primaryEffectPriorStd   Standard deviation for the average outcome effect.
#' @param secondaryEffectPriorStd Standard deviation for the average data-source effect.
#' @param globalExposureEffectPriorMean  Prior mean for the global main exposure effect;
#'                                       can be a multiple entry vector if there are multiple outcomes of interest
#' @param globalExposureEffectPriorStd   Prior standard deviation for the global main exposure effect;
#'                                       can be a multiple entry vector if there are multiple outcomes of interest
#' @param primaryEffectPrecisionPrior  Shape and scale for the gamma prior of the precision term in the
#'                             random effects model (normal) for individual outcome effects.
#' @param secondaryEffectPrecisionPrior Shape and scale for the gamma prior of the precision term in the
#'                             random effects model (normal) for individual data-source effects.
#' @param errorPrecisionPrior  Shape and scale for the gamma prior of the precision term in the
#'                             normal model for random errors.
#' @param errorPrecisionStartValue Initial value for the error distribution's precision term.
#' @param includeSourceEffect   Whether or not to consider the data-source-specific (secondary) random effects. Default is TRUE.
#' @param includeExposureEffect Whether or not to estimate the main effect of interest. Default is TRUE.
#' @param exposureEffectCount   Number of main outcomes of interest to estimate effect for? Default = 1
#' @param separateExposurePrior Use a separable prior on the main exposure effect? Default is FALSE.
#' @param chainLength          Number of MCMC iterations.
#' @param burnIn               Number of MCMC iterations to consider as burn in.
#' @param subSampleFrequency   Subsample ("thinning") frequency for the MCMC.
#'
#' @return A list with all the settings to use in the `computeHierarchicalMetaAnalysis()` function.
#'
#' @seealso [computeHierarchicalMetaAnalysis]
#'
#'
#' @export
generateBayesianHMAsettings <- function(primaryEffectPriorStd = 1.0,
                                        secondaryEffectPriorStd = 1.0,
                                        globalExposureEffectPriorMean = c(0.0),
                                        globalExposureEffectPriorStd = c(2.0),
                                        primaryEffectPrecisionPrior = c(1.0, 1.0),
                                        secondaryEffectPrecisionPrior = c(1.0, 1.0),
                                        errorPrecisionPrior = c(1.0, 1.0),
                                        errorPrecisionStartValue = 1.0,
                                        includeSourceEffect = TRUE,
                                        includeExposureEffect = TRUE,
                                        exposureEffectCount = 1,
                                        separateExposurePrior = FALSE,
                                        chainLength = 1100000,
                                        burnIn = 1e+05,
                                        subSampleFrequency = 100){

  settings = list(
    primaryEffectPriorStd = primaryEffectPriorStd,
    secondaryEffectPriorStd = secondaryEffectPriorStd,
    primaryEffectPrecisionPrior = primaryEffectPrecisionPrior,
    secondaryEffectPrecisionPrior = secondaryEffectPrecisionPrior,
    errorPrecisionPrior = errorPrecisionPrior,
    errorPrecisionStartValue = errorPrecisionStartValue,
    includeSourceEffect = includeSourceEffect,
    includeExposureEffect = includeExposureEffect,
    exposureEffectCount = exposureEffectCount,
    separateExposurePrior = separateExposurePrior,
    chainLength = chainLength,
    burnIn = burnIn,
    subSampleFrequency = subSampleFrequency
  )

  dimExposurePriorMean = length(globalExposureEffectPriorMean)
  dimExposurePriorStd = length(globalExposureEffectPriorStd)

  if(exposureEffectCount == 1){
    if(dimExposurePriorStd != 1 | dimExposurePriorMean != 1){
      stop("Dimensions of `globalExposureEffectPriorMean` and `globalExposureEffectPriorStd` should either be 1 or match `exposureEffectCount`!")
    }
    settings$globalExposureEffectPriorMean = globalExposureEffectPriorMean
    settings$globalExposureEffectPriorStd = globalExposureEffectPriorStd
  }else if(exposureEffectCount < 1){
    settings$globalExposureEffectPriorMean = globalExposureEffectPriorMean
    settings$globalExposureEffectPriorStd = globalExposureEffectPriorStd
  }else{
    if(dimExposurePriorMean == 1){
      settings$globalExposureEffectPriorMean = rep(globalExposureEffectPriorMean, exposureEffectCount)
    }else{
      if(dimExposurePriorMean != exposureEffectCount){
        stop("Dimension of `globalExposureEffectPriorMean` should either be 1 or match `exposureEffectCount`!")
      }
      settings$globalExposureEffectPriorMean = globalExposureEffectPriorMean
    }

    if(dimExposurePriorStd == 1){
      settings$globalExposureEffectPriorStd = rep(globalExposureEffectPriorStd, exposureEffectCount)
    }else{
      if(dimExposurePriorStd != exposureEffectCount){
        stop("Dimension of `globalExposureEffectPriorStd` should either be 1 or match `exposureEffectCount`!")
      }
      settings$globalExposureEffectPriorStd  = globalExposureEffectPriorStd
    }
  }

  return(settings)
}



#' Compute a Bayesian random-effects hierarchical meta-analysis
#'
#' @description
#' Compute a Bayesian hierarchical meta-analysis (two-level model) to learn the global effect
#' with bias correction via negative control outcomes analysis.
#' Bayesian inference is performed using the Markov chain Monte Carlo (MCMC) engine BEAST.
#' Normal priors are used for the global effect, outcome-specific effects, and data-source-specific effects;
#' a half normal prior is used for the standard deviation; a gamma prior is used for the precision parameters.
#'
#' @param data                 A data frame containing either normal, skew-normal, custom parametric,
#'                             or grid likelihood data, with one row per database.
#' @param settings             Model settings list generated by `generateBayesianHMAsettings()`
#' @param alpha                The alpha (expected type I error) used for the credible intervals.
#' @param seed                 Seed for the random number generator.
#' @param showProgressBar      Showing a progress bar for MCMC?
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
#' seed = 666)
#'
#' @export
computeHierarchicalMetaAnalysis <- function(data,
                                            settings = generateBayesianHMAsettings(),
                                            alpha = 0.05,
                                            seed = 1,
                                            showProgressBar = TRUE){
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

  # read settings list
  chainLength = settings$chainLength
  burnIn = settings$burnIn
  subSampleFrequency = settings$subSampleFrequency
  exposureEffectCount = settings$exposureEffectCount

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

  cat("Data model list built!\n\n")

  # convert to Java
  ## configuration
  hmaConfiguration = rJava::.jnew("org.ohdsi.metaAnalysis.HierarchicalMetaAnalysis$HierarchicalMetaAnalysisConfiguration")

  hmaConfiguration$hierarchicalLocationPrimaryHyperStdDev = as.numeric(settings$primaryEffectPriorStd)
  hmaConfiguration$hierarchicalLocationSecondaryHyperStdDev = as.numeric(settings$secondaryEffectPriorStd)
  hmaConfiguration$gammaHyperPrimaryShape = as.numeric(settings$primaryEffectPrecisionPrior[1])
  hmaConfiguration$gammaHyperPrimaryScale = as.numeric(settings$primaryEffectPrecisionPrior[2])
  hmaConfiguration$gammaHyperSecondaryShape = as.numeric(settings$secondaryEffectPrecisionPrior[1])
  hmaConfiguration$gammaHyperSecondaryScale = as.numeric(settings$secondaryEffectPrecisionPrior[2])
  #hmaConfiguration$exposureHyperStdDev = as.numeric(settings$globalExposureEffectPriorStd[1])
  hmaConfiguration$tauShape = as.numeric(settings$errorPrecisionPrior[1])
  hmaConfiguration$tauScale = as.numeric(settings$errorPrecisionPrior[2])
  hmaConfiguration$startingTau = as.numeric(settings$errorPrecisionStartValue)
  hmaConfiguration$includeSecondary = as.logical(settings$includeSourceEffect)
  hmaConfiguration$includeExposure = as.logical(settings$includeExposureEffect)
  hmaConfiguration$separateEffectPrior = as.logical(settings$separateExposurePrior)
  hmaConfiguration$effectCount = as.integer(exposureEffectCount)
  hmaConfiguration$seed = rJava::.jlong(seed)

  ## deal with exposure prior mean: pop out the 0.0 entry first
  hmaConfiguration$exposureHyperLocation$remove(rJava::.jnew("java/lang/Double", 0.0))
  for(i in c(1:exposureEffectCount)){
    hmaConfiguration$exposureHyperLocation$add(rJava::.jnew("java.lang.Double",
                                                            as.numeric(settings$globalExposureEffectPriorMean[i])))
  }
  ## similarly, work on exposure prior std (pop out the default 2.0 entry)
  hmaConfiguration$exposureHyperStdDev$remove(rJava::.jnew("java/lang/Double", 2.0))
  for(i in c(1:exposureEffectCount)){
    hmaConfiguration$exposureHyperStdDev$add(rJava::.jnew("java.lang.Double",
                                                            as.numeric(settings$globalExposureEffectPriorStd[i])))
  }

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
    as.numeric(seed),
    as.logical(showProgressBar)
  )

  # run analysis
  hierarchicalMetaAnalysis$setConsoleWidth(getOption("width"))
  inform("Performing MCMC. This may take a while")
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
  ## add parameter names to the trace matrix
  colnames(traces) = parameterNames

  # cat(sprintf("All parameter names: %s \n\n",
  #             paste(parameterNames, collapse = ",")))

  mainParameters = c("tau", "outcome.mean", "outcome.scale")
  if(settings$includeSourceEffect){
    mainParameters = c(mainParameters, c("source.mean", "source.scale"))
  }
  if(settings$includeExposureEffect){
    # if(exposureEffectCount == 1){
    #   exposureNames = c("exposure")
    # }else if(exposureEffectCount > 1){
    if(exposureEffectCount >= 1){
      exposureNames = paste0("exposure", c(1:exposureEffectCount))
    }else{
      exposureNames = NULL
    }
    mainParameters = c(mainParameters, exposureNames)
    ## TODO: need to modify for multiple exposure effects case!!!
    if(settings$separateExposurePrior){
      if(exposureEffectCount >= 1){
        effectBiasCol = paste0("outcome",
                               c((length(data)-exposureEffectCount+1):length(data))
        )
        effectBiasSamps = traces[,effectBiasCol]
        effectSamps = traces[,exposureNames]
        traces = cbind(traces, effectSamps)
        if(exposureEffectCount == 1){
          unadjEffectNames = c("unadjustedExposure")
        }else{
          unadjEffectNames = paste0("unadjustedExposure", c(1:exposureEffectCount))
        }
        newN = ncol(traces)
        newCols = c((newN-exposureEffectCount+1):newN)
        colnames(traces)[newCols] = unadjEffectNames
        if(settings$includeSourceEffect){
          if(exposureEffectCount > 1){
            for(i in 1:exposureEffectCount){
              effectBiasSamps[,i] = effectBiasSamps[,i] + traces[,"source.mean"]
            }
          }else{
            effectBiasSamps = effectBiasSamps + traces[,"source.mean"]
          }
        }
        traces[,exposureNames] = effectSamps - effectBiasSamps
      }
    }
  }
  mainParameterIndices = which(parameterNames %in% mainParameters)

  estimates = apply(traces[,mainParameterIndices], 2, summarizeChain, alpha)
  estimates = data.frame(t(estimates), row.names = NULL)
  names(estimates) = c("mean", "median", "LB", "UB", "se")

  # print(mainParameters)
  # print(estimates)

  estimates$parameter = mainParameters

  ## add attributes
  type = detectApproximationType(data[[1]])
  attr(estimates, "traces") <- traces
  attr(estimates, "type") <- type
  attr(estimates, "ess") <- coda::effectiveSize(traces)
  attr(estimates, "sourceLabels") <- labelReferences
  attr(estimates, "settings") <- settings

  return(estimates)
}


#' Compute source-specific biases and bias-corrected estimates from hierarchical meta analysis results
#'
#' @description Extract source-specific biases and obtain bias-corrected estimates for each data source,
#' given the results from `computeHierarchicalMetaAnalysis()`.
#'
#' @return A data frame with point estimates, 95% credible intervals and sample standard errors for the
#' effect size after bias correction within each data source.
#'
#' @param estimates A data frame as output from the `computeHierarchicalMetaAnalysis()` function.
#' @param alpha     The alpha (expected type I error) used for the credible intervals.
#'
#'
#' @seealso [computeHierarchicalMetaAnalysis]
#'
#' @export
extractSourceSpecificEffects <- function(estimates, alpha = 0.05){
  if(!"source.mean" %in% estimates$parameter){
    stop("Cannot find source.mean as a main parameter in the fitted model!")
  }
  if(!"traces" %in% names(attributes(estimates))){
    stop("Cannot find MCMC samples from the input model object!")
  }
  traces = attr(estimates, "traces")
  #sourceColumns = which(stringr::str_detect(colnames(traces), "source[0-9]+"))
  sourceColumns = which(grepl("source[0-9]+", colnames(traces)))
  if(length(sourceColumns) < 1){
    stop("Cannot find posterior samples for source-specific effects!")
  }
  sourceColumnNames = colnames(traces)[sourceColumns]

  #effectColumns = which(stringr::str_starts(estimates$parameter, "exposure*"))
  effectColumns = which(grepl("^exposure*",estimates$parameter))
  if(length(effectColumns) < 1){
    stop("Cannot find main exposure effects in main parameters!")
  }
  effectColumnNames = estimates$parameter[effectColumns]

  newTraces = NULL #matrix(nrow = nrow(traces), ncol = length(effectColumns)*length(sourceColumns))
  newColumns = NULL
  for(e in effectColumnNames){
    for(s in sourceColumnNames){
      this.source.effect = traces[,e] + (traces[,s] - traces[,"source.mean"])
      newTraces = cbind(newTraces, this.source.effect)
      newColumns = c(newColumns, sprintf("%s.%s", e, s))
    }
  }

  colnames(newTraces) = newColumns

  newEstimates = apply(newTraces, 2, summarizeChain, alpha)
  newEstimates = data.frame(t(newEstimates), row.names = NULL)
  names(newEstimates) = c("mean", "median", "LB", "UB", "se")
  newEstimates$parameter = newColumns

  attr(newEstimates, "traces") = newTraces
  attr(newEstimates, "ess") <- coda::effectiveSize(newTraces)
  if("sourceLabels" %in% names(attributes(estimates))){
    attr(newEstimates, "sourceLabels") <- attr(estimates, "sourceLabels")
  }else{
    attr(newEstimates, "sourceLabels") <- NULL
  }

  return(newEstimates)

}
