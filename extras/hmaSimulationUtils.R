# @author Fan Bu
# simulation utilities for HMA

# rJava::.jinit(parameters="-Xmx32g", force.init = TRUE)
# options(java.parameters = c("-Xms200g", "-Xmx200g"))

library(EvidenceSynthesis)
library(dplyr)

#### settings of the simulation experiments -----
#setwd("./extras")
cachepath = "cache-3"
if(!dir.exists(cachepath)){dir.create(cachepath)}

# numrep = 50 # number of repetitions of the experiments

#trueRRs = c(1, 2, 4, 6)
trueRRs = c(1, 1.5, 2, 4)

nSites = 10
sitePop = 10000
mNegativeControls = 50
meanBias = 0.5
biasStd = 0.1
meanSiteEffect = 0 # has to be fixed at 0 !
siteEffectStd = 0.15

# seed0 = 71
# set.seed(seed0)

## simulation function
simulateAndApproximate <- function(trueRR, seed = 42, cachePath = NULL, repId = 1){
  metaPopulations = simulateMetaAnalysisWithNegativeControls(meanExposureEffect = log(trueRR),
                                                             nSites = nSites,
                                                             mNegativeControls = mNegativeControls,
                                                             sitePop = sitePop,
                                                             treatedFraction = 0.3,
                                                             minBackgroundHazard = 0.05,
                                                             maxBackgroundHazard = 0.05,
                                                             nStrata = 1,
                                                             meanBias = meanBias,
                                                             biasStd = biasStd,
                                                             meanSiteEffect = meanSiteEffect,
                                                             siteEffectStd = siteEffectStd,
                                                             seed = seed)
  metaLPs = lapply(metaPopulations, createApproximations, "adaptive grid")

  if(!is.null(cachePath)){
    saveRDS(metaLPs, file.path(cachePath, sprintf("metaLPs-%s-%s.rds", trueRR, repId)))
  }

  return(metaLPs)
}

## function to shift likelihood profiles (only do it on the exposure profiles (last entry of meta list))
shiftLP <- function(LPlist, oldRR = 1, newRR = 2, cachePath = NULL, repId = 1){
  shift = log(newRR) - log(oldRR)
  exposureLP = LPlist[[length(LPlist)]]
  for(l in 1:length(exposureLP)){
    exposureLP[[l]]$point = exposureLP[[l]]$point + shift
  }
  LPlist[[length(LPlist)]] = exposureLP

  if(!is.null(cachePath)){
    saveRDS(LPlist, file.path(cachePath, sprintf("metaLPs-%s-%s.rds", newRR, repId)))
  }

  return(LPlist)
}


# summarize estimates in string
estimateAsString <- function(v, digits = 3, exponentiate = TRUE){
  if(exponentiate){v = exp(v)}
  v = round(v, digits = digits)
  s = sprintf("%s (%s, %s)", v[1], v[2], v[3])
  return(s)
}

# summarize estimates in (postMedian, lb95, ub95) tuples, with a method label column
estimatesAsVector <- function(v, label = "standard", exponentiate = TRUE){
  if(exponentiate){v = exp(v)}

  res = as.list(v)
  names(res) = c("Estimate", "LB", "UB")
  res$label = label
  return(res)
}


# fit four different models
fitModels <- function(trueRR, metaLPs, seed = 666, exportAsString = TRUE,
                      cachePath = "cache2", repId = 1,
                      ...){

  if(exportAsString){
    res = list(RR = trueRR)
  }else{
    res = NULL
  }

  print("fit models one by one! ")


  # 1 ## fit a joint model with main exposure effect included (uninformative mean):
  resFile = file.path(cachePath, sprintf("fittedModel-uninformed-%s-%s.rds", trueRR, repId))
  if(file.exists(resFile)){
    maWithExposure = readRDS(resFile)
    #cat("loading from fitted model object...\n")
  }else{
    settings = generateBayesianHMAsettings(globalExposureEffectPriorMean = c(0),
                                           globalExposureEffectPriorStd = c(10.0),
                                           exposureEffectCount = 1,
                                           ...)
    maWithExposure = computeHierarchicalMetaAnalysis(data = metaLPs,
                                                     settings = settings,
                                                     seed = seed)
    saveRDS(maWithExposure, file = resFile)
  }

  # final row is main exposure...
  estRow = nrow(maWithExposure)

  if(exportAsString){
    res$uninformed = estimateAsString(maWithExposure[estRow,2:4])
  }else{
    res = bind_rows(res, as.data.frame(estimatesAsVector(maWithExposure[estRow,2:4], label = "uninformed")))
  }
  rm(maWithExposure)

  # 2. try using a "friendlier" prior for the meta-analytic effect
  resFile = file.path(cachePath, sprintf("fittedModel-informed-%s-%s.rds", trueRR, repId))
  if(file.exists(resFile)){
    maWithExposure = readRDS(resFile)
  }else{
    settings = generateBayesianHMAsettings(globalExposureEffectPriorMean = c(log(trueRR)),
                                           globalExposureEffectPriorStd = c(2.0),
                                           exposureEffectCount = 1,
                                           ...)
    maWithExposure = computeHierarchicalMetaAnalysis(data = metaLPs,
                                                     settings = settings,
                                                     seed = seed)
    saveRDS(maWithExposure, file = resFile)
  }

  estRow = nrow(maWithExposure)
  if(exportAsString){
    res$informed = estimateAsString(maWithExposure[estRow,2:4])
  }else{
    res = bind_rows(res, as.data.frame(estimatesAsVector(maWithExposure[estRow,2:4], label = "informed")))
  }
  rm(maWithExposure)

  #3.(try with the option that does the separate prior on biased effect estimand
  resFile = file.path(cachePath, sprintf("fittedModel-separate-%s-%s.rds", trueRR, repId))
  if(file.exists(resFile)){
    maWithExposure = readRDS(resFile)
  }else{
    settings = generateBayesianHMAsettings(globalExposureEffectPriorMean = c(log(trueRR)),
                                           globalExposureEffectPriorStd = c(10.0),
                                           exposureEffectCount = 1,
                                           separateExposurePrior = TRUE, ...)
    maWithExposure = computeHierarchicalMetaAnalysis(data = metaLPs,
                                                     settings = settings,
                                                     seed = seed)
    saveRDS(maWithExposure, file = resFile)
  }

  estRow = nrow(maWithExposure)
  if(exportAsString){
    res$separate = estimateAsString(maWithExposure[estRow,2:4])
  }else{
    res = bind_rows(res, as.data.frame(estimatesAsVector(maWithExposure[estRow,2:4], label = "separate")))
  }
  rm(maWithExposure)

  #4. two-stage approach
  resFile = file.path(cachePath, sprintf("fittedModel-twoStage-%s-%s.rds", trueRR, repId))
  if(file.exists(resFile)){
    adjustedMainEffectSamps= readRDS(resFile)
  }else{
    settings = generateBayesianHMAsettings(includeExposureEffect = FALSE, ...)
    maNCsOnly = computeHierarchicalMetaAnalysis(data = metaLPs[1:(length(metaLPs)-1)],
                                                seed = seed,
                                                settings = settings)
    traces = attr(maNCsOnly, "traces")

    # newBiasSamps = rnorm(nrow(traces),
    #                      mean = traces[,"outcome.mean"] + traces[,"source.mean"],
    #                      sd = 1/sqrt(traces[,"outcome.scale"]))
    newBiasSamps = traces[,"outcome.mean"] + traces[,"source.mean"] # try the mean bias samples only...

    maExposure = computeBayesianMetaAnalysis(data = metaLPs[[length(metaLPs)]],
                                             seed = seed,
                                             chainLength = 300000,
                                             burnIn = 5e+04)
    tracesExposure = attr(maExposure, "traces")
    mainEffectSamps = tracesExposure[,1]
    adjustedMainEffectSamps = mainEffectSamps - newBiasSamps

    saveRDS(adjustedMainEffectSamps, file = resFile)
  }
  if(exportAsString){
    res$twoStage = estimateAsString(summarizeChain(adjustedMainEffectSamps)[2:4])
  }else{
    res = bind_rows(res, as.data.frame(estimatesAsVector(summarizeChain(adjustedMainEffectSamps)[2:4],
                                                         label = "twoStage")))
  }
  rm(adjustedMainEffectSamps)

  res$id = repId
  if(exportAsString){
    res = data.frame(res)
  }else{
    res$RR = trueRR
  }
  return(res)

}
