# a stupid simple script to test the hierarchical meta analysis function

## construct a list of grid-like profile likelihoods
## default: the last sets of LPs in the constructed list represents the outcome of interest
## example here: 3 outcomes in total, 4 data sources

## read in data
dataModelList = list()
for(i in 1:3){
  dataModelList[[i]] = as.data.frame(readr::read_csv(sprintf("../ForDavid/grids_example_%s.csv", i)))
}

## construct data model for one element in the list
exDataModel = EvidenceSynthesis:::constructDataModel(dataModelList[[1]])

## construct a list of data model objects
allDataModels <- rJava::.jnew("java.util.ArrayList")
for(i in 1:length(dataModelList)){
  thisDataModel = EvidenceSynthesis:::constructDataModel(dataModelList[[i]])
  allDataModels$add(rJava::.jcast(thisDataModel,
                                  "org.ohdsi.metaAnalysis.DataModel"))
}

## cast it to a List (because rJava refuses to do that itself...)
allDataModels <- rJava::.jcast(allDataModels, "java.util.List")

## construct the configuration object
hmaConfiguration = rJava::.jnew("org.ohdsi.metaAnalysis.HierarchicalMetaAnalysis$HierarchicalMetaAnalysisConfiguration")

## construct the `HierarchicalMetaAnalysis` object
hmaObject = rJava::.jnew(
  "org.ohdsi.metaAnalysis.HierarchicalMetaAnalysis",
  allDataModels,
  hmaConfiguration
)

## try the whole function
estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = dataModelList,
                                                               seed = 666)


## try an example list of likelihood profiles extracted from LegendT2dm class CES

#dataModelList = readRDS('data/likelihoodProfileLists.rda')
data("likelihoodProfileLists")

estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = likelihoodProfileLists,
                                                               seed = 666,
                                                               chainLength = 10000,
                                                               burnIn = 100)

## try the option of negative controls only (no exposure effect)
estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = likelihoodProfileLists,
                                                               includeExposureEffect = FALSE,
                                                               seed = 666,
                                                               chainLength = 10000,
                                                               burnIn = 100)


#### Simulate with negative controls across multiple databases
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

createApproximations <- function(populations, approximation) {
  fitModelInDatabase <- function(population, approximation) {
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                              data = population,
                                              modelType = "cox"
    )
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData,
                                           fixedCoefficients = c(approximation != "normal")
    )
    approximation <- approximateLikelihood(cyclopsFit, "x", approximation = approximation)
    return(approximation)
  }
  data <- lapply(populations, fitModelInDatabase, approximation = approximation)
  if (approximation != "adaptive grid") {
    data <- do.call("rbind", data)
  }
  return(data)
}

trueRR = 20
metaPopulations = simulateMetaAnalysisWithNegativeControls(meanExposureEffect = log(trueRR),
                                                           mNegativeControls = 50,
                                                           sitePop = 1000,
                                                           treatedFraction = 0.3)
metaLPs = lapply(metaPopulations, createApproximations, "adaptive grid")

## fit a joint model with main exposure effect included:
maWithExposure = computeHierarchicalMetaAnalysis(data = metaLPs,
                                                 seed = 666)
# trueRR = 2
# exp(maWithExposure[6,1:3])
# mean        LB       UB
# 1.925166 0.4293734 14.54448

# trueRR = 10
# exp(maWithExposure[c(2,6), 1:3])
# mean        LB        UB
# 1.472945 0.9310029  2.303401
# 2.796862 0.6341481 20.977832

# trueRR = 20
# exp(maWithExposure[6, 1:3])
# mean       LB       UB
# 3.825785 1.255261 14.18314


## fit two separate models: (1) NCs only (2) main effect meta analysis
trueRR = 2
metaPopulations = simulateMetaAnalysisWithNegativeControls(meanExposureEffect = log(trueRR),
                                                           mNegativeControls = 50,
                                                           sitePop = 1000,
                                                           treatedFraction = 0.3)
metaLPs = lapply(metaPopulations, createApproximations, "adaptive grid")

maNCsOnly = computeHierarchicalMetaAnalysis(data = metaLPs,
                                            seed = 666,
                                            includeExposureEffect = FALSE)

## probably need to force source.mean = 0???
# > maNCsOnly
# mean         LB        UB        se     parameter
# 1  3.49568277  1.4444395 4.9911600 0.9653108           tau
# 2 -0.01835586 -0.2820391 0.2304304 0.1472115  outcome.mean
# 3  4.61433971  2.1109409 7.6662797 1.6085099 outcome.scale
# 4  0.37279051 -0.1951348 0.8173173 0.2516555   source.mean
# 5  2.95126975  0.6622224 5.6526388 1.4039522  source.scale

traces = attr(maNCsOnly, "traces")
biasParams = traces[,which(colnames(traces) %in% c("outcome.mean", "outcome.scale"))]

maExposure = computeBayesianMetaAnalysis(data = metaLPs[[length(metaLPs)]],
                                         seed = 666)
tracesExposure = attr(maExposure, "traces")
mainEffectSamps = tracesExposure[,1]



