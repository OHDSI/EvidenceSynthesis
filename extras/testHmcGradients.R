# test script for gradients not matching error

library(EvidenceSynthesis)

rJava::.jinit(parameters="-Xmx32g", force.init = TRUE)
options(java.parameters = c("-Xms200g", "-Xmx200g"))

## read in test data 1
data("likelihoodProfileLists")

settings = generateBayesianHMAsettings(chainLength = 50000,
                                       burnIn = 1000)
### change tolerance to 1E-1 (in HmcOptions), this can run
### original setting 1E-3 --> throws error

estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = likelihoodProfileLists,
                                                               settings = settings,
                                                               seed = 666)


## test data 2 (larger test set)
## based on patient-level trajectory simulation
nSites = 10
sitePop = 10000
mNegativeControls = 50
meanBias = 0.5
biasStd = 0.1
meanSiteEffect = 0
siteEffectStd = 0.15
seed = 42

### (1) simulate data
simulateAndApproximate <- function(trueRR){
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

  return(metaLPs)
}

metaLPs = simulateAndApproximate(trueRR = 2)

### (2) fit HMA model
settings = generateBayesianHMAsettings(globalExposureEffectPriorMean = c(0),
                                       globalExposureEffectPriorStd = c(10.0),
                                       exposureEffectCount = 1,
                                       chainLength = 2000, # try just a few iterations
                                       burnIn = 100)
maWithExposure = computeHierarchicalMetaAnalysis(data = metaLPs,
                                                 settings = settings,
                                                 seed = seed)
### Oops...
### Error in .jcall("RJavaTools", "Ljava/lang/Object;", "invokeMethod", cl,  :
### java.lang.RuntimeException: Gradients do not match:
