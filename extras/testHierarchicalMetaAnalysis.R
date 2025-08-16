# a stupid simple script to test the hierarchical meta analysis function

## construct a list of grid-like profile likelihoods
## default: the last sets of LPs in the constructed list represents the outcome of interest
## example here: 3 outcomes in total, 4 data sources


rJava::.jinit(parameters="-Xmx32g", force.init = TRUE)
options(java.parameters = c("-Xms200g", "-Xmx200g"))

#### step-by-step construction test for rJava interface ----

# ## read in data
# dataModelList = list()
# for(i in 1:3){
#   dataModelList[[i]] = as.data.frame(readr::read_csv(sprintf("../ForDavid/grids_example_%s.csv", i)))
# }
#
# ## construct data model for one element in the list
# exDataModel = EvidenceSynthesis:::constructDataModel(dataModelList[[1]])
#
# ## construct a list of data model objects
# allDataModels <- rJava::.jnew("java.util.ArrayList")
# for(i in 1:length(dataModelList)){
#   thisDataModel = EvidenceSynthesis:::constructDataModel(dataModelList[[i]])
#   allDataModels$add(rJava::.jcast(thisDataModel,
#                                   "org.ohdsi.metaAnalysis.DataModel"))
# }
#
# ## cast it to a List (because rJava refuses to do that itself...)
# allDataModels <- rJava::.jcast(allDataModels, "java.util.List")
#
# ## construct the configuration object
# hmaConfiguration = rJava::.jnew("org.ohdsi.metaAnalysis.HierarchicalMetaAnalysis$HierarchicalMetaAnalysisConfiguration")
#
# ## construct the `HierarchicalMetaAnalysis` object
# hmaObject = rJava::.jnew(
#   "org.ohdsi.metaAnalysis.HierarchicalMetaAnalysis",
#   allDataModels,
#   hmaConfiguration
# )
#
# ## try the whole function
# estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = dataModelList,
#                                                                seed = 666)
# print(estimates) # show summary table


#### try an example list of likelihood profiles extracted from LegendT2dm class CES ----
data("likelihoodProfileLists")
settings = generateBayesianHMAsettings(chainLength = 1e6,
                                       burnIn = 1e4,
                                       separateExposurePrior = TRUE)
estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = likelihoodProfileLists,
                                                               settings = settings,
                                                               seed = 666)
print(estimates)


## try the option of negative controls only (no exposure effect)
settings$includeExposureEffect = FALSE
estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = likelihoodProfileLists,
                                                               settings = settings,
                                                               seed = 666)
print(estimates)

## try option of HMC
settings$useHMC = TRUE
settings$includeExposureEffect = TRUE
estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = likelihoodProfileLists,
                                                               settings = settings,
                                                               seed = 666)
print(estimates)

## try option of heteroscedastic variances across data sources
settings$useHMC = FALSE
estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = likelihoodProfileLists,
                                                               settings = settings,
                                                               seed = 666)
print(estimates)



#### Test with another LegendT2dm example -----
ooiLps = readRDS("extras/ooiLpList.rds")[1:2] # first two OOI
ncLps = readRDS("extras/ncLpsList.rds")[1:10] # first 10 NC
metaLPs = c(ncLps, ooiLps)

settings = generateBayesianHMAsettings(globalExposureEffectPriorStd = c(10.0),
                                       exposureEffectCount = 2)

## default setting
metaDefault = computeHierarchicalMetaAnalysis(data = metaLPs,
                                             settings = settings,
                                             seed = 666)

## separable prior (usually performs better)
settings$separateExposurePrior = TRUE
metaSep = computeHierarchicalMetaAnalysis(data = metaLPs,
                                          settings = settings,
                                          seed = 666)

## new model that accepts multiple priors for multiple OOIs
settings = generateBayesianHMAsettings(globalExposureEffectPriorMean = c(1.0, 2.0),
                                       globalExposureEffectPriorStd = c(2.0, 10.0),
                                       exposureEffectCount = 2)
metaTest = computeHierarchicalMetaAnalysis(data = metaLPs,
                                           settings = settings,
                                           seed = 666)

#### Simulate with negative controls across multiple databases ----


## start simulations: joint model
## UPDATE: increase background rates
trueRR = 1
metaPopulations = simulateMetaAnalysisWithNegativeControls(meanExposureEffect = log(trueRR),
                                                           mNegativeControls = 50,
                                                           sitePop = 10000,
                                                           treatedFraction = 0.5,
                                                           minBackgroundHazard = 0.05,
                                                           maxBackgroundHazard = 0.05,
                                                           nStrata = 1,
                                                           meanBias = 0.5, biasStd = 0.1,
                                                           siteEffectStd = 0.15)
metaLPs = lapply(metaPopulations, createApproximations, "adaptive grid")

## fit a joint model with main exposure effect included:
settings = generateBayesianHMAsettings(useHMC = TRUE)
maWithExposure = computeHierarchicalMetaAnalysis(data = metaLPs,
                                                 settings = settings,
                                                 seed = 666)
# summary of results
maWithExposure
# check out 95% CI
exp(maWithExposure[6,2:4])
# extract the MCMC samples (parameters by the columns...)
traces = attr(maWithExposure, "traces")
# check out traceplot for something
plot(traces[,"source.mean"]+traces[,"outcome.mean"], type = "l") ## model weak identifiability, so need to check source+outcome mean
plot(traces[,"exposure1"], type = "l")

## try using an informative prior
settings$globalExposureEffectPriorMean = log(trueRR)
settings$globalExposureEffectPriorStd = 2.0
maWithExposure2 = computeHierarchicalMetaAnalysis(data = metaLPs,
                                                 settings = settings,
                                                 seed = 666)
maWithExposure2
# check out 95% CI
exp(maWithExposure2[6,2:4])
# extract the MCMC samples
traces2 = attr(maWithExposure2, "traces")
# traceplot
plot(traces2[,"exposure1"], type = "l")

## try separable prior
settings$separateExposurePrior = TRUE
maWithExposureSep = computeHierarchicalMetaAnalysis(data = metaLPs,
                                                    settings = settings,
                                                    seed = 42)
maWithExposureSep
#traces = attr(maWithExposureSep, "traces")
exp(maWithExposureSep[6,2:4]) # check out 95% CI
