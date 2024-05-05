# a stupid simple script to test the hierarchical meta analysis function

## construct a list of grid-like profile likelihoods
## default: the last sets of LPs in the constructed list represents the outcome of interest
## example here: 3 outcomes in total, 4 data sources


rJava::.jinit(parameters="-Xmx32g", force.init = TRUE)
options(java.parameters = c("-Xms200g", "-Xmx200g"))

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

settings = generateBayesianHMAsettings(chainLength = 10000,
                                       burnIn = 100)

estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = likelihoodProfileLists,
                                                               settings = settings,
                                                               seed = 666)

## try the option of negative controls only (no exposure effect)

settings$includeExposureEffect = FALSE
estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = likelihoodProfileLists,
                                                               settings = settings,
                                                               seed = 666)

## try with (almost) forcing source.mean to be 0 (with tiny prior std)

settings$secondaryEffectPriorStd = 0.0001
estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = likelihoodProfileLists,
                                                               settings = settings,
                                                               seed = 666)



#### Try with another LegendT2dm example -----
ooiLps = readRDS("extras/ooiLpList.rds")[1:2] # first two OOI
ncLps = readRDS("extras/ncLpsList.rds")[1:10] # first 10 NC
metaLPs = c(ncLps, ooiLps)

settings = generateBayesianHMAsettings(globalExposureEffectPriorStd = c(10.0),
                                       exposureEffectCount = 2)

metaDefault = computeHierarchicalMetaAnalysis(data = metaLPs,
                                             settings = settings,
                                             seed = 666)
# seems to work??

settings$separateExposurePrior = TRUE
metaSep = computeHierarchicalMetaAnalysis(data = metaLPs,
                                          settings = settings,
                                          seed = 666)

## try with new model that accepts multiple priors for multiple OOIs

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
maWithExposure = computeHierarchicalMetaAnalysis(data = metaLPs,
                                                 secondaryEffectPriorStd = 0.0001,
                                                 seed = 666)
maWithExposure

## (try using a "friendlier" prior for the meta-analytic effect)
maWithExposure = computeHierarchicalMetaAnalysis(data = metaLPs,
                                                 # includeSourceEffect = FALSE,
                                                 # secondaryEffectPriorStd = 0.0001,
                                                 globalExposureEffectPrior = c(log(trueRR), 2.0),
                                                 seed = 666)#42)
maWithExposure
exp(maWithExposure[6,2:4])

## (try with the option that does the separate prior on biased effect estimand)
maWithExposureSep = computeHierarchicalMetaAnalysis(data = metaLPs,
                                                 #secondaryEffectPriorStd = 0.0001,
                                                 #includeSourceEffect = FALSE,
                                                 globalExposureEffectPrior = c(log(trueRR), 2.0),
                                                 separateExposurePrior = TRUE,
                                                 seed = 42)#42)
maWithExposureSep
#traces = attr(maWithExposureSep, "traces")
exp(maWithExposureSep[6,2:4])

# trueRR = 2
# > exp(maWithExposure[6,1:3])
# mean        LB       UB
# 2.065413 1.111241 2.818157

## set prior center at trueRR:
# > exp(maWithExposure[6,1:3])
# mean       LB       UB
# 1.8921 1.144289 2.713055

# trueRR = 5
# exp(maWithExposure[6,1:3])
# mean        LB      UB
# 2.586767 0.6052171 9.13415

## set prior center at trueRR:
# > exp(maWithExposure[6,1:3])
# mean       LB       UB
# 4.538813 1.692695 15.75625

# trueRR = 10
# > exp(maWithExposure[6, 1:3])
#    mean       LB       UB
# 8.504463 6.853034 12.31833

## set prior center at trueRR:
# > exp(maWithExposure[6, 1:3])
#   mean       LB       UB
# 7.069216 5.316982 9.481317

# trueRR = 20
# exp(maWithExposure[6, 1:3])
#    mean     LB       UB
# 10.96739 8.1194 14.17135

## set prior center at trueRR:
# exp(maWithExposure[6, 1:3])
#   mean       LB       UB
# 17.30228 12.06975 24.01206


## fit two separate models: (1) NCs only (2) main effect meta analysis
# trueRR = 10
# metaPopulations = simulateMetaAnalysisWithNegativeControls(meanExposureEffect = log(trueRR),
#                                                            mNegativeControls = 50,
#                                                            sitePop = 1000,
#                                                            treatedFraction = 0.3)
# metaLPs = lapply(metaPopulations, createApproximations, "adaptive grid")

maNCsOnly = computeHierarchicalMetaAnalysis(data = metaLPs[1:(length(metaLPs)-1)],
                                            #secondaryEffectPriorStd = 0.0001,
                                            seed = 666,
                                            #includeSourceEffect = FALSE,
                                            includeExposureEffect = FALSE,
                                            showProgressBar = TRUE)
maNCsOnly

## forced source.mean = 0 (numerically)
# > maNCsOnly
#          mean            LB           UB           se     parameter
# 1 3.219338e+00  1.3106086591 5.1953496335 1.0665086487           tau
# 2 4.333239e-01  0.1465749126 0.7495771992 0.1589051920  outcome.mean
# 3 4.944244e+00  2.2291183323 8.1694654437 1.6766489053 outcome.scale
# 4 4.627727e-06 -0.0001986662 0.0001972235 0.0001133314   source.mean
# 5 3.181404e+00  0.8619298496 6.3158341143 1.5088371886  source.scale

traces = attr(maNCsOnly, "traces")
# biasParams = traces[,which(colnames(traces) %in% c("outcome.mean", "outcome.scale"))]
# newBiasSamps = rnorm(nrow(biasParams), mean = biasParams[,1], sd = sqrt(1/biasParams[,2]))
newBiasSamps = rnorm(nrow(traces),
                     mean = traces[,"outcome.mean"] + traces[,"source.mean"],
                     sd = 1/sqrt(traces[,"outcome.scale"]))

maExposure = computeBayesianMetaAnalysis(data = metaLPs[[length(metaLPs)]],
                                         seed = 42)#666) ## Nov 26, 2023: want to have non-zero prior mean/center!!
tracesExposure = attr(maExposure, "traces")
mainEffectSamps = tracesExposure[,1]

## substract using predictive samples of bias
adjustedMainEffectSamps = mainEffectSamps - newBiasSamps

## try another way: subtract mean bias directly...
#adjustedMainEffectSamps = mainEffectSamps - biasParams[,1]

# trueRR = 2
summarizeChain(adjustedMainEffectSamps)
# > EvidenceSynthesis:::summarizeChain(adjustedMainEffectSamps)
#             lower      upper
# 0.7881223 -8.9545905 10.9922600  4.9457848
exp(c(median(adjustedMainEffectSamps), quantile(adjustedMainEffectSamps, c(.025, 0.975))))
#                    2.5%        97.5%
# 2.249790e+00 8.471367e-05 4.174356e+04

# trueRR = 5
EvidenceSynthesis:::summarizeChain(adjustedMainEffectSamps)
# > EvidenceSynthesis:::summarizeChain(adjustedMainEffectSamps)
#              lower    upper
# 0.8886139 -9.6428057 11.5253916  5.2448567
exp(c(median(adjustedMainEffectSamps), quantile(adjustedMainEffectSamps, c(.025, 0.975))))
#                    2.5%        97.5%
#   2.454859e+00 5.652003e-05 9.692100e+04

## try subtracting the mean bias term only (without drawing another bias sample...)
EvidenceSynthesis:::summarizeChain(adjustedMainEffectSamps)
#              lower      upper
# 1.59081439 1.42651281 1.75591734 0.08371495
exp(c(median(adjustedMainEffectSamps), quantile(adjustedMainEffectSamps, c(.025, 0.975))))
#                2.5%    97.5%
#   4.905457 4.172667 5.810921

# trueRR = 10
EvidenceSynthesis:::summarizeChain(adjustedMainEffectSamps)
# > EvidenceSynthesis:::summarizeChain(adjustedMainEffectSamps)
#              lower      upper
# 1.004835 -10.281391  12.076646   5.443380
exp(c(median(adjustedMainEffectSamps), quantile(adjustedMainEffectSamps, c(.025, 0.975))))
#                   2.5%        97.5%
# 2.797849e+00 2.508974e-05 1.466959e+05

## try subtracting the mean bias term only (without drawing another bias sample...)
EvidenceSynthesis:::summarizeChain(adjustedMainEffectSamps)
#              lower     upper
# 2.5917708 1.9895421 3.4382692 0.4147019
exp(c(median(adjustedMainEffectSamps), quantile(adjustedMainEffectSamps, c(.025, 0.975))))
#                 2.5%     97.5%
#   12.263272  7.801114 37.852637


# trueRR = 20
EvidenceSynthesis:::summarizeChain(adjustedMainEffectSamps)
# > EvidenceSynthesis:::summarizeChain(adjustedMainEffectSamps)
#             lower     upper
#  3.116734 -8.214927 15.178261  5.974132
exp(c(median(adjustedMainEffectSamps), quantile(adjustedMainEffectSamps, c(.025, 0.975))))
#                      2.5%        97.5%
#  2.277527e+01 2.117537e-04 3.246991e+06

## try subtracting the mean bias term only (without drawing another bias sample...)
EvidenceSynthesis:::summarizeChain(adjustedMainEffectSamps)
#              lower     upper
# 3.1442300 1.8461987 5.0831462 0.8721749
exp(c(median(adjustedMainEffectSamps), quantile(adjustedMainEffectSamps, c(.025, 0.975))))
#               2.5%      97.5%
# 18.591307   7.511103 256.619320



