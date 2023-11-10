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

