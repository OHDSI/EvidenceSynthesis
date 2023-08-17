# a stupid simple script to test the hierarchical meta analysis function

## construct a list of grid-like profile likelihoods
## default: the last sets of LPs in the constructed list represents the outcome of interest
## example here: 3 outcomes in total, 4 data sources

dataModelList = list()
for(i in 1:3){
  dataModelList[[i]] = as.data.frame(readr::read_csv(sprintf("ForDavid/grids_example_%s.csv", i)))
}


estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = dataModelList,
                                                               seed = 666)
