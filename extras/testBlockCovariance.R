## Code to test block diagonal covariance matrix for HMA

# rJava::.jinit(parameters="-Xmx32g", force.init = TRUE)
# options(java.parameters = c("-Xms200g", "-Xmx200g"))

#### build data model ----
## read in data
dataModelList = list()
for(i in 1:3){
  dataModelList[[i]] = as.data.frame(readr::read_csv(sprintf("extras/DM_example/grids_example_%s.csv", i)))
}

## add another set of LPs (fake example)
# dataModelList[[4]] = dataModelList[[3]]

#### create the settings object ----
## set blockCovariance = TRUE to try block diagonal covariance matrix
testSettings = generateBayesianHMAsettings(blockCovariance = TRUE, # set this to FALSE runs ok
                                           chainLength = 11000, # can increase these
                                           burnIn = 1e3,
                                           subSampleFrequency = 10,
                                           exposureEffectCount = 1)


#### try running HMA ----
estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = dataModelList,
                                                               settings = testSettings,
                                                               seed = 666)
print(estimates) # show summary table

## check out MCMC samples
mcmcs = attr(estimates, "traces")
str(mcmcs)
