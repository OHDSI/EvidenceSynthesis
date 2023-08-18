library(testthat)
#library(EvidenceSynthesis)

## load the example data (list of likelihood profiles across outcomes and data sources)
data("hmaLikelihoodList")

test_that('data model construction function', {
  exDataModel = EvidenceSynthesis:::constructDataModel(hmaLikelihoodList[[1]])

  expect_type(exDataModel, "S4")
})

test_that('construct ArrayList of data models', {
  allDataModels <- rJava::.jnew("java.util.ArrayList")
  for(i in 1:length(hmaLikelihoodList)){
    thisDataModel = EvidenceSynthesis:::constructDataModel(hmaLikelihoodList[[i]])
    allDataModels$add(rJava::.jcast(thisDataModel,
                                    "org.ohdsi.metaAnalysis.DataModel"))
  }

  expect_equal(allDataModels$size(), 3)

  allDataModels <- rJava::.jcast(allDataModels, "java.util.List")

  expect_equal(allDataModels$size(), 3)
})


test_that('run hierarchical meta analysis', {
  estimates = EvidenceSynthesis::computeHierarchicalMetaAnalysis(data = hmaLikelihoodList,
                                                                 chainLength = 110000,
                                                                 burnIn = 1e+04,
                                                                 subSampleFrequency = 10,
                                                                 seed = 666)

  expect_type(estimates, "list")
  expect_named(estimates, c("mean", "LB", "UB", "se", "parameter"))
  expect_type(attr(estimates, "traces"), "double")

})








exDataModel = EvidenceSynthesis:::constructDataModel(hmaLikelihoodList[[1]])
