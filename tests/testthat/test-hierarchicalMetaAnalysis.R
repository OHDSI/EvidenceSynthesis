library(testthat)
# library(EvidenceSynthesis)

## test simulation function
test_that("simulate with negative controls on multiple databases", {
  metaPopulations <- simulateMetaAnalysisWithNegativeControls(
    meanExposureEffect = log(2.0),
    mNegativeControls = 10,
    sitePop = 1000,
    treatedFraction = 0.5,
    minBackgroundHazard = 0.2,
    maxBackgroundHazard = 0.2,
    nStrata = 1,
    meanBias = 0.5, biasStd = 0.1,
    siteEffectStd = 0.15,
    nSites = 5
  )
  expect_type(metaPopulations, "list")
  expect_type(attr(metaPopulations, "hyperParameters"), "list")
})

## load the example data (list of likelihood profiles across outcomes and data sources)
data("hmaLikelihoodList")

test_that("data model construction function", {
  exDataModel <- EvidenceSynthesis:::constructDataModel(hmaLikelihoodList[[1]])

  expect_type(exDataModel, "S4")
})

test_that("construct ArrayList of data models", {
  allDataModels <- rJava::.jnew("java.util.ArrayList")
  for (i in 1:length(hmaLikelihoodList)) {
    thisDataModel <- EvidenceSynthesis:::constructDataModel(hmaLikelihoodList[[i]])
    allDataModels$add(rJava::.jcast(
      thisDataModel,
      "org.ohdsi.metaAnalysis.DataModel"
    ))
  }

  expect_equal(allDataModels$size(), 3)

  allDataModels <- rJava::.jcast(allDataModels, "java.util.List")

  expect_equal(allDataModels$size(), 3)
})

test_that("build hma settings list", {
  settings <- EvidenceSynthesis::generateBayesianHMAsettings(
    globalExposureEffectPriorMean = c(2.0, 3.0),
    globalExposureEffectPriorStd = c(1.0, 10.0),
    exposureEffectCount = 2
  )
  expect_type(settings, "list")
  expect_true(length(settings$globalExposureEffectPriorMean) == 2)
  expect_true(length(settings$globalExposureEffectPriorStd) == 2)
})


test_that("run hierarchical meta analysis", {
  settings <- EvidenceSynthesis::generateBayesianHMAsettings(
    chainLength = 110000,
    burnIn = 1e+04,
    subSampleFrequency = 10
  )
  estimates <- EvidenceSynthesis::computeHierarchicalMetaAnalysis(
    data = hmaLikelihoodList,
    settings = settings,
    seed = 666
  )

  expect_type(estimates, "list")
  expect_named(estimates, c("mean", "median", "LB", "UB", "se", "parameter"))
  expect_type(attr(estimates, "traces"), "double")
})



## load another (bigger & different) list of profile likelihoods
## (from LegendT2dm class CES)
data("likelihoodProfileLists")

test_that("run hierarchical meta analysis on bigger data", {
  settings <- EvidenceSynthesis::generateBayesianHMAsettings(
    chainLength = 110000,
    burnIn = 1e+04,
    subSampleFrequency = 10
  )
  estimates <- EvidenceSynthesis::computeHierarchicalMetaAnalysis(
    data = likelihoodProfileLists,
    settings = settings,
    seed = 666
  )

  expect_type(estimates, "list")
  expect_named(estimates, c("mean", "median", "LB", "UB", "se", "parameter"))
  expect_type(attr(estimates, "traces"), "double")
})


test_that("run hierarchical meta analysis without main exposure effect", {
  settings <- EvidenceSynthesis::generateBayesianHMAsettings(
    chainLength = 110000,
    burnIn = 1e+04,
    subSampleFrequency = 10,
    includeExposureEffect = FALSE
  )

  likelihoodProfileLists_ncOnly <- likelihoodProfileLists[1:10]

  estimates <- EvidenceSynthesis::computeHierarchicalMetaAnalysis(
    data = likelihoodProfileLists_ncOnly,
    settings = settings,
    seed = 666
  )

  expect_type(estimates, "list")
  expect_named(estimates, c("mean", "median", "LB", "UB", "se", "parameter"))
  expect_type(attr(estimates, "traces"), "double")
  expect_false("exposure" %in% estimates$parameter)
})

test_that("run hierarchical meta analysis without two main exposure effects", {
  settings <- EvidenceSynthesis::generateBayesianHMAsettings(
    chainLength = 110000,
    burnIn = 1e+04,
    subSampleFrequency = 10,
    includeExposureEffect = TRUE,
    exposureEffectCount = 2,
    globalExposureEffectPriorMean = c(1.0, 1.5),
    globalExposureEffectPriorStd = c(2.0, 10.0)
  )

  likelihoodProfileLists_twoExposures <- c(
    likelihoodProfileLists[1:10],
    likelihoodProfileLists[11],
    likelihoodProfileLists[11]
  )

  estimates <- EvidenceSynthesis::computeHierarchicalMetaAnalysis(
    data = likelihoodProfileLists_twoExposures,
    settings = settings,
    seed = 666
  )

  expect_type(estimates, "list")
  expect_named(estimates, c("mean", "median", "LB", "UB", "se", "parameter"))
  expect_type(attr(estimates, "traces"), "double")
  expect_true("exposure1" %in% estimates$parameter)
  expect_true("exposure2" %in% estimates$parameter)
})
