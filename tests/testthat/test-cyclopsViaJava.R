library(Cyclops)

library(EvidenceSynthesis)
library(testthat)

test_that("Simple Cyclops example", {

  skip("Requires Cyclops-jni")

  dobson <- data.frame(
    counts = c(18,17,15,20,10,20,25,13,12),
    outcome = gl(3,1,9),
    treatment = gl(3,3))

  data <- createCyclopsData(counts ~ outcome + treatment, data = dobson,
                            modelType = "pr")

  fit <- fitCyclopsModel(data,
                         prior = createPrior("none"),
                         control = createControl(noiseLevel = "silent"))

  instance <- cacheCyclopsModelForJava(fit)

  libraryFileName <- normalizePath(system.file("libs", .Platform$r_arch,
                                               paste0("Cyclops", .Platform$dynlib.ext),
                                               package = "Cyclops"))

  obj <- rJava::.jnew("dr.inference.regression.RegressionInCyclops",
                      libraryFileName, instance)

  tolerance <- 1e-10

  expect_true(fit$log_likelihood == obj$getLogLikelihood())

  obj$setBeta(as.integer(0), as.numeric(0))

  expect_false(fit$log_likelihood == obj$getLogLikelihood())

  grad <- obj$getLogLikelihoodGradient()

  expect_lt(grad[1], -1)

  obj$findMode()

  expect_equal(fit$log_likelihood, obj$getLogLikelihood())

  grad2 <- obj$getLogLikelihoodGradient()

  expect_equal(grad2[1], 0, tolerance = 1E-3)
})
