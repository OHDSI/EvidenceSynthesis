# run this to test if `computeBayesianMetaAnalysis` works
library(EvidenceSynthesis)
library(Cyclops)
library(dplyr)
library(survival)

# Run bladder x 4 example (to match `MultivariateHierarchicalAnalysis.main`) using `MultivariableCoxPartialLikelihood`

reps <- 4

dataModelList <- rJava::.jnew("java.util.ArrayList")
for (i in 1:reps) {
  mode <- coxph(Surv(stop, event) ~ (rx - 1) + size, data = bladder, ties = "breslow")

  data <- rJava::.jnew(
    "org.ohdsi.data.CoxData",
    as.integer(bladder$event),
    as.double(bladder$stop),
    as.double(c(bladder$rx, bladder$size))
  )

  parameter <- rJava::.jnew("dr.inference.model.Parameter$Default", coef(mode))
  likelihood <- rJava::.jnew(
    "org.ohdsi.likelihood.MultivariableCoxPartialLikelihood",
    rJava::.jcast(parameter, "dr.inference.model.Parameter"),
    data$getSortedData()
  )

  dataModelList$add(rJava::.jcast(likelihood, "org.ohdsi.metaAnalysis.DataModel"))
}

cg <- rJava::.jnew("org.ohdsi.metaAnalysis.MultivariableHierarchicalMetaAnalysis$HierarchicalMetaAnalysisConfiguration")

analysis <- rJava::.jnew("org.ohdsi.metaAnalysis.MultivariableHierarchicalMetaAnalysis",
                         rJava::.jcast(dataModelList, "java.util.List"),
                         cg)

chainLength <- 110000
burnIn <- 10000
subSampleFrequency <- 10
seed <- 666
showProgressBar <- TRUE

runner <- rJava::.jnew(
  "org.ohdsi.mcmc.Runner",
  rJava::.jcast(analysis,
                "org.ohdsi.mcmc.Analysis"),
  as.integer(chainLength),
  as.integer(burnIn),
  as.integer(subSampleFrequency),
  as.numeric(seed),
  as.logical(showProgressBar)
)

system.time(
  runner$run()
)

runner$processSamples()

# Run bladder x 4 example (to match `MultivariateHierarchicalAnalysis.main`) using `CyclopsRegressionModel`

cyclopsLibraryFileName <- normalizePath(system.file("libs", .Platform$r_arch, # TODO make Cyclops package function
                                                    paste0("Cyclops", .Platform$dynlib.ext),
                                                    package = "Cyclops"))

dataModelList <- rJava::.jnew("java.util.ArrayList")
handle <- lapply(1:reps, function(i) {

  data <- Cyclops::createCyclopsData(Surv(stop, event) ~ (rx - 1) + size, data = bladder, modelType = "cox")
  fit <- Cyclops::fitCyclopsModel(data)
  instance <- Cyclops::cacheCyclopsModelForJava(fit)

  likelihood <- rJava::.jnew(
    "dr.inference.regression.CyclopsRegressionModel",
    paste0("site", i),
    cyclopsLibraryFileName,
    as.integer(instance),
    as.integer(length(coef(fit))),
    as.logical(TRUE)
  )

  dataModelList$add(rJava::.jcast(likelihood, "org.ohdsi.metaAnalysis.DataModel"))

  return(list(
    data = data,
    fit = fit,
    instance = instance,
    likelihood = likelihood
  ))
})

cg <- rJava::.jnew("org.ohdsi.metaAnalysis.MultivariableHierarchicalMetaAnalysis$HierarchicalMetaAnalysisConfiguration")

analysis <- rJava::.jnew("org.ohdsi.metaAnalysis.MultivariableHierarchicalMetaAnalysis",
                         rJava::.jcast(dataModelList, "java.util.List"),
                         cg)

chainLength <- 110000
burnIn <- 10000
subSampleFrequency <- 10
seed <- 666
showProgressBar <- TRUE

runner <- rJava::.jnew(
  "org.ohdsi.mcmc.Runner",
  rJava::.jcast(analysis,
                "org.ohdsi.mcmc.Analysis"),
  as.integer(chainLength),
  as.integer(burnIn),
  as.integer(subSampleFrequency),
  as.numeric(seed),
  as.logical(showProgressBar)
)

system.time(
  runner$run()
)

runner$processSamples()

# EXTRA BELOW

# parameterNames <- javaAnalysis$getParameterNames()
# trace <- javaAnalysis$getTrace(as.integer(3))
# traces <- matrix(ncol = length(parameterNames) - 2, nrow = length(trace))
# colnames(traces) <- parameterNames[-c(1:2)]
# traces[, 1] <- trace
# for (i in 4:length(parameterNames)) {
#   trace <- javaAnalysis$getTrace(as.integer(i))
#   traces[, i - 2] <- trace
# }
#
# plot(coda::as.mcmc(traces))
# coda::effectiveSize(traces)



