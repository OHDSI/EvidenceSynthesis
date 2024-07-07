library(EvidenceSynthesis)
library(Cyclops)
library(dplyr)
library(survival)

# Read interaction data

data <- readRDS("extras/data/InteractionDataForMarc.rds")
# %>% filter(siteId == 1)

outcomeOfInterest <- 77

cyclopsFits <- lapply(unique(data$siteId), function(id) {
  subset <- data %>% filter(outcomeId == outcomeOfInterest,
                            siteId == id)


  denseData <- createCyclopsData(Surv(survivalTime, y) ~ treatment + subgroup + subgroup * treatment + strata(stratumId),
                                 modelType = "cox", data = subset)

  denseFit <- fitCyclopsModel(denseData, prior = createPrior("none"))
  # TODO Can avoid fitting issues when mode does not exist / is ill-conditioned by:
  #  1. add substantial regularization above (only used to find initial "mode"), or
  #  2. setting iteration steps to 0 or something small

  denseInstance <- cacheCyclopsModelForJava(denseFit)

  # sparseData <- createCyclopsData(Surv(survivalTime, y) ~ strata(stratumId),
  #                                  indicatorFormula = ~ treatment + subgroup + subgroup * treatment,
  #                                  modelType = "cox", data = subset)
  #
  # sparseFit <- fitCyclopsModel(sparseData, prior = createPrior("none"))
  #
  # sparseInstance <- cacheCyclopsModelForJava(sparseFit)

  list(siteId = id,
       data = subset,
       fit =
         denseFit,
       # sparseFit,
       instance =
         denseInstance
       # sparseInstance
  )
})

cyclopsInstances <- unlist(lapply(cyclopsFits, function(fit) { fit$instance} ))

cyclopsLibraryFileName <- normalizePath(system.file("libs", .Platform$r_arch, # TODO make Cyclops package function
                                                    paste0("Cyclops", .Platform$dynlib.ext),
                                                    package = "Cyclops"))

# Create java objects

dataModelList <- rJava::.jnew("java.util.ArrayList")
lapply(cyclopsFits, function(object) {
  likelihood <- rJava::.jnew(
    "dr.inference.regression.CyclopsRegressionModel",
    paste0("site", object$siteId),
    cyclopsLibraryFileName,
    as.integer(object$instance),
    as.integer(length(coef(object$fit))),
    as.logical(TRUE)
  )

  dataModelList$add(rJava::.jcast(likelihood, "org.ohdsi.metaAnalysis.DataModel"))
})

config <- rJava::.jnew("org.ohdsi.metaAnalysis.MultivariableHierarchicalMetaAnalysis$HierarchicalMetaAnalysisConfiguration")

analysis <- rJava::.jnew("org.ohdsi.metaAnalysis.MultivariableHierarchicalMetaAnalysis",
                         rJava::.jcast(dataModelList, "java.util.List"),
                         config)

chainLength <- 110000
burnIn <- 10000
subSampleFrequency <- 100
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
  runner$run() # TODO Use HMC for (substantially?) improved mixing
)

runner$processSamples()

# EXTRA BELOW

parameterNames <- runner$getParameterNames()
trace <- runner$getTrace(as.integer(3))
traces <- matrix(ncol = length(parameterNames) - 2, nrow = length(trace))
colnames(traces) <- parameterNames[-c(1:2)]
traces[, 1] <- trace
for (i in 4:length(parameterNames)) {
  trace <- runner$getTrace(as.integer(i))
  traces[, i - 2] <- trace
}

plot(coda::as.mcmc(traces[,c("mean1", "mean2", "mean3")]))
coda::effectiveSize(traces)



