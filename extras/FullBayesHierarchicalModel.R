# run this to test if `computeBayesianMetaAnalysis` works
library(EvidenceSynthesis)
library(Cyclops)
library(dplyr)
library(survival)

data <- readRDS("extras/data/InteractionDataForMarc.rds") %>% filter(siteId == 1)

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
       fit = denseFit,
       instance =
         denseInstance
         # sparseInstance
       )
  })

cyclopsInstances <- unlist(lapply(cyclopsFits, function(fit) { fit$instance} ))

cyclopsLibraryFileName <- normalizePath(system.file("libs", .Platform$r_arch, # TODO make Cyclops package function
                                                    paste0("Cyclops", .Platform$dynlib.ext),
                                                    package = "Cyclops"))

# Run a single site (silly example -- just uses Metropolis-Hastings and iid prior)

chainLength <- 110000
burnIn <- 10000
subSampleFrequency <- 100
seed <- 666
showProgressBar <- TRUE

priorMean <- 0
priorSd <- 1

siteId <- 1

regressionModel <- rJava::.jnew(
  "dr.inference.regression.CyclopsRegressionModel",
  paste0("site", siteId),
  cyclopsLibraryFileName,
  as.integer(cyclopsInstances[siteId]),
  as.integer(length(coef(cyclopsFits[[siteId]]$fit))),
  as.logical(TRUE)
)

analsysObject <- rJava::.jnew(
  "org.ohdsi.simpleDesign.CyclopsNormalAnalysis",
  regressionModel,
  as.numeric(priorMean),
  as.numeric(priorSd)
)

singleAnalysis <- rJava::.jnew(
  "org.ohdsi.mcmc.Runner",
  rJava::.jcast(analsysObject,
                "org.ohdsi.mcmc.Analysis"),
  as.integer(chainLength),
  as.integer(burnIn),
  as.integer(subSampleFrequency),
  as.numeric(seed),
  as.logical(showProgressBar)
)

singleAnalysis$setConsoleWidth(getOption("width"))

system.time(
  singleAnalysis$run()
)

singleAnalysis$processSamples()

parameterNames <- singleAnalysis$getParameterNames()
trace <- singleAnalysis$getTrace(as.integer(3))
traces <- matrix(ncol = length(parameterNames) - 2, nrow = length(trace))
colnames(traces) <- parameterNames[-c(1:2)]
traces[, 1] <- trace
for (i in 4:length(parameterNames)) {
  trace <- singleAnalysis$getTrace(as.integer(i))
  traces[, i - 2] <- trace
}

plot(coda::as.mcmc(traces))
coda::effectiveSize(traces)

