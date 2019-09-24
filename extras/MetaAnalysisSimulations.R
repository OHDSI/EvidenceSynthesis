library(EvidenceSynthesis)
settings <- createSimulationSettings(randomEffectSd = 1)
simulation <- runSimulation(settings)


x <- evaluateMetaAnalysisApproach(settings, 100, estimateByPooling, useCyclops = FALSE)

x <- evaluateMetaAnalysisApproach(settings, 100, estimateUsingStandardMetaAnalysis, assumeRandomEffect = TRUE)


# Low count simulation -------------------------------

settings <- createSimulationSettings(nSites = 5, n = 5000, nStrata = 5, hazardRatio = 0.5)
set.seed(6)
simulation <- runSimulation(settings)
summary(simulation)

estimateByPooling(simulation, useCyclops = TRUE)
estimateByPooling(simulation, useCyclops = FALSE)
estimateUsingStandardMetaAnalysis(simulation, assumeRandomEffect = FALSE)
estimateUsingStandardMetaAnalysis(simulation, assumeRandomEffect = TRUE)
plotLikelihood(simulation)

pool <- evaluateMetaAnalysisApproach(simulationSettings = settings,
                                     iterations = 100,
                                     metaAnalysisFunction = estimateByPooling)

                                     
metaAnalysis <- evaluateMetaAnalysisApproach(simulationSettings = settings,
                                     iterations = 100,
                                     metaAnalysisFunction = estimateUsingStandardMetaAnalysis,
                                     assumeRandomEffect = FALSE)


