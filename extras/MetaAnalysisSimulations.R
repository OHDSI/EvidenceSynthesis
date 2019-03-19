settings <- createSimulationSettings(randomEffectSd = 1)
simulation <- runSimulation(settings)


x <- evaluateMetaAnalysisApproach(settings, 100, estimateByPooling, useCyclops = FALSE)

x <- evaluateMetaAnalysisApproach(settings, 100, estimateUsingStandardMetaAnalysis, assumeRandomEffect = TRUE)


# Low count simulation -------------------------------

settings <- createSimulationSettings(nSites = 10, n = 5000, nStrata = 5)
set.seed(6)
simulation <- runSimulation(settings)
summary(simulation)

estimateByPooling(simulation, useCyclops = TRUE)
estimateByPooling(simulation, useCyclops = FALSE)
estimateUsingStandardMetaAnalysis(simulation, assumeRandomEffect = FALSE)
estimateUsingStandardMetaAnalysis(simulation, assumeRandomEffect = TRUE)
