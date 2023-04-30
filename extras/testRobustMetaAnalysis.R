# run this to test if `computeBayesianMetaAnalysis` works
library(EvidenceSynthesis)
data('ncLikelihoods')

res <- computeBayesianMetaAnalysis(ncLikelihoods[[5]],
                                   robust = TRUE,
                                   df = 4)
res
