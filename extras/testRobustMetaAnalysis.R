# run this to test if `computeBayesianMetaAnalysis` works
data('ncLikelihoods')

res = computeBayesianMetaAnalysis(ncLikelihoods[[5]], robust = TRUE)
