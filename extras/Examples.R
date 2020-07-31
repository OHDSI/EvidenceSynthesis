

# Assume each site has approximated their local likelihood using the approximateLikelihood function
approximations <- data.frame(sites = c("A", "B", "C", "D"),
                             mu = c(0, 0.1, -0.1, 0),
                             sigma = c(0.1, 0.1, 0.1, 0.1),
                             gamma = c(0, 0, 0, 0.1))

estimate <- computeBayesianMetaAnalysis(approximations)
estimate
# mu     mu95Lb    mu95Ub     muSe        tau      tau95Lb  tau95Ub
# 0.0003129562 -0.1747429 0.1723472 0.089661 0.07759992 0.0002024991 0.301007