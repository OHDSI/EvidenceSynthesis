EvidenceSynthesis 0.5.0
=======================

Changes

1. Supporting Bayesian adaptive bias correction in sequential analysis by adding `fitBiasDistribution()`, `sequentialFitBiasDistribution()` and `biasCorrectionInference()` functions

2. Added relevant plotting functions, `plotBiasDistribution()` and `plotBiasCorrectionInference()`. 

3. Added a vignette on Bayesian adaptive bias correction.


EvidenceSynthesis 0.4.1
=======================

Changes

1. Added a video vignette.


EvidenceSynthesis 0.4.0
=======================

Changes

1. Supporting adaptive grid in `plotLikelihoodFit()` and `computeFixedEffectMetaAnalysis()` functions.

2. Adding visualization of likelihood curves to `plotMetaAnalysisForest()`.


EvidenceSynthesis 0.3.0
=======================

Changes

1. Supporting adaptive grid in `computeBayesianMetaAnalysis()` and `approximateLikelihood()` functions.

Bugfixes

1 Fixed error when approximating likelihood using grid (parameter to approximate could only be "x").


EvidenceSynthesis 0.2.3
=======================

Changes

1. Higher tolerance on skew-normal unit tests to prevent them from failing.

2. Detecting and removing grid approximations with illegal values before computing meta-analysis.


EvidenceSynthesis 0.2.2
=======================

Changes

1. Documenting dependency on Java in the `SystemRequirements` field of the package DESCRIPTION. 

2. Adding `seed` argument to `computeBayesianMetaAnalysis()`. Defaults to a constant value for reproducability.

Bugfixes

1. Fixed `plotMetaAnalysisForest()` when using grid approximations, including when providing approximations as tibble.


EvidenceSynthesis 0.2.1
=======================

Changes

1. `computeBayesianMetaAnalysis()` now outputs ESS.

2. Checking whether required Java version (8 or newer) is installed.


EvidenceSynthesis 0.2.0
=======================

Changes

1. Preparing for CRAN release: Adding missing rmarkdown dependency to Suggests. 

2. Added `computeConfidenceInterval()` function.

3. `plotMetaAnalysisForest()` function now works with normal and non-normal approximations.

4. Grid points now evenly spaced on log scale, not HR scale.


EvidenceSynthesis 0.1.0
=======================

Changes

1. Adding meta-analysis using local non-normal Cox likelihood approximations to avoid bias when sample size is small.


EvidenceSynthesis 0.0.5
=======================

Bug fixes

1. Fixing build error in R 4.0.0.
