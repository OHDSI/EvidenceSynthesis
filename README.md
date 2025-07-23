EvidenceSynthesis
=================

[![Build Status](https://github.com/OHDSI/EvidenceSynthesis/workflows/R-CMD-check/badge.svg)](https://github.com/OHDSI/EvidenceSynthesis/actions?query=workflow%3AR-CMD-check)
[![codecov.io](https://codecov.io/github/OHDSI/EvidenceSynthesis/coverage.svg?branch=main)](https://app.codecov.io/github/OHDSI/EvidenceSynthesis)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/EvidenceSynthesis)](https://cran.r-project.org/package=EvidenceSynthesis)
[![CRAN_Status_Badge](https://cranlogs.r-pkg.org/badges/EvidenceSynthesis)](https://cran.r-project.org/package=EvidenceSynthesis)

EvidenceSynthesis is part of [HADES](https://ohdsi.github.io/Hades/).

Introduction
============

This R package contains routines for combining causal effect estimates and study diagnostics across multiple data sites in a distributed study. This includes functions for performing meta-analysis and forest plots.

Features
========
- Perform a traditional fixed-effects or random-effects meta-analysis, and create a forest plot.
- Use non-normal approximations of the per-data-site likelihood function to avoid bias when facing small and zero counts.

Example
=======

```r
# Simulate some data for this example:
populations <- simulatePopulations()

# Fit a Cox regression at each data site, and approximate likelihood function:
fitModelInDatabase <- function(population) {
  cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                            data = population,
                                            modelType = "cox")
  cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
  approximation <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "grid with gradients")
  return(approximation)
}
approximations <- lapply(populations, fitModelInDatabase)
approximations <- do.call("rbind", approximations)

# At study coordinating center, perform meta-analysis using per-site approximations:
estimate <- computeBayesianMetaAnalysis(approximations)
estimate
#          mu     mu95Lb   mu95Ub      muSe       tau     tau95Lb   tau95Ub     logRr   seLogRr
# 1 0.5770562 -0.2451619 1.382396 0.4154986 0.2733942 0.004919128 0.7913512 0.5770562 0.4152011
```

Technology
==========
This an R package with some parts implemented in Java.

System requirements
===================
Requires R and Java.

Getting Started
===============
1. Make sure your R environment is properly configured. This means that Java must be installed. See [these instructions](https://ohdsi.github.io/Hades/rSetup.html) for how to configure your R environment.
2. In R, use the following commands to download and install EvidenceSynthesis:

    ```r
    install.packages("EvidenceSynthesis")
    ```
  
User Documentation
==================
Documentation can be found on the [package website](https://ohdsi.github.io/EvidenceSynthesis/).

PDF versions of the documentation are also available:

* Vignette: [Effect estimate using non-normal likelihood approximations](https://raw.githubusercontent.com/OHDSI/EvidenceSynthesis/main/extras/NonNormalEffectSynthesis.pdf)
* Package manual: [EvidenceSynthesis.pdf](https://raw.githubusercontent.com/OHDSI/EvidenceSynthesis/main/extras/EvidenceSynthesis.pdf) 

Support
=======
* Developer questions/comments/feedback: <a href="http://forums.ohdsi.org/c/developers">OHDSI Forum</a>
* We use the <a href="https://github.com/OHDSI/EvidenceSynthesis/issues">GitHub issue tracker</a> for all bugs/issues/enhancements

Contributing
============
Read [here](https://ohdsi.github.io/Hades/contribute.html) how you can contribute to this package.
  
License
=======
EvidenceSynthesis is licensed under Apache License 2.0

Development
===========
This package is being developed in RStudio.

### Development status

Beta
