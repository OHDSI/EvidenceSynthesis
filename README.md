EvidenceSynthesis
=================

[![Build Status](https://travis-ci.org/OHDSI/EvidenceSynthesis.svg?branch=master)](https://travis-ci.org/OHDSI/EvidenceSynthesis)
[![codecov.io](https://codecov.io/github/OHDSI/EvidenceSynthesis/coverage.svg?branch=master)](https://codecov.io/github/OHDSI/EvidenceSynthesis?branch=master)

EvidenceSynthesis is part of [HADES](https://ohdsi.github.io/Hades/).

Introduction
============

This R package contains routines for combining causal effect estimates and diagnostics across multiple data sites in a distributed study. This includes functions for performing meta-analysis and forest plots.

Features
========
- Perform a meta-analysis
- Create a forest plot

Screenshots and examples
========================
to do

Technology
==========
This an R package with some parts implemented in Java.

System requirements
===================
Requires R and Java.

Getting Started
===============
1. Make sure your R environment is properly configured. This means that Java must be installed. See [these instructions](https://ohdsi.github.io/MethodsLibrary/rSetup.html) for how to configure your R environment.
2. In R, use the following commands to download and install EvidenceSynthesis:

    ```r
    install.packages("remotes")
    library(remotes)
    install_github("ohdsi/EvidenceSynthesis")
    ```
  
User Documentation
==================
Documentation can be found on the [package website](https://ohdsi.github.io/EvidenceSynthesis/).

PDF versions of the documentation are also available:
* Vignette: [Effect estimate using non-normal likelihood approximations](https://raw.githubusercontent.com/OHDSI/CaseControl/master/inst/doc/SingleStudies.pdf)
* Package manual: [EvidenceSynthesis.pdf](https://raw.githubusercontent.com/OHDSI/EvidenceSynthesis/master/extras/EvidenceSynthesis.pdf) 

Support
=======
* Developer questions/comments/feedback: <a href="http://forums.ohdsi.org/c/developers">OHDSI Forum</a>
* We use the <a href="../../issues">GitHub issue tracker</a> for all bugs/issues/enhancements

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

Under development. Do not use.
