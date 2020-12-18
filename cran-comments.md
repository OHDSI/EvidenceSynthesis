We are confident these changes address the issue in the ATLAS test environment, which we think were caused by a different random number generator in that environment. As a consequence, the unit tests take slightly longer to execute. (46 seconds on my laptop). Note that there does not appear to be an ATLAS equivalent environment on rhub, which would be helpful to debug.

We are happy to see Java has now been installed on CRAN's 'r-oldrel-windows-ix86+x86_64' testing environment, resolving the issue noted by professor Ripley. For future reference, we have explicitly documented the dependency on Java in the SystemRequirements field of the DESCRIPTION. 

Includes 2 changes and 1 bugfix (see NEWS.md)

---

## Test environments
* Ubuntu 16.04.6 LTS (Travis), R 4.0.2
* Windows 10, R 4.0.3

## R CMD check results

There were no ERRORs or WARNINGs. 

## Downstream dependencies

There are no downstream dependencies.