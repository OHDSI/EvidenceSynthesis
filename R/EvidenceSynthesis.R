# Copyright 2021 Observational Health Data Sciences and Informatics
#
# This file is part of EvidenceSynthesis
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' @keywords internal
"_PACKAGE"

#' @importFrom grDevices rgb
#' @importFrom stats density dnorm qnorm quantile runif coef confint median nlm optim pnorm
#' printCoefmat qchisq rexp rnorm
#' @importFrom rlang .data abort warn inform
#' @importFrom methods is
#' @import BeastJar
#' @import survival
#'
NULL

.onLoad <- function(libname, pkgname) {
  beastLocation <- system.file("java/beast.jar", package = "BeastJar")
  rJava::.jpackage(pkgname, lib.loc = libname, morePaths = beastLocation)
}
