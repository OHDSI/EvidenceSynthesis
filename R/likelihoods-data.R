# Copyright 2023 Observational Health Data Sciences and Informatics
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

#' Example profile likelihoods for negative control outcomes
#'
#' A list that contain profile likelihoods a large set of negative control outcomes.
#' They are extracted from a real-world observational healthcare database, with the
#' likelihoods profiled using adaptive grids using the `Cyclops` package.
#'
#' @docType data
#'
#' @format An object of class `list` containing 12 lists, where each list includes
#' several dataframes ith column `point` and `value` for adaptive grid profile likelihoods.
#'
#' @keywords datasets
#'
#' @references Schuemie et al. (2022). Vaccine safety surveillance using routinely collected
#' healthcare data—an empirical evaluation of epidemiological designs. Frontiers in Pharmacology.
#'
#' @examples
#' data("ncLikelihoods")
#' ncLikEx = ncLikelihoods[["5"]][[1]]
#' \donttest{plot(value ~ point, data = ncLikEx)}
"ncLikelihoods"


#' Example profile likelihoods for a synthetic outcome of interest
#'
#' A list that contain profile likelihoods for a synthetic outcome of interest.
#' They are extracted from a real-world observational healthcare database, with the
#' likelihoods profiled using adaptive grids using the `Cyclops` package.
#'
#' @docType data
#'
#' @format An objects of class `list`; the list contains 12 lists,
#' where each list includes several dataframes with column `point` and `value`
#' for adaptive grid profile likelihoods.
#'
#' @keywords datasets
#'
#' @references Schuemie et al. (2022). Vaccine safety surveillance using routinely collected
#' healthcare data—an empirical evaluation of epidemiological designs. Frontiers in Pharmacology.
#'
#' @examples
#' data("ooiLikelihoods")
#' ooiLikEx = ooiLikelihoods[["5"]][[1]]
#' \donttest{plot(value ~ point, data = ooiLikEx)}
"ooiLikelihoods"



