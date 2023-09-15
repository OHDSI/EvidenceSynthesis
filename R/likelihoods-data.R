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
#' ncLikEx <- ncLikelihoods[["5"]][[1]]
#' \donttest{
#' plot(value ~ point, data = ncLikEx)
#' }
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
#' ooiLikEx <- ooiLikelihoods[["5"]][[1]]
#' \donttest{
#' plot(value ~ point, data = ooiLikEx)
#' }
"ooiLikelihoods"


#' Example profile likelihoods for hierarchical meta analysis with bias correction
#'
#' A list that contains profile likelihoods for two negative control outcomes and
#' a synthetic outcome of interest, across four data sources.
#' Each element of the list contains profile likelihoods for one outcome,
#' where each row provides profile likelihood values (over a grid) from one data source.
#'
#' @docType data
#'
#' @format An objects of class `list`; the list contains 3 dataframes,
#' where each dataframe includes four rows of likelihood function values
#' corresponding to the points in the column names.
#'
#' @keywords datasets
#'
#' @examples
#' data("hmaLikelihoodList")
#' hmaLikEx <- hmaLikelihoodList[[1]]
#' \donttest{
#' plot(as.numeric(hmaLikEx[2,]) ~ as.numeric(names(hmaLikEx)))
#' }
"hmaLikelihoodList"



#' A bigger example of profile likelihoods for hierarchical meta analysis with bias correction
#'
#' A list that contains profile likelihoods for 10 negative control outcomes and
#' an outcome of interest, across data sources.
#' Each element of the list contains a named list of profile likelihoods for one outcome,
#' where each element is a data frame that provides likelihood values over a grid of
#' parameter values, the element name corresponding to data source name.
#'
#' @docType data
#'
#' @format An objects of class `list`; the list contains 11 named lists, each list for one outcome.
#' Each list contains data frames that record profile likelihoods from different data sources.
#' The first 10 list corresponds to 10 negative control outcomes,
#' whereas the last list the outcome of interest.
#'
#' @keywords datasets
#'
#' @examples
#' data("likelihoodProfileLists")
#' exLP <- likelihoodProfileLists[[1]][[1]]
#' \donttest{
#' plot(value ~ point, data = exLP)
#' }
"likelihoodProfileLists"
