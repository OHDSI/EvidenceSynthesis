#' Example profile likelihoods
#'
#' Two lists, `ooiLikelihoods` and `ncLikelihoods`, that contain profile likelihoods
#' for a synthetic outcome of interest and a large set of negative control outcomes.
#' They are extracted from a real-world observational healthcare database, with the
#' likelihoods profiled using adaptive grids using the `Cyclops` package.
#'
#' @docType data
#'
#' @format Two objects of class `list`; each list contains 12 lists,
#' where each list includes several dataframes with column `point` and `value`
#' for adaptive grid profile likelihoods.
#'
#' @keywords datasets
#'
#' @references Schuemie et al. (2022). Vaccine safety surveillance using routinely collected
#' healthcare data—an empirical evaluation of epidemiological designs. Frontiers in Pharmacology.
#'
#' @examples
#' data("likelihoods")
#' ncLikEx = ncLikelihoods[["5"]][[1]]
#' \donttest{plot(value ~ point, data = ncLikEx)}
"likelihoods"


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



