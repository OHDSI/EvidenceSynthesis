# Copyright 2022 Observational Health Data Sciences and Informatics
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

# Copied from https://github.com/beast-dev/BeastJar/blob/master/R/Utilities.R

#' Determine if Java virtual machine supports Java
#'
#' @description 
#' Tests Java virtual machine (JVM) java.version system property to check if version >= 8.
#'
#' @return
#' Returns TRUE if JVM supports Java >= 8.
#'
#' @examples
#' supportsJava8()
#'
#' @export
supportsJava8 <- function() {
  javaVersionText <-
    rJava::.jcall("java/lang/System", "S", "getProperty", "java.version")
  majorVersion <- as.integer(regmatches(
    javaVersionText,
    regexpr(pattern = "^\\d+", text = javaVersionText)
  ))
  if (majorVersion == 1) {
    twoDigitVersion <- regmatches(javaVersionText,
                                  regexpr(pattern = "^\\d+\\.\\d+", text = javaVersionText))
    majorVersion <- as.integer(regmatches(twoDigitVersion,
                                          regexpr("\\d+$", text = twoDigitVersion)))
  }
  support <- majorVersion >= 8
  # message(paste0("Using JVM version ",
  #                javaVersionText,
  #                " (>= 8? ", support, ")"))
  return(support)
}
