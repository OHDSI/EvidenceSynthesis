# @file PackageMaintenance
#
# Copyright 2017 Observational Health Data Sciences and Informatics
#
# This file is part of EmpiricalCalibration
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

# Format and check code:
OhdsiRTools::formatRFolder()
OhdsiRTools::checkUsagePackage("EmpiricalCalibration")
OhdsiRTools::updateCopyrightYearFolder()

# Create manual and vignettes:
shell("rm extras/EmpiricalCalibration.pdf")
shell("R CMD Rd2pdf ./ --output=extras/EmpiricalCalibration.pdf")

rmarkdown::render("vignettes/EmpiricalPCalibrationVignette.Rmd",
                  output_file = "../inst/doc/EmpiricalPCalibrationVignette.pdf",
                  rmarkdown::pdf_document(latex_engine = "pdflatex",
                                          toc = TRUE,
                                          number_sections = TRUE))

rmarkdown::render("vignettes/EmpiricalCICalibrationVignette.Rmd",
                  output_file = "../inst/doc/EmpiricalCiCalibrationVignette.pdf",
                  rmarkdown::pdf_document(latex_engine = "pdflatex",
                                          toc = TRUE,
                                          number_sections = TRUE))

# Release package:
devtools::build_win()

devtools::release()
