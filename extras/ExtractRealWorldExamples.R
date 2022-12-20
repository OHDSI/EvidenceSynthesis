# Find example in depression study that had most assymetric confidence intervals, while still having estimates for all 4 DBs
# Extracts data from the S3 archive.
dbs <- data.frame(s3Folders = c("mschuemi/Epi421/PopEstDepression_Ccae/",
                                "mschuemi/Epi421/PopEstDepression_Mdcd/",
                                "mschuemi/Epi421/PopEstDepression_Mdcr/",
                                "mschuemi/Epi421/PopEstDepression_Optum/"),
                  database = c("CCAE",
                               "MDCD",
                               "MDCR",
                               "Optum"),
                  stringsAsFactors = FALSE)
localFolder <- "d:/SmallCountMetaAnalysis/Example"

# Fetch full set of estimates from S3 storage -----------------------------------------------------------
getEstimates <- function(database) {
  writeLines(paste("Getting file for", database))
  s3Folder <- dbs$s3Folders[dbs$database == database]
  ArchiveR::restoreFile(paste0(s3Folder, "calibratedEstimates.csv"), localFolder = localFolder)
  estimates <- read.csv(file.path(localFolder, "calibratedEstimates.csv"))
  estimates$database <- database
  unlink(file.path(localFolder, "calibratedEstimates.csv"))
  return(estimates)
}

estimates <- lapply(dbs$database, getEstimates)
estimates <- do.call(rbind, estimates)

exposures <- data.frame(
  cohortId = c(750982, 725131, 4030840, 4327941, 703547, 717607, 715259, 743670, 797617, 715939, 755695, 722031, 739138, 40234834, 710062, 738156, 721724),
  cohortName = c("Bupropion", "Mirtazapine", "Electroconvulsive therapy", "Psychotherapy", "Trazodone", "Desvenlafaxine", "duloxetine", "venlafaxine", "Citalopram", "Escitalopram", "Fluoxetine", "Paroxetine", "Sertraline", "vilazodone", "Amitriptyline", "Doxepin", "Nortriptyline")
)
outcomes <- data.frame(
  outcomeId = c(2557, 2556, 2551, 2559, 2820, 2821, 2822, 2823, 2824, 2826, 2828, 2829, 2830, 2831, 2832, 2833, 2834, 2835, 2836, 2837, 2838, 2839),
  outcomeName = c("Acute myocardial infarction", "Gastrointestinal hemhorrage", "Open-angle glaucoma", "Stroke", "Suicide and suicidal ideation", "Insomnia", "Diarrhea", "Nausea", "Hypothyroidism", "Constipation", "Seizure", "Delirium", "Alopecia", "Tinnitus", "Vertigo", "Hyponatremia", "Decreased libido", "Fracture", "Hypotension", "Acute liver injury", "Hyperprolactinemia", "Ventricular arrhythmia and sudden cardiac death")
)

library(dplyr)
estimates <- estimates %>%
  mutate(targetId = floor(targetId / 1000),
         comparatorId = floor(comparatorId / 1000)) %>%
  inner_join(
    exposures %>%
      select(targetId = "cohortId",
             targetName = "cohortName"),
    by = "targetId"
  ) %>%
  inner_join(
    exposures %>%
      select(comparatorId = "cohortId",
             comparatorName = "cohortName"),
    by = "comparatorId"
  ) %>%
  inner_join(
    outcomes,
    by = "outcomeId"
  )
saveRDS(estimates, file.path(localFolder, "allEstimates.rds"))

# Some more preprocessing --------------------------------------------------

# Exclude psychotherapy and ECT
exclude <- c("Psychotherapy", "Electroconvulsive therapy")

results <- readRDS(file.path(localFolder, "allEstimates.rds"))
results <- results[results$analysisId == 3, ] # On-treatment, matching
results <- results[!(results$targetName %in% exclude) & !(results$comparatorName %in% exclude), ]
# results <- results[results$outcomeType == "hoi", ]
results$ratio <- 1 - (log(results$rr) - log(results$ci95lb)) / (log(results$ci95ub) - log(results$rr))
results$absRatio <- abs(results$ratio)
results$validCount <- 0
results$validCount[!is.na(results$seLogRr)] <- 1
validCount <- aggregate(validCount ~ targetId + comparatorId + outcomeId, results, sum)
results$anyCount <- 1
anyCount <- aggregate(anyCount ~ targetId + comparatorId + outcomeId, results, sum)
results$targetZero <- results$eventsTreated == 0
targetZeroCount <- aggregate(targetZero ~ targetId + comparatorId + outcomeId, results, sum)
results$comparatorCount <- results$eventsComparator
comparatorCount <- aggregate(comparatorCount ~ targetId + comparatorId + outcomeId, results, mean)

examples <- data.frame()
# Example 1: Estimates in all DBs, but maximum mean skew
meanRatio <- aggregate(absRatio ~ targetId + comparatorId + outcomeId, results, mean, na.action = na.exclude)
meanRatio <- merge(meanRatio, validCount)
meanRatio <- meanRatio[meanRatio$validCount == 4, ]
meanRatio <- meanRatio[order(-meanRatio$absRatio), ]

index <- 1
example1 <- results[results$targetId == meanRatio$targetId[index] & results$comparatorId == meanRatio$comparatorId[index] & results$outcomeId == meanRatio$outcomeId[index], ]
example1
example1$exampleId <- 1
examples <- rbind(examples, example1)

# Example 2: Estimates in two DBs
meanRatio <- aggregate(absRatio ~ targetId + comparatorId + outcomeId, results, mean, na.action = na.exclude)
meanRatio <- merge(meanRatio, validCount)
meanRatio <- merge(meanRatio, anyCount)
meanRatio <- meanRatio[meanRatio$validCount == 2, ]
meanRatio <- meanRatio[meanRatio$anyCount == 4, ]
meanRatio <- meanRatio[order(-meanRatio$absRatio), ]

index <- 1
example2 <- results[results$targetId == meanRatio$targetId[index] & results$comparatorId == meanRatio$comparatorId[index] & results$outcomeId == meanRatio$outcomeId[index], ]
example2
example2$exampleId <- 2
examples <- rbind(examples, example2)

# Example 3: Estimates in two DBs, missing estimates have 0 count in target
meanRatio <- aggregate(absRatio ~ targetId + comparatorId + outcomeId, results, mean, na.action = na.exclude)
meanRatio <- merge(meanRatio, validCount)
meanRatio <- merge(meanRatio, anyCount)
meanRatio <- merge(meanRatio, targetZeroCount)
meanRatio <- merge(meanRatio, comparatorCount)
meanRatio <- meanRatio[meanRatio$validCount == 2, ]
meanRatio <- meanRatio[meanRatio$anyCount == 4, ]
meanRatio <- meanRatio[meanRatio$targetZero == 2, ]
meanRatio <- meanRatio[order(-meanRatio$comparatorCount), ]

index <- 1
example3 <- results[results$targetId == meanRatio$targetId[index] & results$comparatorId == meanRatio$comparatorId[index] & results$outcomeId == meanRatio$outcomeId[index], ]
example3
example3$exampleId <- 3
examples <- rbind(examples, example3)

write.csv(examples, file.path(localFolder, "Examples.csv"), row.names = FALSE)

# sum(is.na(results$calSeLogRr))
# temp <- example$targetId
# example$targetId <- example$comparatorId
# example$comparatorId <- temp
#
# example <- example2


# Fetch patient-level data from S3 storage ----------------------------------
examples <- read.csv(file.path(localFolder, "Examples.csv"))
# ArchiveR::listFiles("mschuemi/Epi421")

for (database in dbs$database) {
  writeLines(paste("Getting files for", database))
  s3Folder <- dbs$s3Folders[dbs$database == database]
  ArchiveR::restoreFile(paste0(s3Folder, "cmOutput/", "outcomeModelReference.rds"), localFolder = localFolder)
  omr <- readRDS(file.path(localFolder, "outcomeModelReference.rds"))
  unlink(file.path(localFolder, "outcomeModelReference.rds"))
  for (exampleId in unique(examples$exampleId)) {
    writeLines(paste("Getting file for example ", exampleId))
    example <- examples[examples$exampleId == exampleId, ]
    strataFile <- omr$strataFile[floor(omr$targetId/1000) == example$targetId[1] &
                                   floor(omr$comparatorId/1000) == example$comparatorId[1] &
                                   omr$outcomeId == example$outcomeId[1] &
                                   omr$analysisId == 3]
    switchTandC <- FALSE
    if (length(strataFile) == 0) {
      temp <- example$targetId
      example$targetId <- example$comparatorId
      example$comparatorId <- temp
      temp <- example$targetName
      example$targetName <- example$comparatorName
      example$comparatorName <- temp
      strataFile <- omr$strataFile[floor(omr$targetId/1000) == example$targetId[1] &
                                     floor(omr$comparatorId/1000) == example$comparatorId[1] &
                                     omr$outcomeId == example$outcomeId[1] &
                                     omr$analysisId == 3]
      switchTandC <- TRUE

    }
    strataFile <- basename(strataFile)
    localFile <- file.path(localFolder, sprintf("StratPop_%s_%s.rds", database, exampleId))
    ArchiveR::restoreFile(paste0(s3Folder, "cmOutput/", strataFile), localFile = localFile)

    # Cleanup
    x <- readRDS(localFile)
    x <- data.frame(rowId = x$rowId,
                    stratumId = x$stratumId,
                    x = x$treatment,
                    y = as.integer(x$outcomeCount != 0),
                    time = x$survivalTime)
    if (switchTandC) {
      x$x <- 1 - x$x
    }
    saveRDS(x, localFile)
  }
}

