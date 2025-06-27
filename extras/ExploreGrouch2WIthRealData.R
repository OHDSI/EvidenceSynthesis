library(DatabaseConnector)
library(dplyr)

connectionDetails <- createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("assureServer"), keyring::key_get("assureDatabase"), sep = "/"),
  user = keyring::key_get("assureUser"),
  password = keyring::key_get("assurePassword")
)

profilesDatabaseSchema <- "dp_temp"
profilesTable <- "db_profile_results"
vocabularyDatabaseSchema <- "vocabulary_20240229"

# Fetch data -------------------------------------------------------------------
# recordCountAnalysisIds <- c(401, 601, 701, 801, 901, 1801, 2101)
recordCountAnalysisIds <- c(401) # Just conditions for now
numberOfPersonsAnalysisId <- 1
sexAnalysisId <- 2
ageAnalysisId <- 101
lengthOfOpAnalysisId <- 108
numberOfOpAnalysisId <- 113

connection <- connect(connectionDetails)
sql <- "
SELECT cdm_source_name,
  analysis_id,
  CAST(stratum_1 AS INT) AS stratum_id,
  count_value
FROM @profiles_database_schema.@profiles_table
WHERE analysis_id IN (@analysis_ids);
"
demographics <- renderTranslateQuerySql(
  connection = connection,
  sql = sql,
  profiles_database_schema = profilesDatabaseSchema,
  profiles_table = profilesTable,
  analysis_ids = c(numberOfPersonsAnalysisId,
                   sexAnalysisId,
                   ageAnalysisId,
                   lengthOfOpAnalysisId,
                   numberOfOpAnalysisId),
  snakeCaseToCamelCase = TRUE
)
sql <- "
SELECT cdm_source_name,
  ancestor_concept_id AS concept_id,
  concept_name,
  SUM(count_value) AS record_count
FROM @profiles_database_schema.@profiles_table
INNER JOIN @vocabulary_database_schema.concept_ancestor
  ON CAST(stratum_1 AS INT) = descendant_concept_id
INNER JOIN (
	SELECT concept_id,
	  concept_name
	FROM @vocabulary_database_schema.concept
	INNER JOIN (
	  SELECT *
	  FROM @vocabulary_database_schema.concept_ancestor
	  WHERE ancestor_concept_id = 441840 /* SNOMED clinical finding */
	  AND (min_levels_of_separation > 2
		OR descendant_concept_id IN (433736, 433595, 441408, 72404, 192671, 137977, 434621, 437312, 439847, 4171917, 438555, 4299449, 375258, 76784, 40483532, 4145627, 434157, 433778, 258449, 313878)
		)
	) temp
	  ON concept_id = descendant_concept_id
	WHERE concept_name NOT LIKE '%finding'
		AND concept_name NOT LIKE 'Disorder of%'
		AND concept_name NOT LIKE 'Finding of%'
		AND concept_name NOT LIKE 'Disease of%'
		AND concept_name NOT LIKE 'Injury of%'
		AND concept_name NOT LIKE '%by site'
		AND concept_name NOT LIKE '%by body site'
		AND concept_name NOT LIKE '%by mechanism'
		AND concept_name NOT LIKE '%of body region'
		AND concept_name NOT LIKE '%of anatomical site'
		AND concept_name NOT LIKE '%of specific body structure%'
		AND domain_id = 'Condition'
) valid_groups
	ON ancestor_concept_id = valid_groups.concept_id
WHERE analysis_id IN (@analysis_ids)
GROUP BY cdm_source_name,
  ancestor_concept_id,
  concept_name
"
conceptCounts <- renderTranslateQuerySql(
  connection = connection,
  sql = sql,
  profiles_database_schema = profilesDatabaseSchema,
  profiles_table = profilesTable,
  vocabulary_database_schema = vocabularyDatabaseSchema,
  analysis_ids = c(recordCountAnalysisIds),
  snakeCaseToCamelCase = TRUE
)
# executeSql(connection, "COMMIT;")

disconnect(connection)

conceptSummary <- conceptCounts |>
  group_by(conceptId, conceptName) |>
  summarise(maxCount = max(recordCount),
            databaseCount = n(),
            .groups = "drop") |>
  arrange(desc(maxCount))

saveRDS(demographics, "e:/temp/demographics.rds")
saveRDS(conceptCounts, "e:/temp/conceptCounts.rds")

# Perform analysis on a single concept -----------------------------------------
conceptId <- 436665 # Bipolar disorder

demographics <- readRDS("e:/temp/demographics.rds")
conceptCounts <- readRDS("e:/temp/conceptCounts.rds")

conceptCounts <- conceptCounts |>
  filter(conceptId == !!conceptId)

nInternalKnots <- 1
splineDegree <- 2
nAge <- demographics |>
  filter(analysisId == ageAnalysisId) |>
  select(age = stratumId, count= countValue) |>
  group_by(age) |>
  summarise(count = sum(count)) |>
  arrange(age) |>
  mutate(cumCount = cumsum(count))
nTotal <- sum(nAge$count)
knotsAtCumCount <- (1:nInternalKnots / (nInternalKnots+1)) * nTotal
internalKnots <- sapply(knotsAtCumCount, function(x) nAge$age[max(which(nAge$cumCount <= x))])
boundaryKnots <- c(min(nAge$age), max(nAge$age))
designMatrix <- splines::bs(x = boundaryKnots[1]:boundaryKnots[2],
                            knots = internalKnots,
                            Boundary.knots = boundaryKnots,
                            degree = splineDegree)

# cdmSourceName = "merative_mdcr_v3466"
prepareDatabaseData <- function(cdmSourceName) {
  nPersons <- demographics |>
    filter(cdmSourceName == !!cdmSourceName, analysisId == numberOfPersonsAnalysisId) |>
    pull(countValue)

  # Length of observation is not in the model, and since we assume it is independent from age and sex
  # we just need the total days observed:
  meanDaysPerObservationPeriod <- demographics |>
    filter(cdmSourceName == !!cdmSourceName, analysisId == lengthOfOpAnalysisId) |>
    mutate(daysObserved = (15 + stratumId * 30) * countValue) |>
    summarise(daysObserved = sum(daysObserved) / sum(countValue)) |>
    pull()
  nObservationPeriod <- demographics |>
    filter(cdmSourceName == !!cdmSourceName, analysisId == numberOfOpAnalysisId)  |>
    mutate(nObservationPeriod = stratumId * countValue) |>
    summarise(nObservationPeriod = sum(nObservationPeriod)) |>
    pull()
  daysObserved <- meanDaysPerObservationPeriod * nObservationPeriod

  # Proportion females?
  pFemale <- demographics |>
    filter(cdmSourceName == !!cdmSourceName, analysisId == sexAnalysisId) |>
    filter(stratumId == 8532) |>
    mutate(proportion = countValue / nPersons) |>
    pull(proportion)

  # Age variables
  nAge <- demographics |>
    filter(cdmSourceName == !!cdmSourceName, analysisId == ageAnalysisId) |>
    select(age = stratumId, count = countValue) |>
    right_join(tibble(age = boundaryKnots[1]:boundaryKnots[2]), by = join_by(age)) |>
    mutate(count = if_else(is.na(count), 0, count)) |>
    arrange(age)
  ageVars <- apply(nAge, 1, function(x) designMatrix[x[1] - boundaryKnots[1] + 1, ] * x[2] / nPersons, simplify = FALSE)
  ageVars <- do.call(rbind, ageVars) |>
    apply(2, sum) |>
    t() |>
    as_tibble()
  # ageVars <- ageVars[, -1]
  colnames(ageVars) <- paste0("age", colnames(ageVars))

  # Concept record count
  nEvents <- conceptCounts |>
    filter(cdmSourceName == !!cdmSourceName) |>
    pull(recordCount)

  row <- bind_cols(
    ageVars,
    tibble(
      female = pFemale,
      daysObserved = daysObserved,
      cdmSourceName = cdmSourceName,
      nEvents = nEvents
    )
  )
  return(row)
}
data <- lapply(conceptCounts$cdmSourceName, prepareDatabaseData)
data <- bind_rows(data)
ageVarNames <- colnames(data)[grep("^age", colnames(data))]
formula <- as.formula(sprintf("nEvents ~ %s + female + offset(log(daysObserved))",
                              paste(ageVarNames, collapse = " + ")))


plot(data$female, data$nEvents / data$daysObserved)
plot(data$age1, data$nEvents / data$daysObserved)
plot(data$age2, data$nEvents / data$daysObserved)
plot(data$age3, data$nEvents / data$daysObserved)

# Create prediction intervals using bootstrap ----------------------------------
# Note: Unable to fit a model when using 14 internal databases

# data$rate <- 100 * 365.25 * data$nEvents / data$daysObserved

doPrediction <- function(dummy, data) {
  print(dummy)
  sampledData <- data[sample.int(nrow(data), nrow(data), replace = TRUE), ]
  # sampledData <- data
  cyclopsData <- Cyclops::createCyclopsData(formula, modelType = "pr", data = sampledData)
  fit <- Cyclops::fitCyclopsModel(cyclopsData, prior =  Cyclops::createPrior("laplace", 0.01))
  if (fit$return_flag != "SUCCESS") {
    return()
  }
  # Cyclops predict() does not support predicting for new data, except when using sparse representation,
  # so using own implementation instead:
  coefs <- coef(fit)
  prediction <- exp(coefs[1] + apply(t((t(data[ageVarNames]) * coefs[ageVarNames])), 1, sum) + coefs["female"]*data$female) * data$daysObserved
  # data$prediction = prediction
  return(prediction)
}

bootstrap <- lapply(1:100, doPrediction, data = data)
bootstrap <- do.call(rbind, bootstrap)

alpha <- 0.01 / nrow(data)
cis <- apply(bootstrap, 2, function(x) quantile(x, c(0.5, alpha/2, 1 - alpha/2)))
rownames(cis) <- c("predictedEvents", "lb", "ub")
dataWithPredictions <- data |>
  bind_cols(as_tibble(t(cis)))

dataWithPredictions <- dataWithPredictions |>
  mutate(outlier = nEvents < lb | nEvents > ub)




# # Some code to understand analysis IDs:
# analysisIds <- renderTranslateQuerySql(
#   connection = connection,
#   sql = "SELECT analysis_id, COUNT(*) FROM @schema.@table GROUP BY analysis_id;",
#   schema = schema,
#   table = table,
#   snakeCaseToCamelCase = TRUE
# )
# analyses <- readr::read_csv("https://raw.githubusercontent.com/OHDSI/Achilles/refs/heads/main/inst/csv/achilles/achilles_analysis_details.csv") |>
#   SqlRender::snakeCaseToCamelCaseNames()
# analyses <- analyses |>
#   inner_join(analysisIds, by = join_by(analysisId))
# subset <- analyses[grepl("records", analyses$analysisName) & !grepl("source", analyses$analysisName), ]
