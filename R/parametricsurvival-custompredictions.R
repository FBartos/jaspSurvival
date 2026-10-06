# Custom predictions share selected fits with the existing prediction output.
# Source validation and fitting-data defaults live in custompredictions-input.R.
.sapCustomPredictionDependencies <- c(
  "customPredictions", "customPredictionSource", "customPredictionData", "customPredictionFile",
  "customPredictionMean", "customPredictionPercentiles", "customPredictionPercentileValues",
  "customPredictionConfidenceInterval", "customPredictionConfidenceLevels",
  "customPredictionTable", "customPredictionExport", "customPredictionColumnPrefix"
)

.sapCustomPredictions <- function(jaspResults, original, fitting, options) {

  export <- options[["customPredictionSource"]] == "dataset" && options[["customPredictionExport"]]
  if (!options[["customPredictionTable"]] && !export)
    return()
  if (!is.null(jaspResults[["customPredictionsResult"]]))
    return()

  dependencies <- c(.sapGetDependencies(options), .sapSelectedOutputDependencies,
                    .sapSimulationDependencies, .sapCustomPredictionDependencies)
  container <- createJaspContainer(title = gettext("Custom Predictions"))
  container$dependOn(dependencies)
  container$position <- 3.5
  jaspResults[["customPredictionsResult"]] <- container
  if (!.saSurvivalReady(options)) {
    .sapCustomPredictionPlaceholder(container, options)
    return()
  }

  fits <- .sapFlattenFit(.sapExtractFit(jaspResults, options, type = "selected"), options)
  required <- unique(c(unlist(lapply(fits, function(fit) attr(fit, "modelTerms")[["components"]]), use.names = FALSE),
                       unlist(options[["subgroup"]], use.names = FALSE)))
  setup <- try(list(
    specification = .sapCustomPredictionSpecification(options),
    input         = .sapCustomPredictionInput(original, fitting, options, required)
  ), silent = TRUE)
  if (inherits(setup, "try-error")) {
    .sapCustomPredictionError(container, "inputError", gettext("Predictions"), setup)
    return()
  }

  specification <- setup[["specification"]]
  records <- lapply(fits, function(fit)
    try(.sapCustomPredictionRecord(setup[["input"]], original, fit, specification, options), silent = TRUE))
  messages <- character(0)
  for (i in seq_along(records)) {
    record <- records[[i]]
    key    <- paste0("model", i)
    if (inherits(record, "try-error")) {
      .sapCustomPredictionError(container, key, attr(fits[[i]], "label"), record)
    } else if (options[["customPredictionTable"]]) {
      container[[key]] <- .sapCustomPredictionTable(record, specification)
    } else {
      notes <- .sapCustomPredictionMessages(record, specification)
      label <- attr(fits[[i]], "label")
      if (length(notes) > 0L)
        messages <- c(messages, paste0(if (nzchar(label)) paste0(label, ": ") else "", notes))
    }
  }
  if (length(messages) > 0L) {
    table <- createJaspTable(title = gettext("Prediction Notes"))
    container[["messages"]] <- table
    .saFillMessageTable(table, unique(messages))
  }
  if (export)
    .sapCustomPredictionExport(container, records, specification, original, options, dependencies)

  return()
}

.sapCustomPredictionPlaceholder <- function(container, options) {

  table <- createJaspTable(title = gettext("Predictions"))
  table$addColumnInfo(name = "row", title = gettext("Observation"), type = "integer")
  container[["predictions"]] <- table
  specification <- try(.sapCustomPredictionSpecification(options), silent = TRUE)
  if (inherits(specification, "try-error"))
    table$setError(conditionMessage(attr(specification, "condition")))
  else
    .sapCustomPredictionTableColumns(table, specification)

  return()
}

.sapCustomPredictionError <- function(container, key, title, error) {

  table <- createJaspTable(title = title)
  container[[key]] <- table
  table$setError(conditionMessage(attr(error, "condition")))

  return()
}

.sapCustomPredictionRecord <- function(input, original, fit, specification, options) {

  if (inherits(fit, "try-error"))
    stop(attr(fit, "condition"))
  rows       <- .sapCustomPredictionRows(input, fit, options, original)
  prepared   <- .sapCustomPredictionData(input, rows, fit, options)
  prediction <- .sapCustomPredictionValues(fit, prepared[["data"]], specification, options)

  return(c(list(fit = fit), prepared, prediction))
}

.sapCustomPredictionValues <- function(fit, data, specification, options) {

  levels <- specification[["levels"]]
  if (nrow(data) == 0L) {
    empty <- data.frame(estimate = numeric(0))
    for (i in seq_along(levels)) {
      empty[[paste0("lower", i)]] <- numeric(0)
      empty[[paste0("upper", i)]] <- numeric(0)
    }
    return(list(
      values   = rep(list(empty), length(specification[["measures"]])),
      warnings = gettext("No observations with supported predictor values are available for this model.")
    ))
  }

  # Reuse draws across quantities and widths, so wider intervals contain narrower
  # ones even without a user-selected seed. Clear plotting-only design overrides.
  predictionSeed <- NULL
  if (options[["setSeed"]])
    predictionSeed <- options[["seed"]]
  else if (length(levels) > 0L)
    predictionSeed <- sample.int(.Machine$integer.max, 1L)
  attr(fit, "predictionX")      <- NULL
  attr(fit, "predictionLabels") <- NULL

  values   <- list()
  warnings <- character(0)
  for (measure in specification[["measures"]]) {
    arguments <- list(
      fit     = fit,
      newdata = data,
      type    = measure[["type"]],
      ci      = length(levels) > 0L,
      B       = options[["confidenceIntervalSimulationDraws"]],
      seed    = predictionSeed
    )
    if (measure[["type"]] == "quantile")
      arguments[["quantiles"]] <- measure[["probability"]]
    record <- list()
    for (i in seq_len(max(1L, length(levels)))) {
      arguments[["cl"]] <- if (length(levels) > 0L) levels[i] else 0.95
      result   <- do.call(.sapSummaryPredictions, arguments)
      warnings <- c(warnings, attr(result, "predictionWarnings"))
      if (i == 1L)
        record[["estimate"]] <- .sapCustomPredictionVector(result, "est", fit, nrow(data))
      if (length(levels) > 0L) {
        record[[paste0("lower", i)]] <- .sapCustomPredictionVector(result, "lcl", fit, nrow(data))
        record[[paste0("upper", i)]] <- .sapCustomPredictionVector(result, "ucl", fit, nrow(data))
      }
    }
    if (any(lengths(record) != nrow(data)))
      stop(gettext("The model predictions could not be matched to the input observations."))
    values[[length(values) + 1L]] <- as.data.frame(record)
  }

  return(list(values = values, warnings = unique(warnings)))
}

.sapCustomPredictionVector <- function(result, column, fit, rows) {

  values <- unlist(lapply(result, function(row) row[[column]]), use.names = FALSE)
  # flexsurv returns one prediction for intercept-only fits, even for multiple rows.
  if (fit[["ncovs"]] == 0L && length(values) == 1L)
    values <- rep(values, rows)

  return(values)
}

.sapCustomPredictionTableColumns <- function(table, specification) {

  for (i in seq_along(specification[["measures"]])) {
    label  <- specification[["measures"]][[i]][["label"]]
    prefix <- paste0("measure", i, ".")
    table$addColumnInfo(name = paste0(prefix, "estimate"), title = gettext("Estimate"), overtitle = label, type = "number")
    for (j in seq_along(specification[["levels"]])) {
      level <- 100 * specification[["levels"]][j]
      table$addColumnInfo(name = paste0(prefix, "lower", j), title = gettextf("%1$s%% CI Lower", level), overtitle = label, type = "number")
      table$addColumnInfo(name = paste0(prefix, "upper", j), title = gettextf("%1$s%% CI Upper", level), overtitle = label, type = "number")
    }
  }

  return()
}

.sapCustomPredictionTable <- function(record, specification) {

  table <- createJaspTable(title = attr(record[["fit"]], "label"))
  table$addColumnInfo(name = "row", title = gettext("Observation"), type = "integer")
  for (j in seq_along(record[["data"]])) {
    variable <- names(record[["data"]])[j]
    table$addColumnInfo(
      name  = paste0("predictor", j),
      title = jaspBase::decodeColNames(variable),
      type  = if (is.factor(record[["data"]][[variable]])) "string" else "number"
    )
  }
  .sapCustomPredictionTableColumns(table, specification)

  data <- data.frame(row = record[["rows"]])
  for (j in seq_along(record[["data"]]))
    data[[paste0("predictor", j)]] <- record[["data"]][[j]]
  for (i in seq_along(record[["values"]])) {
    values        <- record[["values"]][[i]]
    names(values) <- paste0("measure", i, ".", names(values))
    data          <- cbind(data, values)
  }
  table$setData(data)
  table$showSpecifiedColumnsOnly <- TRUE
  for (message in .sapCustomPredictionMessages(record, specification))
    table$addFootnote(message)

  return(table)
}

.sapCustomPredictionMessages <- function(record, specification) {

  messages <- record[["warnings"]]
  if (record[["omitted"]] > 0L)
    messages <- c(messages, gettextf("%1$i observations have missing or unsupported predictor values and were not predicted.", record[["omitted"]]))
  if (.sapConstraintActive(record[["fit"]]) && length(specification[["levels"]]) > 0L)
    messages <- c(messages, gettext("Confidence intervals are unavailable because the fitted model is on the minimum-spread boundary."))

  return(messages)
}

.sapCustomPredictionColumns <- function(records, specification, original, options) {

  columns       <- list()
  prefixOptions <- options
  prefixOptions[["exportColumnPrefix"]] <- options[["customPredictionColumnPrefix"]]
  for (record in records) {
    if (inherits(record, "try-error"))
      next
    prefix <- .saExportModelPrefix(record[["fit"]], prefixOptions, length(records) > 1L)
    for (i in seq_along(record[["values"]])) {
      title <- specification[["measures"]][[i]][["label"]]
      for (field in names(record[["values"]][[i]])) {
        suffix <- if (field == "estimate") title else {
          level <- as.integer(sub("lower|upper", "", field))
          gettextf("%1$s %2$s%% CI %3$s", title, 100 * specification[["levels"]][level],
                    if (startsWith(field, "lower")) gettext("lower") else gettext("upper"))
        }
        values  <- record[["values"]][[i]][[field]]
        aligned <- rep(NA_real_, nrow(original))
        aligned[record[["rows"]]] <- ifelse(is.finite(values), values, NA_real_)
        columns[[paste0(prefix, suffix)]] <- aligned
      }
    }
  }

  return(columns)
}

.sapCustomPredictionExport <- function(container, records, specification, original, options, dependencies) {

  columns <- .sapCustomPredictionColumns(records, specification, original, options)
  # Preflight every name before writing any columns; another analysis may own it.
  for (name in names(columns))
    if (jaspBase:::columnExists(name) && !jaspBase:::columnIsMine(name))
      .quitAnalysis(gettextf("The column '%1$s' already exists. Specify a different prediction column prefix.", name))
  for (name in names(columns)) {
    container[[name]] <- createJaspColumn(columnName = name, dependencies = dependencies)
    container[[name]]$setScale(columns[[name]])
  }

  return()
}
