# Observation-level predictions use the original predictor scale and the fitted
# flexsurv terms, retaining fitted scaling constants and polynomial bases.
.sapCustomPredictionDependencies <- c(
  "customPredictions", "customPredictionSource", "customPredictionData", "customPredictionFile",
  "customPredictionMean", "customPredictionPercentiles", "customPredictionPercentileValues",
  "customPredictionConfidenceInterval", "customPredictionConfidenceLevels",
  "customPredictionTable", "customPredictionExport", "customPredictionColumnPrefix"
)

.sapCustomPredictionVariables <- function(options) {
  variables <- unique(unlist(options[c("covariates", "factors", "subgroup")], use.names = FALSE))
  return(variables[nzchar(variables)])
}

.sapCustomPredictionNumbers <- function(expression, label, lower, upper, inclusive = TRUE) {
  text <- gsub("[;|\r\n\t]+", ",", trimws(expression))
  values <- try(eval(parse(text = paste0("c(", text, ")")), envir = baseenv()), silent = TRUE)
  valid <- is.numeric(values) && !is.complex(values) && length(values) > 0L && all(is.finite(values))
  if (valid)
    valid <- if (inclusive) all(values >= lower & values <= upper) else all(values > lower & values < upper)
  if (!valid)
    stop(gettextf("%1$s must contain numeric values between %2$s and %3$s%4$s, separated by commas or semicolons.",
                  label, lower, upper, if (inclusive) "" else gettext(" (exclusive)")))
  return(unique(unname(values)))
}

.sapCustomPredictionSpecification <- function(options) {
  measures <- if (options[["customPredictionMean"]]) list(list(type = "mean", probability = NA_real_, label = gettext("Mean"))) else list()
  if (options[["customPredictionPercentiles"]]) {
    percentiles <- .sapCustomPredictionNumbers(options[["customPredictionPercentileValues"]], gettext("Percentiles"), 0, 100)
    measures <- c(measures, lapply(percentiles, function(percentile)
      list(type = "quantile", probability = percentile / 100, label = gettextf("Percentile %1$s", percentile))))
  }
  if (length(measures) == 0L)
    stop(gettext("Select a mean or at least one percentile for custom predictions."))
  levels <- if (options[["customPredictionConfidenceInterval"]])
    .sapCustomPredictionNumbers(options[["customPredictionConfidenceLevels"]], gettext("Confidence levels"), 0, 100, inclusive = FALSE) / 100 else numeric(0)
  return(list(measures = measures, levels = levels))
}

.sapCustomPredictionInput <- function(original, fitting, options, required = .sapCustomPredictionVariables(options)) {
  source    <- options[["customPredictionSource"]]
  variables <- .sapCustomPredictionVariables(options)
  if (source == "dataset") {
    data <- original
  } else if (source == "csv") {
    path <- trimws(options[["customPredictionFile"]])
    if (!nzchar(path) || !file.exists(path))
      stop(gettext("Select an existing CSV file for custom predictions."))
    data <- try(utils::read.csv(path, check.names = FALSE, colClasses = "character", na.strings = "",
                               fileEncoding = "UTF-8-BOM", strip.white = TRUE), silent = TRUE)
    if (inherits(data, "try-error"))
      stop(gettext("The prediction CSV file could not be read. Use a comma-separated file with a header row."))
    labels <- jaspBase::decodeColNames(required)
    if (anyDuplicated(names(data)))
      stop(gettext("The prediction CSV file contains duplicate column names."))
    missing <- setdiff(labels, names(data))
    if (length(missing) > 0L)
      stop(gettextf("The prediction CSV file is missing these columns: %1$s.", paste(missing, collapse = ", ")))
    data <- data[labels]
    names(data) <- required
  } else {
    columns <- options[["customPredictionData"]]
    if (length(columns) == 0L)
      stop(gettext("Add at least one observation to the prediction input table."))
    rows <- length(columns[[1L]][["values"]])
    if (length(variables) > 0L) {
      if (length(columns) != length(variables) || any(vapply(columns, function(column) length(column[["values"]]), integer(1)) != rows))
        stop(gettext("The prediction input table does not match the predictors. Reset the table and enter the observations again."))
      data <- as.data.frame(lapply(columns, function(column) as.character(unlist(column[["values"]], use.names = FALSE))), stringsAsFactors = FALSE)
      names(data) <- variables
    } else {
      data <- data.frame(row.names = seq_len(rows))
    }
  }
  if (source %in% c("dataset", "manual")) {
    # Subgroup choices must be resolved before selecting the applicable fit.
    for (variable in intersect(variables, options[["subgroup"]])) {
      blank <- is.na(data[[variable]]) | !nzchar(trimws(data[[variable]]))
      data[[variable]][blank] <- as.character(fitting[[variable]][1L])
    }
  }
  if (source == "dataset")
    return(data)
  if (nrow(data) == 0L)
    stop(gettext("The prediction input contains no observations."))
  if (.sapHasSubgroups(options)) {
    for (variable in options[["subgroup"]]) {
      invalid <- is.na(data[[variable]]) | !data[[variable]] %in% levels(fitting[[variable]])
      if (any(invalid))
        stop(gettextf("Prediction row %1$i: '%2$s' must be one of these subgroup levels: %3$s.", which(invalid)[1L],
                      jaspBase::decodeColNames(variable), paste(levels(fitting[[variable]]), collapse = ", ")))
    }
    combinations <- unique(fitting[options[["subgroup"]]])
    supported <- vapply(seq_len(nrow(data)), function(row)
      any(Reduce(`&`, lapply(options[["subgroup"]], function(variable)
        as.character(combinations[[variable]]) == as.character(data[[variable]][row])))), logical(1))
    if (any(!supported))
      stop(gettextf("Prediction row %1$i specifies a subgroup combination absent from the fitting data.", which(!supported)[1L]))
  }
  return(data)
}

.sapCustomPredictionDatasetRows <- function(dataset, options) {
  predictors <- .sapCustomPredictionVariables(options)
  available <- Reduce(`|`, lapply(predictors, function(variable)
    !is.na(dataset[[variable]]) & nzchar(trimws(as.character(dataset[[variable]])))), init = rep(FALSE, nrow(dataset)))
  if (options[["censoringType"]] == "interval") {
    # A single open bound still records a censored outcome.
    missing <- !is.finite(dataset[[options[["intervalStart"]]]]) & !is.finite(dataset[[options[["intervalEnd"]]]])
  } else {
    outcomes <- c(if (options[["censoringType"]] == "counting") c(options[["intervalStart"]], options[["intervalEnd"]]) else options[["timeToEvent"]],
                  options[["eventStatus"]])
    missing <- !stats::complete.cases(dataset[outcomes])
  }
  return(which(available & missing))
}

.sapCustomPredictionRows <- function(input, fit, options, original = input) {
  rows <- if (options[["customPredictionSource"]] == "dataset") .sapCustomPredictionDatasetRows(original, options) else seq_len(nrow(input))
  training <- attr(fit, "dataset")
  if (.sapHasSubgroups(options) && attr(fit, "subgroupLabel") != gettext("Full dataset")) {
    # Match the actual subgroup values, avoiding ambiguous concatenated labels.
    for (variable in options[["subgroup"]])
      rows <- rows[!is.na(input[[variable]][rows]) & as.character(input[[variable]][rows]) == as.character(training[[variable]][1L])]
  }
  return(rows)
}

.sapCustomPredictionData <- function(input, rows, fit, options) {
  variables <- unique(unlist(attr(fit, "modelTerms")[["components"]], use.names = FALSE))
  training  <- attr(fit, "dataset")
  source    <- options[["customPredictionSource"]]
  data      <- input[rows, variables, drop = FALSE]
  valid     <- rep(TRUE, length(rows))
  for (variable in variables) {
    label <- jaspBase::decodeColNames(variable)
    value <- as.character(data[[variable]])
    blank <- is.na(value) | !nzchar(trimws(value))
    if (is.factor(training[[variable]])) {
      levels <- levels(training[[variable]])
      if (source %in% c("dataset", "manual")) value[blank] <- levels[1L]
      invalid <- is.na(value) | !value %in% levels
      if (source != "dataset" && any(invalid))
        stop(gettextf("Prediction row %1$i: '%2$s' must be one of these levels: %3$s.", rows[which(invalid)[1L]], label, paste(levels, collapse = ", ")))
      data[[variable]] <- factor(value, levels = levels, ordered = is.ordered(training[[variable]]))
    } else {
      if (source == "manual") {
        numeric <- rep(mean(training[[variable]]), length(value))
        numeric[!blank] <- vapply(value[!blank], function(expression) {
          result <- try(eval(parse(text = expression), envir = baseenv()), silent = TRUE)
          if (!is.numeric(result) || is.complex(result) || length(result) != 1L || !is.finite(result)) return(NA_real_)
          return(as.numeric(result))
        }, numeric(1))
      } else {
        numeric <- suppressWarnings(as.numeric(value))
        if (source == "dataset") numeric[blank] <- mean(training[[variable]])
      }
      invalid <- !is.finite(numeric)
      if (source != "dataset" && any(invalid))
        stop(gettextf("Prediction row %1$i: '%2$s' requires one finite numeric value%3$s.", rows[which(invalid)[1L]], label,
                      if (source == "manual") gettext(" or expression") else ""))
      data[[variable]] <- numeric
    }
    valid <- valid & !invalid
  }
  # Validate transformations on the new values without refitting data-dependent bases.
  transformations <- attr(fit, "modelTerms")[["transformations"]]
  for (variable in names(transformations)) {
    transformation <- transformations[[variable]][["transformation"]]
    value <- suppressWarnings(as.numeric(as.character(data[[variable]])))
    invalid <- !is.finite(value)
    if (transformation %in% c("power", "inversePower", "customPower", "customInversePower")) invalid <- invalid | value <= 0
    if (transformation == "arrhenius") {
      kelvin <- value + if (transformations[[variable]][["temperatureUnit"]] == "celsius") 273.15 else 0
      invalid <- invalid | kelvin <= 0
    }
    if (transformation == "reciprocal") invalid <- invalid | value == 0
    if (transformation == "squareRoot") invalid <- invalid | value < 0
    if (transformation %in% c("square", "customExponent", "polynomial")) {
      exponent <- switch(transformation, square = 2, customExponent = transformations[[variable]][["powerExponent"]], polynomial = transformations[[variable]][["polynomialDegree"]])
      invalid <- invalid | !is.finite(suppressWarnings(value^exponent))
    }
    invalid[is.na(invalid)] <- TRUE
    if (source != "dataset" && any(invalid))
      stop(gettextf("Prediction row %1$i: the value of '%2$s' is outside the domain of its model transformation.", rows[which(invalid)[1L]], jaspBase::decodeColNames(variable)))
    valid <- valid & !invalid
  }
  return(list(data = data[valid, , drop = FALSE], rows = rows[valid], omitted = sum(!valid)))
}

.sapCustomPredictionValues <- function(fit, data, specification, options) {
  if (nrow(data) == 0L) {
    empty <- data.frame(estimate = numeric(0))
    for (i in seq_along(specification[["levels"]])) {
      empty[[paste0("lower", i)]] <- numeric(0)
      empty[[paste0("upper", i)]] <- numeric(0)
    }
    return(list(values = rep(list(empty), length(specification[["measures"]])),
                warnings = gettext("No observations with supported predictor values are available for this model.")))
  }
  # Reuse parameter draws across quantities and interval widths, so wider
  # intervals contain narrower ones even when no reproducible seed is selected.
  predictionSeed <- if (options[["setSeed"]]) options[["seed"]] else
    if (length(specification[["levels"]]) > 0L) sample.int(.Machine$integer.max, 1L)
  values   <- list()
  warnings <- character(0)
  for (measure in specification[["measures"]]) {
    arguments <- list(fit = fit, newdata = data, type = measure[["type"]], ci = length(specification[["levels"]]) > 0L,
                      B = options[["confidenceIntervalSimulationDraws"]], seed = predictionSeed)
    if (measure[["type"]] == "quantile") arguments[["quantiles"]] <- measure[["probability"]]
    # Clear plotting-only design overrides: custom observations use fitted terms.
    attr(arguments[["fit"]], "predictionX") <- NULL
    attr(arguments[["fit"]], "predictionLabels") <- NULL
    record <- list()
    for (i in seq_len(max(1L, length(specification[["levels"]])))) {
      arguments[["cl"]] <- if (length(specification[["levels"]]) > 0L) specification[["levels"]][i] else 0.95
      result <- do.call(.sapSummaryPredictions, arguments)
      warnings <- c(warnings, attr(result, "predictionWarnings"))
      extract <- function(column) {
        out <- unlist(lapply(result, function(row) row[[column]]), use.names = FALSE)
        # Intercept-only fits return one prediction, regardless of newdata rows.
        return(if (fit[["ncovs"]] == 0L && length(out) == 1L) rep(out, nrow(data)) else out)
      }
      if (i == 1L) record[["estimate"]] <- extract("est")
      if (length(specification[["levels"]]) > 0L) {
        record[[paste0("lower", i)]] <- extract("lcl")
        record[[paste0("upper", i)]] <- extract("ucl")
      }
    }
    if (any(lengths(record) != nrow(data)))
      stop(gettext("The model predictions could not be matched to the input observations."))
    values[[length(values) + 1L]] <- as.data.frame(record)
  }
  return(list(values = values, warnings = unique(warnings)))
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

.sapCustomPredictionTable <- function(record, specification, options) {
  table <- createJaspTable(title = attr(record[["fit"]], "label"))
  table$addColumnInfo(name = "row", title = gettext("Observation"), type = "integer")
  for (variable in names(record[["data"]]))
    table$addColumnInfo(name = paste0("predictor", match(variable, names(record[["data"]]))), title = jaspBase::decodeColNames(variable),
                        type = if (is.factor(record[["data"]][[variable]])) "string" else "number")
  .sapCustomPredictionTableColumns(table, specification)
  data <- data.frame(row = record[["rows"]])
  for (j in seq_along(record[["data"]])) data[[paste0("predictor", j)]] <- record[["data"]][[j]]
  for (i in seq_along(record[["values"]])) {
    values <- record[["values"]][[i]]
    names(values) <- paste0("measure", i, ".", names(values))
    data <- cbind(data, values)
  }
  table$setData(data)
  table$showSpecifiedColumnsOnly <- TRUE
  for (message in .sapCustomPredictionMessages(record, specification)) table$addFootnote(message)
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
  columns <- list()
  prefixOptions <- options
  prefixOptions[["exportColumnPrefix"]] <- options[["customPredictionColumnPrefix"]]
  for (record in records) {
    if (inherits(record, "try-error")) next
    prefix <- .saExportModelPrefix(record[["fit"]], prefixOptions, length(records) > 1L)
    for (i in seq_along(record[["values"]])) {
      title <- specification[["measures"]][[i]][["label"]]
      for (field in names(record[["values"]][[i]])) {
        suffix <- if (field == "estimate") title else {
          level <- as.integer(sub("lower|upper", "", field))
          gettextf("%1$s %2$s%% CI %3$s", title, 100 * specification[["levels"]][level],
                    if (startsWith(field, "lower")) gettext("lower") else gettext("upper"))
        }
        values <- record[["values"]][[i]][[field]]
        aligned <- rep(NA_real_, nrow(original))
        aligned[record[["rows"]]] <- ifelse(is.finite(values), values, NA_real_)
        columns[[paste0(prefix, suffix)]] <- aligned
      }
    }
  }
  return(columns)
}

.sapCustomPredictions <- function(jaspResults, original, fitting, options) {
  if (!options[["customPredictionTable"]] && !(options[["customPredictionSource"]] == "dataset" && options[["customPredictionExport"]])) return()
  if (!is.null(jaspResults[["customPredictionsResult"]])) return()
  dependencies <- c(.sapGetDependencies(options), .sapSelectedOutputDependencies, .sapSimulationDependencies, .sapCustomPredictionDependencies)
  container <- createJaspContainer(title = gettext("Custom Predictions"))
  container$dependOn(dependencies)
  container$position <- 3.5
  jaspResults[["customPredictionsResult"]] <- container
  if (!.saSurvivalReady(options)) {
    table <- createJaspTable(title = gettext("Predictions"))
    table$addColumnInfo(name = "row", title = gettext("Observation"), type = "integer")
    container[["predictions"]] <- table
    specification <- try(.sapCustomPredictionSpecification(options), silent = TRUE)
    if (inherits(specification, "try-error")) table$setError(conditionMessage(attr(specification, "condition"))) else
      .sapCustomPredictionTableColumns(table, specification)
    return()
  }
  fits <- .sapFlattenFit(.sapExtractFit(jaspResults, options, type = "selected"), options)
  required <- unique(c(unlist(lapply(fits, function(fit) attr(fit, "modelTerms")[["components"]]), use.names = FALSE),
                       unlist(options[["subgroup"]], use.names = FALSE)))
  setup <- try(list(specification = .sapCustomPredictionSpecification(options), input = .sapCustomPredictionInput(original, fitting, options, required)), silent = TRUE)
  if (inherits(setup, "try-error")) {
    table <- createJaspTable(title = gettext("Predictions"))
    container[["inputError"]] <- table
    table$setError(conditionMessage(attr(setup, "condition")))
    return()
  }
  records <- list()
  messages <- character(0)
  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    record <- try({
      if (inherits(fit, "try-error")) stop(attr(fit, "condition"))
      rows <- .sapCustomPredictionRows(setup[["input"]], fit, options, original)
      prepared <- .sapCustomPredictionData(setup[["input"]], rows, fit, options)
      prediction <- .sapCustomPredictionValues(fit, prepared[["data"]], setup[["specification"]], options)
      c(list(fit = fit), prepared, prediction)
    }, silent = TRUE)
    records[[i]] <- record
    if (inherits(record, "try-error")) {
      table <- createJaspTable(title = attr(fit, "label"))
      container[[paste0("model", i)]] <- table
      table$setError(conditionMessage(attr(record, "condition")))
    } else if (options[["customPredictionTable"]]) {
      container[[paste0("model", i)]] <- .sapCustomPredictionTable(record, setup[["specification"]], options)
    } else {
      notes <- .sapCustomPredictionMessages(record, setup[["specification"]])
      if (length(notes) > 0L)
        messages <- c(messages, paste0(if (nzchar(attr(fit, "label"))) paste0(attr(fit, "label"), ": ") else "", notes))
    }
  }
  if (length(messages) > 0L) {
    table <- createJaspTable(title = gettext("Prediction Notes"))
    container[["messages"]] <- table
    .saFillMessageTable(table, unique(messages))
  }
  if (options[["customPredictionSource"]] == "dataset" && options[["customPredictionExport"]]) {
    columns <- .sapCustomPredictionColumns(records, setup[["specification"]], original, options)
    for (name in names(columns))
      if (jaspBase:::columnExists(name) && !jaspBase:::columnIsMine(name))
        .quitAnalysis(gettextf("The column '%1$s' already exists. Specify a different prediction column prefix.", name))
    for (name in names(columns)) {
      container[[name]] <- createJaspColumn(columnName = name, dependencies = dependencies)
      container[[name]]$setScale(columns[[name]])
    }
  }
  return()
}
