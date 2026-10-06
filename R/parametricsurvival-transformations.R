# A variable has one transformation per model; interactions share that choice.
.sapPrepareModelTerms <- function(options) {

  components <- max(.sapComponents(options))
  mixture    <- options[["analysisType"]] == "mixture"
  for (i in seq_along(options[["modelTerms"]])) {
    model           <- options[["modelTerms"]][[i]]
    transformations <- list()
    restrictions    <- rep(list(list()), components)
    if (options[["advancedSpecification"]]) {
      modifierIndex <- which(vapply(options[["modelTermModifiers"]], function(modifier)
        modifier[["model"]] == model[["name"]], logical(1)))
      rows <- unlist(lapply(options[["modelTermModifiers"]][modifierIndex], function(modifier)
        modifier[["predictors"]]), recursive = FALSE, use.names = FALSE)
      for (row in rows) {
        variable <- row[["variable"]][[1L]]
        values <- .sapModelTermRestrictions(row, variable, model[["title"]], components, mixture)
        for (k in seq_along(values))
          if (length(values[[k]]) > 0L)
            restrictions[[k]][[variable]] <- values[[k]]

        specification <- .sapTransformationSpecification(row)
        if (!is.null(specification))
          transformations[[variable]] <- specification
      }
    }
    model[["transformations"]]       <- transformations
    model[["parameterRestrictions"]] <- restrictions
    options[["modelTerms"]][[i]]      <- model
  }

  return(options)
}

.sapTransformationSpecification <- function(row) {

  transformation <- row[["transformation"]]
  if (transformation == "none")
    return(NULL)

  specification <- list(transformation = transformation, temperatureUnit = row[["temperatureUnit"]])
  if (transformation %in% c("customPower", "customInversePower")) {
    if (row[["customLogBase"]] == 1)
      .quitAnalysis(gettext("The custom logarithm base must be different from 1."))
    specification[["logBase"]] <- row[["customLogBase"]]
  }
  if (transformation == "customExponent") {
    if (row[["powerExponent"]] == 0)
      .quitAnalysis(gettext("The custom exponent must be different from 0."))
    specification[["powerExponent"]] <- row[["powerExponent"]]
  }
  if (transformation == "polynomial")
    specification[["polynomialDegree"]] <- row[["polynomialDegree"]]

  return(specification)
}

.sapVariableExpression <- function(variable, transformation, dataset) {

  if (is.null(transformation))
    return(variable)

  numericExpression <- paste0("as.numeric(as.character(", variable, "))")
  number <- function(value, scientific = FALSE)
    format(value, digits = 17, scientific = scientific, trim = TRUE)
  expression <- switch(transformation[["transformation"]],
    # Embed fitting-data constants, so predictions do not recompute scaling.
    scale        = paste0("I((", variable, " - ", number(mean(dataset[[variable]]), scientific = TRUE),
                          ") / ", number(stats::sd(dataset[[variable]]), scientific = TRUE), ")"),
    power        = paste0("log(", numericExpression, ")"),
    customPower  = paste0("log(", numericExpression, ", base = ", number(transformation[["logBase"]]), ")"),
    inversePower = paste0("I(-log(", numericExpression, "))"),
    customInversePower = paste0("I(-log(", numericExpression, ", base = ", number(transformation[["logBase"]]), "))"),
    arrhenius    = paste0("I(1/(", numericExpression, if (transformation[["temperatureUnit"]] == "celsius") " + 273.15", "))"),
    reciprocal   = paste0("I(1/", numericExpression, ")"),
    square       = paste0("I(", numericExpression, "^2)"),
    customExponent = paste0("I(", numericExpression, "^", number(transformation[["powerExponent"]]), ")"),
    squareRoot   = paste0("sqrt(", numericExpression, ")"),
    polynomial   = paste0("stats::poly(", numericExpression, ", degree = ", transformation[["polynomialDegree"]], ", raw = FALSE)")
  )

  return(paste(deparse(str2lang(expression), width.cutoff = 500L), collapse = ""))
}

# Per-row domain checks are shared by fitting validation and custom predictions.
# Sample-level checks (variance and polynomial degree) belong to fitting only.
.sapTransformationInvalid <- function(value, specification) {

  domain <- switch(specification[["transformation"]],
    power =, inversePower =, customPower =, customInversePower = value <= 0,
    arrhenius = (value + if (specification[["temperatureUnit"]] == "celsius") 273.15 else 0) <= 0,
    reciprocal = value == 0,
    squareRoot = value < 0,
    square = !is.finite(suppressWarnings(value^2)),
    customExponent = !is.finite(suppressWarnings(value^specification[["powerExponent"]])),
    polynomial = !is.finite(suppressWarnings(value^specification[["polynomialDegree"]])),
    rep(FALSE, length(value))
  )
  invalid <- !is.finite(value) | domain
  invalid[is.na(invalid)] <- TRUE

  return(invalid)
}

.sapCheckTransformationData <- function(dataset, modelTerms) {

  for (variable in names(modelTerms[["transformations"]])) {
    specification  <- modelTerms[["transformations"]][[variable]]
    transformation <- specification[["transformation"]]
    value          <- suppressWarnings(as.numeric(as.character(dataset[[variable]])))
    label          <- jaspBase::decodeColNames(variable)
    invalid        <- .sapTransformationInvalid(value, specification)
    if (any(!is.finite(value)))
      .quitAnalysis(gettextf("The transformation of '%1$s' requires finite numeric values or numeric factor-level labels.", label))
    if (transformation == "scale") {
      deviation <- stats::sd(dataset[[variable]])
      if (!is.finite(deviation) || deviation == 0)
        .quitAnalysis(gettextf("Scaling '%1$s' requires a finite, nonzero standard deviation.", label))
    }
    if (transformation %in% c("power", "inversePower", "customPower", "customInversePower") && any(invalid))
      .quitAnalysis(gettextf("The logarithmic transformation of '%1$s' requires strictly positive stress values.", label))
    if (transformation == "arrhenius") {
      if (any(invalid))
        .quitAnalysis(gettextf("The Arrhenius transformation of '%1$s' requires temperatures above absolute zero. Check the selected temperature units.", label))
    }
    if (transformation == "reciprocal" && any(invalid))
      .quitAnalysis(gettextf("The reciprocal transformation of '%1$s' requires nonzero stress values.", label))
    if (transformation == "squareRoot" && any(invalid))
      .quitAnalysis(gettextf("The square-root transformation of '%1$s' requires nonnegative stress values.", label))
    if (transformation == "square" && any(invalid))
      .quitAnalysis(gettextf("The squared stress values of '%1$s' exceed the numeric range.", label))
    if (transformation == "customExponent") {
      exponent <- specification[["powerExponent"]]
      if (exponent < 0 && any(value == 0))
        .quitAnalysis(gettextf("A negative exponent for '%1$s' requires nonzero values.", label))
      if (exponent != round(exponent) && any(value < 0))
        .quitAnalysis(gettextf("A fractional exponent for '%1$s' requires nonnegative values.", label))
      if (any(invalid))
        .quitAnalysis(gettextf("The power-transformed values of '%1$s' exceed the numeric range.", label))
    }
    if (transformation == "polynomial") {
      degree <- specification[["polynomialDegree"]]
      if (length(unique(value)) <= degree)
        .quitAnalysis(gettextf("The polynomial transformation of '%1$s' requires at least %2$s distinct values.", label, degree + 1))
      if (any(invalid))
        .quitAnalysis(gettextf("The polynomial transformation of '%1$s' exceeds the numeric range.", label))
    }
  }

  return()
}

.sapFittedFactors <- function(fit, options) {

  variables <- unlist(attr(fit, "modelTerms")[["components"]], use.names = FALSE)
  return(intersect(variables, unlist(options[["factors"]], use.names = FALSE)))
}

.sapTransformationNames <- function(coefficients, fit) {

  transformations <- attr(fit, "modelTerms")[["transformations"]]
  variables       <- unique(unlist(attr(fit, "modelTerms")[["components"]], use.names = FALSE))
  for (i in seq_along(coefficients)) {
    # Keep namespace operators intact; only split interaction separators.
    originalTerms <- strsplit(coefficients[i], "(?<!:):(?!:)", perl = TRUE)[[1]]
    displayTerms  <- originalTerms
    for (variable in names(transformations)) {
      specification <- transformations[[variable]]
      expression    <- .sapVariableExpression(variable, specification, attr(fit, "dataset"))
      variableLabel <- jaspBase::decodeColNames(variable)
      if (specification[["transformation"]] == "polynomial") {
        if (specification[["polynomialDegree"]] == 1)
          displayTerms <- gsub(expression, gettextf("Polynomial(%1$s, %2$s)", variableLabel, 1), displayTerms, fixed = TRUE)
        else
          for (degree in seq_len(specification[["polynomialDegree"]]))
            displayTerms <- gsub(paste0(expression, degree), gettextf("Polynomial(%1$s, %2$s)", variableLabel, degree), displayTerms, fixed = TRUE)
      } else {
        label <- .sapTransformationLabel(variableLabel, specification)
        displayTerms <- gsub(expression, label, displayTerms, fixed = TRUE)
      }
    }
    ordinaryTerms <- displayTerms == originalTerms
    if (!all(ordinaryTerms)) {
      displayTerms[ordinaryTerms] <- vapply(originalTerms[ordinaryTerms], .saTermNames, character(1), variables = variables)
      coefficients[i] <- paste(jaspBase::decodeColNames(displayTerms), collapse = jaspBase::interactionSymbol)
    }
  }

  return(coefficients)
}

.sapTransformationLabel <- function(variable, specification) {

  transformation <- specification[["transformation"]]
  label <- switch(transformation,
    scale              = gettext("Scale"),
    power              = gettext("Power"),
    customPower        = gettext("Power"),
    inversePower       = gettext("Inverse power"),
    customInversePower = gettext("Inverse power"),
    customExponent     = gettext("Exponent"),
    arrhenius          = gettext("Arrhenius"),
    reciprocal         = gettext("Reciprocal"),
    square             = gettext("Square"),
    squareRoot         = gettext("Square root")
  )
  if (transformation %in% c("customPower", "customInversePower"))
    return(gettextf("%1$s(%2$s, base %3$s)", label, variable, specification[["logBase"]]))
  if (transformation == "customExponent")
    return(gettextf("%1$s(%2$s, %3$s)", label, variable, specification[["powerExponent"]]))

  return(paste0(label, "(", variable, ")"))
}
