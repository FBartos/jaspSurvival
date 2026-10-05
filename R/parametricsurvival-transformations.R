# A variable has one transformation per model; interactions share that choice.
.sapPrepareModelTerms <- function(options) {

  for (i in seq_along(options[["modelTerms"]])) {
    model <- options[["modelTerms"]][[i]]
    transformations <- list()
    restrictions    <- rep(list(list()), max(.sapComponents(options)))
    if (options[["advancedSpecification"]]) {
      modifierIndex <- which(vapply(options[["modelTermModifiers"]], function(modifier)
        modifier[["model"]] == model[["name"]], logical(1)))
      for (modifier in options[["modelTermModifiers"]][modifierIndex]) {
        for (row in modifier[["predictors"]]) {
          variable <- row[["variable"]][[1L]]
          if (options[["analysisType"]] == "mixture") {
            componentRows <- row[["componentRestrictions"]]
            componentKeys <- vapply(componentRows, function(component) component[["component"]][[1L]], character(1))
            for (k in seq_along(restrictions)) {
              component <- componentRows[[match(as.character(k), componentKeys)]]
              restriction <- .sapParseParameterRestriction(component[["parameterRestriction"]], variable, gettextf("%1$s, component %2$i", model[["title"]], k))
              if (length(restriction) > 0L)
                restrictions[[k]][[variable]] <- restriction
            }
          } else {
            restriction <- .sapParseParameterRestriction(row[["parameterRestriction"]], variable, model[["title"]])
            if (length(restriction) > 0L)
              restrictions[[1L]][[variable]] <- restriction
          }
          if (row[["transformation"]] != "none") {
            specification <- list(
              transformation = row[["transformation"]], temperatureUnit = row[["temperatureUnit"]]
            )
            if (row[["transformation"]] %in% c("customPower", "customInversePower")) {
              base <- row[["customLogBase"]]
              if (base == 1)
                .quitAnalysis(gettext("The custom logarithm base must be different from 1."))
              specification[["logBase"]] <- base
            }
            if (row[["transformation"]] == "customExponent") {
              if (row[["powerExponent"]] == 0)
                .quitAnalysis(gettext("The custom exponent must be different from 0."))
              specification[["powerExponent"]] <- row[["powerExponent"]]
            }
            if (row[["transformation"]] == "polynomial")
              specification[["polynomialDegree"]] <- row[["polynomialDegree"]]
            transformations[[variable]] <- specification
          }
        }
      }
    }
    model[["transformations"]] <- transformations
    model[["parameterRestrictions"]] <- restrictions
    options[["modelTerms"]][[i]] <- model
  }

  return(options)
}

.sapVariableExpression <- function(variable, transformation, dataset) {

  if (is.null(transformation))
    return(variable)

  numeric <- paste0("as.numeric(as.character(", variable, "))")
  expression <- switch(transformation[["transformation"]],
    scale        = paste0("I((", variable, " - ", format(mean(dataset[[variable]]), digits = 17, scientific = TRUE, trim = TRUE),
                          ") / ", format(stats::sd(dataset[[variable]]), digits = 17, scientific = TRUE, trim = TRUE), ")"),
    power        = paste0("log(", numeric, ")"),
    customPower  = paste0("log(", numeric, ", base = ", format(transformation[["logBase"]], digits = 17, scientific = FALSE, trim = TRUE), ")"),
    inversePower = paste0("I(-log(", numeric, "))"),
    customInversePower = paste0("I(-log(", numeric, ", base = ", format(transformation[["logBase"]], digits = 17, scientific = FALSE, trim = TRUE), "))"),
    arrhenius    = paste0("I(1/(", numeric, if (transformation[["temperatureUnit"]] == "celsius") " + 273.15", "))"),
    reciprocal   = paste0("I(1/", numeric, ")"),
    square       = paste0("I(", numeric, "^2)"),
    customExponent = paste0("I(", numeric, "^", format(transformation[["powerExponent"]], digits = 17, scientific = FALSE, trim = TRUE), ")"),
    squareRoot   = paste0("sqrt(", numeric, ")"),
    polynomial   = paste0("stats::poly(", numeric, ", degree = ", transformation[["polynomialDegree"]], ", raw = FALSE)")
  )

  return(paste(deparse(str2lang(expression), width.cutoff = 500L), collapse = ""))
}

.sapCheckTransformationData <- function(dataset, modelTerms) {

  for (variable in names(modelTerms[["transformations"]])) {
    specification  <- modelTerms[["transformations"]][[variable]]
    transformation <- specification[["transformation"]]
    value          <- suppressWarnings(as.numeric(as.character(dataset[[variable]])))
    label          <- jaspBase::decodeColNames(variable)
    if (any(!is.finite(value)))
      .quitAnalysis(gettextf("The transformation of '%1$s' requires finite numeric values or numeric factor-level labels.", label))
    if (transformation == "scale") {
      deviation <- stats::sd(dataset[[variable]])
      if (!is.finite(deviation) || deviation == 0)
        .quitAnalysis(gettextf("Scaling '%1$s' requires a finite, nonzero standard deviation.", label))
    }
    if (transformation %in% c("power", "inversePower", "customPower", "customInversePower") && any(value <= 0))
      .quitAnalysis(gettextf("The logarithmic transformation of '%1$s' requires strictly positive stress values.", label))
    if (transformation == "arrhenius") {
      kelvin <- value + if (specification[["temperatureUnit"]] == "celsius") 273.15 else 0
      if (any(kelvin <= 0))
        .quitAnalysis(gettextf("The Arrhenius transformation of '%1$s' requires temperatures above absolute zero. Check the selected temperature units.", label))
    }
    if (transformation == "reciprocal" && any(value == 0))
      .quitAnalysis(gettextf("The reciprocal transformation of '%1$s' requires nonzero stress values.", label))
    if (transformation == "squareRoot" && any(value < 0))
      .quitAnalysis(gettextf("The square-root transformation of '%1$s' requires nonnegative stress values.", label))
    if (transformation == "square" && any(!is.finite(value^2)))
      .quitAnalysis(gettextf("The squared stress values of '%1$s' exceed the numeric range.", label))
    if (transformation == "customExponent") {
      exponent <- specification[["powerExponent"]]
      if (exponent < 0 && any(value == 0))
        .quitAnalysis(gettextf("A negative exponent for '%1$s' requires nonzero values.", label))
      if (exponent != round(exponent) && any(value < 0))
        .quitAnalysis(gettextf("A fractional exponent for '%1$s' requires nonnegative values.", label))
      if (any(!is.finite(value^exponent)))
        .quitAnalysis(gettextf("The power-transformed values of '%1$s' exceed the numeric range.", label))
    }
    if (transformation == "polynomial") {
      degree <- specification[["polynomialDegree"]]
      if (length(unique(value)) <= degree)
        .quitAnalysis(gettextf("The polynomial transformation of '%1$s' requires at least %2$s distinct values.", label, degree + 1))
      if (any(!is.finite(value^degree)))
        .quitAnalysis(gettextf("The polynomial transformation of '%1$s' exceeds the numeric range.", label))
    }
  }

  return()
}

.sapFittedFactors <- function(fit, options) {

  variables <- unlist(attr(fit, "modelTerms")[["components"]], use.names = FALSE)
  return(intersect(variables, unlist(options[["factors"]], use.names = FALSE)))
}

.sapTransformationNames <- function(names, fit) {

  transformations <- attr(fit, "modelTerms")[["transformations"]]
  variables <- unique(unlist(attr(fit, "modelTerms")[["components"]], use.names = FALSE))
  for (i in seq_along(names)) {
    # Keep namespace operators intact; only split interaction separators.
    originalTerms <- strsplit(names[i], "(?<!:):(?!:)", perl = TRUE)[[1]]
    displayTerms <- originalTerms
    for (variable in names(transformations)) {
      specification <- transformations[[variable]]
      expression <- .sapVariableExpression(variable, specification, attr(fit, "dataset"))
      variableLabel <- jaspBase::decodeColNames(variable)
      if (specification[["transformation"]] == "polynomial") {
        if (specification[["polynomialDegree"]] == 1)
          displayTerms <- gsub(expression, gettextf("Polynomial(%1$s, %2$s)", variableLabel, 1), displayTerms, fixed = TRUE)
        else
          for (degree in seq_len(specification[["polynomialDegree"]]))
            displayTerms <- gsub(paste0(expression, degree), gettextf("Polynomial(%1$s, %2$s)", variableLabel, degree), displayTerms, fixed = TRUE)
      } else {
        label <- switch(specification[["transformation"]],
          scale = gettext("Scale"),
          power = gettext("Power"), customPower = gettext("Power"),
          inversePower = gettext("Inverse power"), customInversePower = gettext("Inverse power"),
          customExponent = gettext("Exponent"),
          arrhenius = gettext("Arrhenius"), reciprocal = gettext("Reciprocal"),
          square = gettext("Square"), squareRoot = gettext("Square root")
        )
        label <- if (specification[["transformation"]] %in% c("customPower", "customInversePower"))
          gettextf("%1$s(%2$s, base %3$s)", label, variableLabel, specification[["logBase"]]) else
          if (specification[["transformation"]] == "customExponent")
            gettextf("%1$s(%2$s, %3$s)", label, variableLabel, specification[["powerExponent"]]) else
            paste0(label, "(", variableLabel, ")")
        displayTerms <- gsub(expression, label, displayTerms, fixed = TRUE)
      }
    }
    ordinaryTerms <- displayTerms == originalTerms
    if (!all(ordinaryTerms)) {
      displayTerms[ordinaryTerms] <- vapply(originalTerms[ordinaryTerms], .saTermNames, character(1), variables = variables)
      names[i] <- paste(jaspBase::decodeColNames(displayTerms), collapse = jaspBase::interactionSymbol)
    }
  }

  return(names)
}

.sapModelsNested <- function(fit0, fit1) {

  if (fit1[["npars"]] <= fit0[["npars"]])
    return(FALSE)

  spaces0 <- .sapRegressionSpaces(fit0)
  spaces1 <- .sapRegressionSpaces(fit1)
  if (length(spaces0) != length(spaces1)) return(FALSE)
  normalize <- function(design) {
    norm <- sqrt(colSums(design^2))
    norm[norm == 0] <- 1
    return(sweep(design, 2, norm, "/"))
  }
  rank <- function(design) if (ncol(design) == 0L) 0L else qr(normalize(design))[["rank"]]
  for (k in seq_along(spaces0)) {
    if (!identical(rownames(spaces0[[k]][["design"]]), rownames(spaces1[[k]][["design"]]))) return(FALSE)
    design0 <- spaces0[[k]][["design"]]
    design1 <- spaces1[[k]][["design"]]
    difference <- spaces0[[k]][["offset"]] - spaces1[[k]][["offset"]]
    # Fixed coefficients create an affine offset, which must also be representable.
    if (rank(cbind(design1, design0, difference)) != rank(design1)) return(FALSE)
  }
  return(TRUE)
}
