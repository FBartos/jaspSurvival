# A variable has one transformation per model; interactions share that choice.
.sapPrepareModelTerms <- function(options) {

  for (i in seq_along(options[["modelTerms"]])) {
    model <- options[["modelTerms"]][[i]]
    rows  <- model[["components"]]
    model[["components"]] <- lapply(rows, function(row) row[["components"]])
    transformations <- list()
    if (options[["advancedSpecification"]]) {
      for (row in rows) {
        if (length(row[["components"]]) == 1 && row[["transformation"]] != "none")
          transformations[[row[["components"]][[1]]]] <- list(
            transformation = row[["transformation"]], temperatureUnit = row[["temperatureUnit"]]
          )
      }
    }
    model[["transformations"]] <- transformations
    options[["modelTerms"]][[i]] <- model
  }

  return(options)
}

.sapVariableExpression <- function(variable, transformation) {

  if (is.null(transformation))
    return(variable)

  numeric <- paste0("as.numeric(as.character(", variable, "))")
  expression <- switch(transformation[["transformation"]],
    exponential  = paste0("I(", numeric, ")"),
    power        = paste0("log(", numeric, ")"),
    inversePower = paste0("I(-log(", numeric, "))"),
    arrhenius    = paste0("I(1/(", numeric, if (transformation[["temperatureUnit"]] == "celsius") " + 273.15", "))"),
    reciprocal   = paste0("I(1/", numeric, ")"),
    square       = paste0("I(", numeric, "^2)"),
    squareRoot   = paste0("sqrt(", numeric, ")")
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
    if (transformation %in% c("power", "inversePower") && any(value <= 0))
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
  }

  return()
}

.sapFittedFactors <- function(fit, options) {

  variables <- unlist(attr(fit, "modelTerms")[["components"]], use.names = FALSE)
  return(intersect(variables, unlist(options[["factors"]], use.names = FALSE)))
}

.sapTransformationNames <- function(names, fit) {

  transformations <- attr(fit, "modelTerms")[["transformations"]]
  for (variable in names(transformations)) {
    specification <- transformations[[variable]]
    expression    <- .sapVariableExpression(variable, specification)
    label <- switch(specification[["transformation"]],
      exponential = gettext("Exponential"), power = gettext("Power"), inversePower = gettext("Inverse power"),
      arrhenius = gettext("Arrhenius"), reciprocal = gettext("Reciprocal"),
      square = gettext("Square"), squareRoot = gettext("Square root")
    )
    names <- gsub(expression, paste0(label, "(", jaspBase::decodeColNames(variable), ")"), names, fixed = TRUE)
  }

  return(names)
}

.sapModelsNested <- function(fit0, fit1) {

  if (fit1[["npars"]] <= fit0[["npars"]])
    return(FALSE)

  design0 <- cbind(1, stats::model.matrix(fit0))
  design1 <- cbind(1, stats::model.matrix(fit1))
  if (!identical(rownames(design0), rownames(design1)))
    return(FALSE)

  normalize <- function(design) {
    norm <- sqrt(colSums(design^2))
    norm[norm == 0] <- 1
    return(sweep(design, 2, norm, "/"))
  }
  design0 <- normalize(design0)
  design1 <- normalize(design1)
  return(qr(cbind(design1, design0))[["rank"]] == qr(design1)[["rank"]])
}
