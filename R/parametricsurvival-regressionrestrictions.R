# Restrictions use the coefficients displayed for each main-effect term.
.sapParseParameterRestriction <- function(expression, variable, modelTitle) {

  if (is.character(expression) && length(expression) == 1L && trimws(expression) == "")
    return(numeric(0))

  value <- .sapNumericExpression(expression, multiple = TRUE, allowSpaces = TRUE)
  if (!.sapFiniteNumeric(value))
    .quitAnalysis(gettextf("The parameter restriction for '%1$s' in '%2$s' must contain finite numeric values or expressions, separated by commas or semicolons.",
                          jaspBase::decodeColNames(variable), modelTitle))

  return(unname(value))
}

.sapModelTermRestrictions <- function(row, variable, modelTitle, components, mixture) {

  if (!mixture)
    return(list(.sapParseParameterRestriction(row[["parameterRestriction"]], variable, modelTitle)))

  rows <- row[["componentRestrictions"]]
  keys <- vapply(rows, function(component) component[["component"]][[1L]], character(1))
  return(lapply(seq_len(components), function(k) {
    component <- rows[[match(as.character(k), keys)]]
    title     <- gettextf("%1$s, component %2$i", modelTitle, k)
    return(.sapParseParameterRestriction(component[["parameterRestriction"]], variable, title))
  }))
}
.sapValidateParameterRestrictions <- function(dataset, options) {

  models <- options[["modelTerms"]]
  if (!any(vapply(models, function(model) any(lengths(model[["parameterRestrictions"]]) > 0L), logical(1))))
    return()

  validate <- function(data) {
    for (model in models)
      for (k in seq_len(max(.sapComponents(options))))
        .sapRegressionFixed(data, model, .sapGetFormula(options, model, data), k)
  }
  if (!.sapHasSubgroups(options) || options[["includeFullDatasetInSubgroupAnalysis"]])
    validate(dataset)
  if (.sapHasSubgroups(options)) {
    subgroups <- .sapSubgroupFactor(dataset, options)
    for (group in unique(subgroups))
      validate(droplevels(dataset[subgroups == group, , drop = FALSE]))
  }

  return()
}
.sapRegressionFixed <- function(dataset, modelTerms, formula, component = 1L) {

  if (length(modelTerms[["parameterRestrictions"]]) == 0L)
    return(numeric(0))

  restrictions <- modelTerms[["parameterRestrictions"]][[component]]
  if (length(restrictions) == 0L)
    return(numeric(0))

  .sapCheckTransformationData(dataset, modelTerms)
  frame    <- stats::model.frame(formula, dataset)
  design   <- stats::model.matrix(formula, frame)
  labels   <- attr(stats::terms(formula), "term.labels")
  assigned <- attr(design, "assign")
  fixed    <- numeric(0)

  modelTitle <- if (length(modelTerms[["parameterRestrictions"]]) > 1L) gettextf("%1$s, component %2$i", modelTerms[["title"]], component) else modelTerms[["title"]]
  for (variable in names(restrictions)) {
    transformation <- modelTerms[["transformations"]][[variable]]
    expression <- if (is.null(transformation)) variable else .sapVariableExpression(variable, transformation, dataset)
    term       <- match(expression, labels)
    index      <- which(assigned == term)
    value      <- restrictions[[variable]]
    label      <- jaspBase::decodeColNames(variable)
    if (length(index) == 0L)
      .quitAnalysis(gettextf("No main-effect regression coefficient for '%1$s' is present in '%2$s'. Add its main effect or leave its parameter restriction empty.", label, modelTitle))
    if (length(value) != 1L && length(value) != length(index)) {
      if (is.factor(dataset[[variable]]) && is.null(transformation))
        .quitAnalysis(gettextf("Factor '%1$s' in '%2$s' has %3$i regression coefficients (%4$i levels). Supply one value or %3$i values; %5$i were supplied.",
                              label, modelTitle, length(index), nlevels(dataset[[variable]]), length(value)))
      else
        .quitAnalysis(gettextf("Term '%1$s' in '%2$s' has %3$i regression coefficients. Supply one value or %3$i values; %4$i were supplied.",
                              label, modelTitle, length(index), length(value)))
    }
    fixed[colnames(design)[index]] <- rep(value, length.out = length(index))
  }

  attr(fixed, "parameterNames") <- colnames(design)[-1L]
  return(fixed)
}
.sapModelTermsForComponents <- function(modelTerms, components) {

  if (length(modelTerms[["parameterRestrictions"]]) > components)
    modelTerms[["parameterRestrictions"]] <- modelTerms[["parameterRestrictions"]][seq_len(components)]
  return(modelTerms)
}
.sapmRegressionFixed <- function(regression, family, components) {

  if (length(regression) == 0L)
    return(numeric(0))

  return(unlist(lapply(seq_len(components), function(k) {
    if (length(regression[[k]]) == 0L)
      return(numeric(0))
    parameters <- .sapmRegressionParameterNames(names(regression[[k]]), family, k)
    return(stats::setNames(as.numeric(regression[[k]]), parameters))
  }), use.names = TRUE))
}
.sapmParameterOrder <- function(mixture, family, components, covariates) {

  if (length(covariates) == 0L)
    return(mixture[["dlist"]][["pars"]])

  regression <- unlist(lapply(seq_len(components), function(k)
    .sapmRegressionParameterNames(covariates, family, k)), use.names = FALSE)
  return(c(mixture[["dlist"]][["pars"]], regression))
}

.sapmRegressionParameterNames <- function(covariates, family, component) {

  # flexsurv leaves the first component's covariate names bare; later components
  # are ancillary formulas and use locationK(covariate) names.
  if (component == 1L)
    return(covariates)
  return(paste0(family[["location"]], component, "(", covariates, ")"))
}
.sapRegressionSpaces <- function(fit) {

  components <- attr(fit, "components")
  locations  <- if (components > 1L) paste0(.sapmFamily(attr(fit, "family"))[["location"]], seq_len(components)) else fit[["dlist"]][["location"]]
  return(lapply(locations, function(location) {
    design <- fit[["data"]][["mml"]][[location]]
    if (is.null(design))
      design <- matrix(1, nrow(fit[["data"]][["m"]]), 1L, dimnames = list(rownames(fit[["data"]][["m"]]), "(Intercept)"))
    parameters <- c(match(location, rownames(fit[["res.t"]])), fit[["covpars"]][fit[["mx"]][[location]]])
    fixed  <- which(parameters %in% fit[["fixedpars"]])
    offset <- if (length(fixed) > 0L) as.vector(design[, fixed, drop = FALSE] %*% fit[["res.t"]][parameters[fixed], "est"]) else rep(0, nrow(design))
    free   <- setdiff(seq_len(ncol(design)), fixed)
    return(list(design = design[, free, drop = FALSE], offset = offset))
  }))
}

.sapModelsNested <- function(fit0, fit1) {

  if (fit1[["npars"]] <= fit0[["npars"]])
    return(FALSE)

  spaces0 <- .sapRegressionSpaces(fit0)
  spaces1 <- .sapRegressionSpaces(fit1)
  if (length(spaces0) != length(spaces1))
    return(FALSE)

  normalize <- function(design) {
    norm <- sqrt(colSums(design^2))
    norm[norm == 0] <- 1
    return(sweep(design, 2, norm, "/"))
  }
  rank <- function(design) if (ncol(design) == 0L) 0L else qr(normalize(design))[["rank"]]
  for (k in seq_along(spaces0)) {
    if (!identical(rownames(spaces0[[k]][["design"]]), rownames(spaces1[[k]][["design"]])))
      return(FALSE)
    design0    <- spaces0[[k]][["design"]]
    design1    <- spaces1[[k]][["design"]]
    difference <- spaces0[[k]][["offset"]] - spaces1[[k]][["offset"]]
    # Fixed coefficients create an affine offset. Both that offset and every
    # free design column must be representable in the larger model's space.
    if (rank(cbind(design1, design0, difference)) != rank(design1))
      return(FALSE)
  }

  return(TRUE)
}
