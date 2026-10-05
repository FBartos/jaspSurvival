# Restrictions use the coefficients displayed for each main-effect term.
.sapParseParameterRestriction <- function(expression, variable, modelTitle) {

  if (is.character(expression) && length(expression) == 1L && trimws(expression) == "")
    return(numeric(0))

  value <- expression
  if (is.character(expression) && length(expression) == 1L) {
    text  <- gsub("[;|\r\n\t]+", ",", trimws(expression))
    value <- try(eval(parse(text = paste0("c(", text, ")")), envir = baseenv()), silent = TRUE)
    if (inherits(value, "try-error")) {
      text  <- paste(strsplit(text, "[[:space:]]+")[[1L]], collapse = ",")
      value <- try(eval(parse(text = paste0("c(", text, ")")), envir = baseenv()), silent = TRUE)
    }
  }

  if (!is.numeric(value) || is.complex(value) || length(value) == 0L || any(!is.finite(value)))
    .quitAnalysis(gettextf("The parameter restriction for '%1$s' in '%2$s' must contain finite numeric values or expressions, separated by commas or semicolons.",
                          jaspBase::decodeColNames(variable), modelTitle))

  return(unname(value))
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
    if (length(regression[[k]]) == 0L) return(numeric(0))
    names <- if (k == 1L) names(regression[[k]]) else paste0(family[["location"]], k, "(", names(regression[[k]]), ")")
    return(stats::setNames(as.numeric(regression[[k]]), names))
  }), use.names = TRUE))
}
.sapmParameterOrder <- function(mixture, family, components, covariates) {

  if (length(covariates) == 0L)
    return(mixture[["dlist"]][["pars"]])

  return(c(mixture[["dlist"]][["pars"]], unlist(lapply(seq_len(components), function(k) {
    if (k == 1L) return(covariates)
    return(paste0(family[["location"]], k, "(", covariates, ")"))
  }), use.names = FALSE)))
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
