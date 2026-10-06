# Fixed native-scale baseline distribution parameters, shared by both analyses.
.sapFixedParameters <- function(options, distribution, components = 1L) {

  fixed <- rep(list(numeric(0)), components)
  if (!options[["restrictParameters"]])
    return(fixed)

  rows   <- options[["fixedParameters"]]
  family <- flexsurv::flexsurv.dists[[distribution]]
  key    <- names(.sapDistributions)[match(distribution, .sapDistributions)]
  keys   <- vapply(rows, function(x) x[["restriction"]][[1L]], character(1))

  for (k in seq_len(components)) {
    values <- rows[[match(paste0(key, ":", k), keys)]]
    for (i in seq_along(family[["pars"]])) {
      parameter  <- family[["pars"]][i]
      expression <- values[[parameter]]
      if (is.character(expression) && length(expression) == 1L && trimws(expression) == "")
        next

      value <- .sapNumericExpression(expression)
      valid <- .sapFiniteNumeric(value, expectedLength = 1L)
      if (valid)
        valid <- is.finite(suppressWarnings(family[["transforms"]][[i]](value)))
      if (!valid)
        stop(gettextf("The fixed value for %1$s in %2$s, component %3$i, must be a single finite number within the parameter's range.",
                      parameter, .sapOption2DistributionName(distribution), k))
      fixed[[k]][parameter] <- value
    }
  }

  return(fixed)
}
.sapCheckFixedConstraint <- function(fixed, constraint) {

  if (is.null(constraint))
    return()
  for (k in seq_along(fixed)) {
    value <- fixed[[k]][constraint[["parameter"]]]
    if (length(value) == 0L || is.na(value))
      next
    feasible <- if (constraint[["direction"]] == "lower") value >= constraint[["naturalBound"]] else value <= constraint[["naturalBound"]]
    if (!feasible)
      stop(gettextf("The fixed %1$s in component %2$i conflicts with the minimum log-time standard deviation. Change the fixed value or the minimum spread.",
                    constraint[["parameter"]], k))
  }

  return()
}
.sapFixedDistribution <- function(dlist, fixed, parameterNames) {

  initialize <- dlist[["inits"]]
  dlist[["inits"]] <- function(t, mf, mml, aux) {
    arguments <- list(t = t, mf = mf, mml = mml, aux = aux)
    inits <- do.call(initialize, arguments[intersect(names(arguments), names(formals(initialize)))])
    if (length(inits) < length(parameterNames))
      inits <- c(inits, rep(0, length(parameterNames) - length(inits)))
    inits[match(names(fixed), parameterNames)] <- fixed
    return(inits)
  }

  return(dlist)
}
.sapFixedFitCall <- function(fitCall, fixed, parameterNames = NULL) {

  if (length(fixed) == 0L)
    return(fitCall)

  dlist <- fitCall[["dist"]]
  if (is.character(dlist))
    dlist <- flexsurv::flexsurv.dists[[dlist]]
  if (is.null(parameterNames))
    parameterNames <- dlist[["pars"]]
  index <- match(names(fixed), parameterNames)
  fitCall[["fixedpars"]] <- sort(index)
  if (!is.null(fitCall[["inits"]]))
    fitCall[["inits"]][index] <- fixed
  else
    dlist <- .sapFixedDistribution(dlist, fixed, parameterNames)
  fitCall[["dist"]] <- dlist

  # flexsurv passes only the free parameters to optim; bounds and steps must match.
  for (name in c("lower", "upper")) {
    if (length(fitCall[[name]]) > 1L)
      fitCall[[name]] <- fitCall[[name]][-index]
  }
  if (length(fitCall[["control"]][["ndeps"]]) > 1L)
    fitCall[["control"]][["ndeps"]] <- fitCall[["control"]][["ndeps"]][-index]

  return(fitCall)
}
.sapCompleteFixedFit <- function(fit, fixed, level) {

  if (length(fixed) == 0L || inherits(fit, "try-error"))
    return(fit)

  fit[["fixedpars"]] <- match(names(fixed), rownames(fit[["res"]]))
  fit[["optpars"]]   <- setdiff(seq_len(nrow(fit[["res"]])), fit[["fixedpars"]])
  attr(fit, "fixedParameters") <- fixed

  # flexsurv's fully fixed fit returns estimates only, with no optimizer or intervals.
  if (fit[["npars"]] == 0L) {
    columns <- c("est", paste0(c("L", "U"), round(100 * level), "%"), "se")
    for (name in c("res", "res.t")) {
      result <- matrix(NA_real_, nrow(fit[[name]]), length(columns), dimnames = list(rownames(fit[[name]]), columns))
      result[, "est"] <- fit[[name]][, "est"]
      fit[[name]] <- result
    }
    fit[["cl"]]           <- level
    fit[["coefficients"]] <- fit[["res.t"]][, "est"]
    fit[["cov"]]          <- matrix(numeric(0), 0L, 0L)
    fit[["opt"]]          <- list(convergence = 0L, par = numeric(0), hessian = matrix(numeric(0), 0L, 0L))
  }

  return(fit)
}
.sapParameterCovariance <- function(fit) {

  parameters <- rownames(fit[["res.t"]])
  covariance <- matrix(0, length(parameters), length(parameters), dimnames = list(parameters, parameters))
  free <- setdiff(seq_along(parameters), fit[["fixedpars"]])
  if (length(free) > 0L)
    covariance[free, free] <- fit[["cov"]]

  return(covariance)
}
.sapFitSingle <- function(dataset, options, distribution, modelTerms, hessian = TRUE) {

  fixed      <- .sapFixedParameters(options, distribution)[[1L]]
  formula    <- .sapGetFormula(options, modelTerms, dataset)
  regression <- .sapRegressionFixed(dataset, modelTerms, formula)
  fixed      <- c(fixed, regression)
  fitCall <- .sapFixedFitCall(list(
    formula = formula,
    data    = dataset,
    dist    = distribution,
    hessian = hessian,
    cl      = options[["coefficientsConfidenceIntervalLevel"]]
  ), fixed, c(flexsurv::flexsurv.dists[[distribution]][["pars"]], attr(regression, "parameterNames")))
  if (options[["weights"]] != "")
    fitCall[["weights"]] <- dataset[[options[["weights"]]]]

  fit <- do.call(flexsurv::flexsurvreg, fitCall)
  return(.sapCompleteFixedFit(fit, fixed, options[["coefficientsConfidenceIntervalLevel"]]))
}
