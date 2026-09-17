#
# Copyright (C) 2013-2018 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# the parametric mixture survival analysis shares the fitting, selection, and output of the parametric survival analysis
# (see .sapRun), this file contains the mixture estimator and the mixture specific output
.sapmDependencies <- c(
  "mixtureComponents", "mixtureMaximumComponents", "mixtureInitialization", "mixtureRestarts", "mixtureMaximumIterations",
  "setSeed", "seed", "compareModelsAcrossComponents"
)

.sapmCheckDataset               <- function(dataset, options) {

  # the initialization of the EM algorithm clusters log event times
  if (options[["censoringType"]] == "interval") {
    exact <- !is.na(dataset[[options[["intervalStart"]]]]) & !is.na(dataset[[options[["intervalEnd"]]]]) & dataset[[options[["intervalStart"]]]] == dataset[[options[["intervalEnd"]]]]
    time  <- dataset[[options[["intervalStart"]]]][exact]
  } else {
    time  <- dataset[[if (options[["censoringType"]] == "counting") options[["intervalEnd"]] else options[["timeToEvent"]]]]
    time  <- time[dataset[[options[["eventStatus"]]]]]
  }

  if (any(time <= 0))
    .quitAnalysis(gettext("The mixture model requires all event times to be positive."))

  return()
}

# mixture estimator
# the mixture is estimated with an EM algorithm (weighted M-steps for each component) and the converged
# solution is wrapped as a custom flexsurvreg distribution so that all flexsurvreg methods remain available
.sapmFitModel                   <- function(dataset, options, distribution, modelTerms, components) {

  # seeding each model makes the fits independent of the order in which they are fitted
  jaspBase::.setSeedJASP(options)

  fit <- try(.sapmFitMixture(dataset, options, distribution, modelTerms, components))

  return(fit)
}
.sapmFitMixture                 <- function(dataset, options, distribution, modelTerms, components) {

  family      <- .sapmFamily(distribution)
  formula     <- .sapGetFormula(options, modelTerms)
  survObject  <- .saGetSurvObject(options, dataset)
  caseWeights <- if (options[["weights"]] != "") dataset[[options[["weights"]]]] else rep(1, nrow(dataset))

  # the M-steps always estimate the intercept (flexsurvreg ignores its removal as well)
  termLabels  <- attr(stats::terms(formula), "term.labels")
  emFormula   <- stats::reformulate(if (length(termLabels) > 0) termLabels else "1", response = formula[[2]])
  covariates  <- stats::model.matrix(emFormula, stats::model.frame(emFormula, dataset))[, -1, drop = FALSE]

  # EM from the selected start and from additional starts
  initializations <- c(options[["mixtureInitialization"]], .sapmRestartInitializations(options))
  emFits          <- lapply(initializations, function(initialization) try(.sapmEm(
    emFormula         = emFormula,
    dataset           = dataset,
    survObject        = survObject,
    covariates        = covariates,
    family            = family,
    components        = components,
    caseWeights       = caseWeights,
    initialization    = initialization,
    maximumIterations = options[["mixtureMaximumIterations"]]
  ), silent = TRUE))
  emValid <- !vapply(emFits, jaspBase::isTryError, logical(1))
  if (!any(emValid))
    stop(gettextf("The EM algorithm failed: %1$s", .sapmCleanError(emFits[[1]])))
  em      <- emFits[emValid][[which.max(vapply(emFits[emValid], function(x) x[["logLik"]], numeric(1)))]]

  # wrap the solution with components ordered by their baseline median lifetime
  wrap <- .sapmWrap(formula, dataset, options, family, components, em, caseWeights)
  if (jaspBase::isTryError(wrap[["fit"]]))
    stop(gettextf("The mixture model could not be finalized: %1$s", .sapmCleanError(wrap[["fit"]])))

  fit <- wrap[["fit"]]
  if (!is.finite(fit[["loglik"]]))
    stop(gettext("The log-likelihood of the mixture model is not finite."))

  # the constructed call contains the data and the distribution functions
  fit[["call"]] <- NULL
  estimates     <- .sapmComponentEstimates(fit, family, components)

  attr(fit, "mixture") <- list(
    family         = family[["family"]],
    components     = components,
    emIterations   = em[["iterations"]],
    emConverged    = em[["converged"]],
    emLogLik       = em[["logLik"]],
    emGap          = abs(fit[["loglik"]] - em[["logLik"]]) > 1e-3,
    restarts       = length(initializations) - 1,
    hessianWarning = wrap[["hessianWarning"]] || any(!is.finite(fit[["cov"]])),
    collapsed      = which(estimates[["probabilities"]] < 1e-3 | estimates[["collapsed"]]),
    posterior      = .sapmPosterior(fit, family, components, survObject)
  )

  return(fit)
}
.sapmRestartInitializations     <- function(options) {

  if (options[["mixtureRestarts"]] == 0)
    return(character(0))

  # the first restart uses a deterministic alternative start, the remaining restarts use random partitions
  firstRestart <- if (options[["mixtureInitialization"]] == "quantiles") "kmeans" else "quantiles"

  return(c(firstRestart, rep("random", options[["mixtureRestarts"]] - 1)))
}
.sapmEm                         <- function(emFormula, dataset, survObject, covariates, family, components, caseWeights, initialization, maximumIterations) {

  nObs          <- nrow(survObject)
  posterior     <- .sapmInitialPosterior(survObject, components, initialization)
  probabilities <- colSums(posterior * caseWeights) / sum(caseWeights)
  mSteps        <- NULL
  logLik        <- -Inf
  iterations    <- 0
  converged     <- FALSE

  for (iteration in seq_len(maximumIterations)) {

    # M-step: weighted fit of each component
    newMSteps <- lapply(seq_len(components), function(k) .sapmMStep(
      emFormula  = emFormula,
      dataset    = dataset,
      survObject = survObject,
      covariates = covariates,
      family     = family,
      weights    = posterior[, k] * caseWeights,
      previous   = mSteps[[k]]
    ))

    # E-step: posterior probabilities of component membership
    likelihood <- vapply(newMSteps, function(mStep) .sapmComponentLikelihood(family, survObject, mStep[["parameters"]]), numeric(nObs))
    likelihood <- matrix(pmax(likelihood, .Machine$double.xmin), nrow = nObs)
    joint      <- sweep(likelihood, 2, probabilities, "*")
    marginal   <- rowSums(joint)
    newLogLik  <- sum(caseWeights * log(marginal)) - .sapmTruncationLogLik(family, survObject, newMSteps, probabilities, caseWeights)

    if (!is.finite(newLogLik))
      stop(gettext("The log-likelihood is not finite."))

    # numerical safeguard: the EM iterations cannot decrease the likelihood
    if (iteration > 1 && newLogLik < logLik - 1e-6 * abs(logLik))
      break

    mSteps        <- newMSteps
    posterior     <- joint / marginal
    probabilities <- colSums(posterior * caseWeights) / sum(caseWeights)
    iterations    <- iteration

    if (is.finite(logLik) && abs(newLogLik - logLik) < 1e-8 * abs(newLogLik)) {
      logLik    <- newLogLik
      converged <- TRUE
      break
    }
    logLik <- newLogLik
  }

  return(list(
    logLik        = logLik,
    iterations    = iterations,
    converged     = converged,
    probabilities = probabilities,
    mSteps        = mSteps
  ))
}
.sapmMStep                      <- function(emFormula, dataset, survObject, covariates, family, weights, previous) {

  # posterior probabilities can underflow to zero which is not allowed as a weight
  weights <- pmax(weights, 1e-10)

  if (!is.null(family[["survreg"]]) && attr(survObject, "type") != "counting") {

    # survreg evaluates weights non-standardly, the call needs to be constructed
    fit <- suppressWarnings(do.call(survival::survreg, list(
      formula = emFormula,
      data    = dataset,
      weights = weights,
      dist    = family[["survreg"]]
    )))

    coefficients <- stats::coef(fit)
    base         <- switch(
      family[["family"]],
      "exp"     = c(rate    = exp(-coefficients[[1]])),
      "lnorm"   = c(meanlog = coefficients[[1]],  sdlog = fit[["scale"]]),
      "llogis"  = c(shape   = 1 / fit[["scale"]], scale = exp(coefficients[[1]])),
      "weibull" = c(shape   = 1 / fit[["scale"]], scale = exp(coefficients[[1]]))
    )
    # survreg models the log time whereas the exponential rate is its reciprocal
    beta         <- coefficients[-1] * if (family[["family"]] == "exp") -1 else 1

  } else {

    fitCall <- list(
      formula = emFormula,
      data    = dataset,
      weights = weights,
      dist    = family[["family"]]
    )
    # start from the previous M-step estimates
    if (!is.null(previous))
      fitCall[["inits"]] <- c(previous[["base"]], previous[["beta"]])

    fit  <- suppressWarnings(suppressMessages(do.call(flexsurv::flexsurvreg, fitCall)))
    base <- fit[["res"]][family[["pars"]], "est"]
    beta <- fit[["res"]][fit[["covpars"]], "est"]
  }

  if (any(!is.finite(base)) || any(!is.finite(beta)))
    stop(gettext("A mixture component degenerated during the estimation. Consider fewer components."))

  return(list(
    base       = base,
    beta       = beta,
    parameters = .sapmParameters(family, base, beta, covariates)
  ))
}
.sapmParameters                 <- function(family, base, beta, covariates) {

  # natural parameters of a component for each observation (covariates act on the transformed location parameter)
  parameters <- as.list(base)

  if (length(beta) > 0) {
    location <- which(family[["pars"]] == family[["location"]])
    parameters[[location]] <- family[["inv.transforms"]][[location]](
      family[["transforms"]][[location]](base[[location]]) + as.vector(covariates %*% beta)
    )
  }

  return(parameters)
}
.sapmComponentLikelihood        <- function(family, survObject, parameters, at = "observed") {

  subsetParameters <- function(index) lapply(parameters, function(x) if (length(x) > 1) x[index] else x)
  density          <- function(x, index) do.call(family[["d"]], c(list(x), subsetParameters(index)))
  distribution     <- function(x, index, lowerTail = TRUE) do.call(family[["p"]], c(list(x), subsetParameters(index), list(lower.tail = lowerTail)))

  type <- attr(survObject, "type")
  out  <- numeric(nrow(survObject))

  if (at == "start") {
    # survival at the left-truncation time
    index       <- survObject[, "start"] > 0
    out[!index] <- 1
    if (any(index))
      out[index] <- distribution(survObject[index, "start"], index, lowerTail = FALSE)
    return(out)
  }

  if (type %in% c("right", "counting")) {

    time  <- survObject[, if (type == "right") "time" else "stop"]
    event <- survObject[, "status"] == 1

    if (any(event))
      out[event]  <- density(time[event], event)
    if (any(!event))
      out[!event] <- distribution(time[!event], !event, lowerTail = FALSE)

  } else if (type == "interval") {

    # status: 0 = right censored, 1 = exact, 2 = left censored, 3 = interval censored
    status <- survObject[, "status"]
    time1  <- survObject[, "time1"]
    time2  <- survObject[, "time2"]

    for (s in 0:3) {
      index <- status == s
      if (!any(index))
        next
      out[index] <- switch(
        as.character(s),
        "0" = distribution(time1[index], index, lowerTail = FALSE),
        "1" = density(time1[index], index),
        "2" = distribution(time1[index], index),
        "3" = distribution(time2[index], index) - distribution(time1[index], index)
      )
    }

  } else {
    stop(gettextf("Censoring type '%1$s' is not supported by the mixture model.", type))
  }

  return(out)
}
.sapmTruncationLogLik           <- function(family, survObject, mSteps, probabilities, caseWeights) {

  if (attr(survObject, "type") != "counting" || all(survObject[, "start"] <= 0))
    return(0)

  survival <- vapply(mSteps, function(mStep) .sapmComponentLikelihood(family, survObject, mStep[["parameters"]], at = "start"), numeric(nrow(survObject)))
  survival <- matrix(survival, nrow = nrow(survObject))

  return(sum(caseWeights * log(pmax(as.vector(survival %*% probabilities), .Machine$double.xmin))))
}
.sapmInitialPosterior           <- function(survObject, components, initialization) {

  type <- attr(survObject, "type")
  if (type == "right") {
    time <- survObject[, "time"]
  } else if (type == "counting") {
    time <- survObject[, "stop"]
  } else {
    # interval censored observations are clustered at their midpoints and left censored observations at half of their upper limit
    status <- survObject[, "status"]
    time   <- survObject[, "time1"]
    time[status == 3] <- (survObject[status == 3, "time1"] + survObject[status == 3, "time2"]) / 2
    time[status == 2] <- survObject[status == 2, "time1"] / 2
  }

  logTime <- log(pmax(time, min(time[time > 0])))
  nObs    <- length(logTime)

  if (initialization == "kmeans" && length(unique(logTime)) <= components)
    initialization <- "quantiles"

  membership <- switch(
    initialization,
    "kmeans"    = {
      clusters <- stats::kmeans(logTime, centers = components, nstart = 20)
      match(clusters[["cluster"]], order(clusters[["centers"]][, 1]))
    },
    "quantiles" = as.integer(cut(rank(logTime, ties.method = "first"), components, labels = FALSE)),
    "random"    = sample.int(components, nObs, replace = TRUE)
  )

  # soft start avoids empty components
  posterior <- matrix(0.05 / (components - 1), nObs, components)
  posterior[cbind(seq_len(nObs), membership)] <- 0.95

  return(posterior)
}
.sapmWrap                       <- function(formula, dataset, options, family, components, em, caseWeights) {

  mixture <- .sapmMixtureDistribution(family, components)
  base    <- lapply(em[["mSteps"]], function(mStep) mStep[["base"]])
  beta    <- lapply(em[["mSteps"]], function(mStep) mStep[["beta"]])
  order   <- .sapmComponentOrder(family, base)

  wrapFit <- function(inits, method = "BFGS") {

    hessianWarning <- FALSE
    fitCall        <- list(
      formula = formula,
      data    = dataset,
      dist    = mixture[["dlist"]],
      dfns    = mixture[["dfns"]],
      inits   = inits,
      method  = method,
      cl      = options[["coefficientsConfidenceIntervalLevel"]]
    )
    # the covariates enter the location parameter of each component
    if (length(beta[[1]]) > 0)
      fitCall[["anc"]] <- stats::setNames(rep(list(formula[-2]), components - 1), paste0(family[["location"]], 2:components))
    if (options[["weights"]] != "")
      fitCall[["weights"]] <- caseWeights

    fit <- try(withCallingHandlers(
      suppressMessages(do.call(flexsurv::flexsurvreg, fitCall)),
      warning = function(w) {
        if (grepl("Hessian", conditionMessage(w), fixed = TRUE))
          hessianWarning <<- TRUE
        invokeRestart("muffleWarning")
      }
    ), silent = TRUE)

    if (!jaspBase::isTryError(fit) && fit[["opt"]][["convergence"]] != 0)
      fit <- try(stop(gettext("The optimizer did not converge.")), silent = TRUE)

    return(list(fit = fit, hessianWarning = hessianWarning))
  }

  inits <- .sapmInits(mixture, base[order], beta[order], em[["probabilities"]][order])
  wrap  <- wrapFit(inits)
  if (jaspBase::isTryError(wrap[["fit"]]))
    wrap <- wrapFit(inits, method = "Nelder-Mead")
  if (jaspBase::isTryError(wrap[["fit"]]))
    return(wrap)

  # the direct maximization might change the ordering of the components
  estimates <- .sapmComponentEstimates(wrap[["fit"]], family, components)
  newOrder  <- .sapmComponentOrder(family, estimates[["base"]])
  if (!identical(newOrder, seq_len(components))) {
    reordered <- wrapFit(.sapmInits(mixture, estimates[["base"]][newOrder], estimates[["beta"]][newOrder], estimates[["probabilities"]][newOrder]))
    if (!jaspBase::isTryError(reordered[["fit"]]))
      wrap <- reordered
  }

  return(wrap)
}
.sapmInits                      <- function(mixture, base, beta, probabilities) {

  # flexsurvreg order: component parameters, stick-breaking weights, location effects of the first component, location effects of the remaining components
  components <- length(base)
  weights    <- numeric(0)
  remaining  <- 1
  for (k in seq_len(components - 1)) {
    weights   <- c(weights, min(max(probabilities[k] / remaining, 1e-6), 1 - 1e-6))
    remaining <- remaining - probabilities[k]
  }

  return(c(
    stats::setNames(unlist(base, use.names = FALSE), mixture[["componentPars"]]),
    stats::setNames(weights, mixture[["weightPars"]]),
    unlist(beta, use.names = FALSE)
  ))
}
.sapmComponentOrder             <- function(family, base) {

  # ascending baseline median lifetime, ties are broken by the remaining parameters
  medians <- vapply(base, function(x) .sapmComponentMedian(family, x), numeric(1))
  ties    <- lapply(seq_along(family[["pars"]]), function(i) vapply(base, function(x) x[[i]], numeric(1)))

  return(do.call(order, c(list(medians), ties)))
}
.sapmComponentMedian            <- function(family, base) {
  return(do.call(family[["q"]], c(list(0.5), as.list(base))))
}
.sapmBaseParameters             <- function(estimates, family, components) {

  # natural baseline parameters (covariates at zero) and mixing probabilities from the transformed estimates
  base    <- lapply(seq_len(components), function(k) {
    stats::setNames(vapply(seq_along(family[["pars"]]), function(i) {
      family[["inv.transforms"]][[i]](estimates[[paste0(family[["pars"]][i], k)]])
    }, numeric(1)), family[["pars"]])
  })
  weights <- if (components > 1) stats::plogis(estimates[paste0("v", seq_len(components - 1))])

  return(list(
    base          = base,
    probabilities = as.vector(.sapmStickBreaking(if (components > 1) matrix(weights, nrow = 1), 1))
  ))
}
.sapmComponentEstimates         <- function(fit, family, components) {

  estimates  <- fit[["res.t"]][, "est"]
  parameters <- .sapmBaseParameters(estimates, family, components)
  isLog      <- vapply(family[["transforms"]], function(f) identical(f, log), logical(1))

  return(list(
    base          = parameters[["base"]],
    probabilities = parameters[["probabilities"]],
    beta          = lapply(seq_len(components), function(k) {
      index <- fit[["mx"]][[paste0(family[["location"]], k)]]
      if (length(index) == 0)
        return(numeric(0))
      return(estimates[fit[["covpars"]][index]])
    }),
    # a component collapses if its log-scale parameters diverge
    collapsed     = vapply(seq_len(components), function(k) {
      any(abs(estimates[paste0(family[["pars"]], k)][isLog]) > 15)
    }, logical(1))
  ))
}
.sapmPosterior                  <- function(fit, family, components, survObject) {

  estimates  <- .sapmComponentEstimates(fit, family, components)
  covariates <- fit[["data"]][["mml"]][[paste0(family[["location"]], 1)]]
  if (!is.null(covariates))
    covariates <- covariates[, -1, drop = FALSE]

  likelihood <- vapply(seq_len(components), function(k) {
    .sapmComponentLikelihood(family, survObject, .sapmParameters(family, estimates[["base"]][[k]], estimates[["beta"]][[k]], covariates))
  }, numeric(nrow(survObject)))
  likelihood <- matrix(pmax(likelihood, .Machine$double.xmin), nrow = nrow(survObject))
  joint      <- sweep(likelihood, 2, estimates[["probabilities"]], "*")

  return(joint / rowSums(joint))
}
.sapmStickBreaking              <- function(weights, nObs) {

  # stick-breaking weights v1, ..., v(K-1) to mixing probabilities p1, ..., pK
  if (is.null(weights) || length(weights) == 0)
    return(matrix(1, nObs, 1))

  weights       <- as.matrix(weights)
  probabilities <- matrix(NA_real_, nrow(weights), ncol(weights) + 1)
  remaining     <- rep(1, nrow(weights))
  for (k in seq_len(ncol(weights))) {
    probabilities[, k] <- remaining * weights[, k]
    remaining          <- remaining - probabilities[, k]
  }
  probabilities[, ncol(weights) + 1] <- remaining

  return(probabilities)
}
.sapmMixtureDistribution        <- function(family, components) {

  componentPars <- as.vector(outer(family[["pars"]], seq_len(components), paste0))
  weightPars    <- if (components > 1) paste0("v", seq_len(components - 1)) else character(0)

  # flexsurv passes the parameters as scalars without covariates and as vectors with covariates
  splitArguments <- function(arguments, n) {
    n         <- max(n, lengths(arguments))
    arguments <- lapply(arguments, rep_len, length.out = n)
    return(list(
      n             = n,
      components    = lapply(seq_len(components), function(k) stats::setNames(arguments[paste0(family[["pars"]], k)], family[["pars"]])),
      probabilities = .sapmStickBreaking(if (components > 1) do.call(cbind, arguments[weightPars]), n)
    ))
  }
  logSumExp      <- function(x) {
    maximum <- Reduce(pmax, lapply(seq_len(ncol(x)), function(k) x[, k]))
    out     <- maximum + log(rowSums(exp(x - maximum)))
    out[!is.finite(maximum)] <- maximum[!is.finite(maximum)]
    return(out)
  }

  dMixture <- function(x, ..., log = FALSE) {
    arguments <- splitArguments(list(...), length(x))
    x         <- rep_len(x, arguments[["n"]])
    logs      <- vapply(seq_len(components), function(k) {
      log(arguments[["probabilities"]][, k]) + do.call(family[["d"]], c(list(x), arguments[["components"]][[k]], list(log = TRUE)))
    }, numeric(arguments[["n"]]))
    out       <- logSumExp(matrix(logs, nrow = arguments[["n"]]))
    return(if (log) out else exp(out))
  }
  pMixture <- function(q, ..., lower.tail = TRUE, log.p = FALSE) {
    arguments <- splitArguments(list(...), length(q))
    q         <- rep_len(q, arguments[["n"]])
    if (!log.p) {
      return(Reduce(`+`, lapply(seq_len(components), function(k) {
        arguments[["probabilities"]][, k] * do.call(family[["p"]], c(list(q), arguments[["components"]][[k]], list(lower.tail = lower.tail)))
      })))
    }
    logs      <- vapply(seq_len(components), function(k) {
      log(arguments[["probabilities"]][, k]) + do.call(family[["p"]], c(list(q), arguments[["components"]][[k]], list(lower.tail = lower.tail, log.p = TRUE)))
    }, numeric(arguments[["n"]]))
    return(logSumExp(matrix(logs, nrow = arguments[["n"]])))
  }
  qMixture <- function(p, ..., lower.tail = TRUE, log.p = FALSE) {
    if (log.p)
      p <- exp(p)
    if (!lower.tail)
      p <- 1 - p
    arguments <- splitArguments(list(...), length(p))
    p         <- rep_len(p, arguments[["n"]])
    # the mixture quantile lies between the smallest and the largest component quantile (found by bisection)
    bounds    <- matrix(vapply(seq_len(components), function(k) {
      do.call(family[["q"]], c(list(p), arguments[["components"]][[k]]))
    }, numeric(arguments[["n"]])), nrow = arguments[["n"]])
    lower     <- Reduce(pmin, lapply(seq_len(components), function(k) bounds[, k]))
    upper     <- Reduce(pmax, lapply(seq_len(components), function(k) bounds[, k]))
    for (i in seq_len(100)) {
      middle <- (lower + upper) / 2
      below  <- do.call(pMixture, c(list(middle), list(...))) < p
      lower  <- ifelse(below, middle, lower)
      upper  <- ifelse(below, upper, middle)
    }
    out <- (lower + upper) / 2
    out[p <= 0] <- 0
    out[p >= 1] <- Inf
    return(out)
  }
  # the restricted mean survival time and the mean are mixtures of the component quantities (without left-truncation)
  rmstMixture <- function(t, ..., start = 0) {
    if (any(start > 0))
      return(flexsurv::rmst_generic(pMixture, t = t, start = start, ...))
    arguments <- splitArguments(list(...), length(t))
    t         <- rep_len(t, arguments[["n"]])
    return(Reduce(`+`, lapply(seq_len(components), function(k) {
      arguments[["probabilities"]][, k] * do.call(family[["rmst"]], c(list(t), arguments[["components"]][[k]]))
    })))
  }
  meanMixture <- function(..., start = 0) {
    if (any(start > 0))
      return(flexsurv::rmst_generic(pMixture, t = Inf, start = start, ...))
    arguments <- splitArguments(list(...), 1)
    return(Reduce(`+`, lapply(seq_len(components), function(k) {
      arguments[["probabilities"]][, k] * do.call(family[["mean"]], arguments[["components"]][[k]])
    })))
  }
  rMixture <- function(n, ...) {
    arguments  <- splitArguments(list(...), n)
    membership <- vapply(seq_len(n), function(i) sample.int(components, 1, prob = arguments[["probabilities"]][i, ]), integer(1))
    uniform    <- stats::runif(n)
    out        <- numeric(n)
    for (k in seq_len(components)) {
      index <- membership == k
      if (any(index))
        out[index] <- do.call(family[["q"]], c(list(uniform[index]), lapply(arguments[["components"]][[k]], function(x) x[index])))
    }
    return(out)
  }

  return(list(
    dlist         = list(
      name           = paste0("mixture.", family[["family"]], ".", components),
      pars           = c(componentPars, weightPars),
      location       = paste0(family[["location"]], 1),
      transforms     = c(rep(family[["transforms"]], components),     rep(list(stats::qlogis), components - 1)),
      inv.transforms = c(rep(family[["inv.transforms"]], components), rep(list(stats::plogis), components - 1))
    ),
    dfns          = list(d = dMixture, p = pMixture, q = qMixture, r = rMixture, rmst = rmstMixture, mean = meanMixture),
    componentPars = componentPars,
    weightPars    = weightPars
  ))
}
.sapmFamily                     <- function(distribution) {

  specification <- flexsurv::flexsurv.dists[[distribution]]

  return(list(
    family         = distribution,
    pars           = specification[["pars"]],
    location       = specification[["location"]],
    transforms     = specification[["transforms"]],
    inv.transforms = specification[["inv.transforms"]],
    # families with a survreg equivalent use the faster weighted survreg M-step
    survreg        = switch(distribution, "exp" = "exponential", "lnorm" = "lognormal", "llogis" = "loglogistic", "weibull" = "weibull", NULL),
    d              = switch(
      distribution,
      "exp"           = stats::dexp,
      "gamma"         = stats::dgamma,
      "genf"          = flexsurv::dgenf,
      "gengamma"      = flexsurv::dgengamma,
      "gompertz"      = flexsurv::dgompertz,
      "llogis"        = flexsurv::dllogis,
      "lnorm"         = stats::dlnorm,
      "weibull"       = stats::dweibull,
      "gengamma.orig" = flexsurv::dgengamma.orig,
      "genf.orig"     = flexsurv::dgenf.orig
    ),
    p              = switch(
      distribution,
      "exp"           = stats::pexp,
      "gamma"         = stats::pgamma,
      "genf"          = flexsurv::pgenf,
      "gengamma"      = flexsurv::pgengamma,
      "gompertz"      = flexsurv::pgompertz,
      "llogis"        = flexsurv::pllogis,
      "lnorm"         = stats::plnorm,
      "weibull"       = stats::pweibull,
      "gengamma.orig" = flexsurv::pgengamma.orig,
      "genf.orig"     = flexsurv::pgenf.orig
    ),
    q              = switch(
      distribution,
      "exp"           = stats::qexp,
      "gamma"         = stats::qgamma,
      "genf"          = flexsurv::qgenf,
      "gengamma"      = flexsurv::qgengamma,
      "gompertz"      = flexsurv::qgompertz,
      "llogis"        = flexsurv::qllogis,
      "lnorm"         = stats::qlnorm,
      "weibull"       = stats::qweibull,
      "gengamma.orig" = flexsurv::qgengamma.orig,
      "genf.orig"     = flexsurv::qgenf.orig
    ),
    h              = switch(
      distribution,
      "exp"           = flexsurv::hexp,
      "gamma"         = flexsurv::hgamma,
      "genf"          = flexsurv::hgenf,
      "gengamma"      = flexsurv::hgengamma,
      "gompertz"      = flexsurv::hgompertz,
      "llogis"        = flexsurv::hllogis,
      "lnorm"         = flexsurv::hlnorm,
      "weibull"       = flexsurv::hweibull,
      "gengamma.orig" = flexsurv::hgengamma.orig,
      "genf.orig"     = flexsurv::hgenf.orig
    ),
    rmst           = switch(
      distribution,
      "exp"           = flexsurv::rmst_exp,
      "gamma"         = flexsurv::rmst_gamma,
      "genf"          = flexsurv::rmst_genf,
      "gengamma"      = flexsurv::rmst_gengamma,
      "gompertz"      = flexsurv::rmst_gompertz,
      "llogis"        = flexsurv::rmst_llogis,
      "lnorm"         = flexsurv::rmst_lnorm,
      "weibull"       = flexsurv::rmst_weibull,
      "gengamma.orig" = flexsurv::rmst_gengamma.orig,
      "genf.orig"     = flexsurv::rmst_genf.orig
    ),
    mean           = switch(
      distribution,
      "exp"           = flexsurv::mean_exp,
      "gamma"         = flexsurv::mean_gamma,
      "genf"          = flexsurv::mean_genf,
      "gengamma"      = flexsurv::mean_gengamma,
      "gompertz"      = flexsurv::mean_gompertz,
      "llogis"        = flexsurv::mean_llogis,
      "lnorm"         = flexsurv::mean_lnorm,
      "weibull"       = flexsurv::mean_weibull,
      "gengamma.orig" = flexsurv::mean_gengamma.orig,
      "genf.orig"     = flexsurv::mean_genf.orig
    )
  ))
}
.sapmCleanError                 <- function(error) {
  return(conditionMessage(attr(error, "condition")))
}

# mixture output
.sapmComponentsTable            <- function(jaspResults, options) {

  if (!is.null(jaspResults[["mixtureComponentsTable"]]))
    return()

  # the extract function automatically groups models by subgroup / distribution / components
  fit <- .sapExtractFit(jaspResults, options, type = "selected")
  fit <- .sapmFilterMixtures(.sapFlattenFit(fit, options), options)
  if (.saSurvivalReady(options) && length(fit) == 0)
    return()

  # output dependencies
  outputDependencies <- c(.sapGetDependencies(options), "compareModelsAcrossDistributions", "interpretModel", "alwaysDisplayModelInformation",
                          "mixtureComponentsTable", "coefficientsConfidenceInterval", "coefficientsConfidenceIntervalLevel")

  .sapSectionWrapper(
    jaspResults   = jaspResults,
    options       = options,
    fit           = fit,
    tableFunction = .sapmComponentsTableFun,
    name          = "mixtureComponentsTable",
    title         = gettext("Mixture Components"),
    dependencies  = outputDependencies,
    position      = 2.2
  )

  return()
}
.sapmClassificationTable        <- function(jaspResults, options) {

  if (!is.null(jaspResults[["mixtureClassificationTable"]]))
    return()

  fit <- .sapExtractFit(jaspResults, options, type = "selected")
  fit <- .sapmFilterMixtures(.sapFlattenFit(fit, options), options)
  if (.saSurvivalReady(options) && length(fit) == 0)
    return()

  outputDependencies <- c(.sapGetDependencies(options), "compareModelsAcrossDistributions", "interpretModel", "alwaysDisplayModelInformation",
                          "mixtureClassificationTable")

  .sapSectionWrapper(
    jaspResults   = jaspResults,
    options       = options,
    fit           = fit,
    tableFunction = .sapmClassificationTableFun,
    name          = "mixtureClassificationTable",
    title         = gettext("Mixture Classification"),
    dependencies  = outputDependencies,
    position      = 2.3
  )

  return()
}
.sapmComponentPlot              <- function(jaspResults, options) {

  if (!is.null(jaspResults[["mixtureComponentPlot"]]))
    return()

  fit <- .sapExtractFit(jaspResults, options, type = "selected")
  fit <- .sapmFilterMixtures(.sapFlattenFit(fit, options), options)
  if (.saSurvivalReady(options) && length(fit) == 0)
    return()
  fit <- .sapNestFit(fit)

  outputDependencies <- c(.sapGetDependencies(options), "compareModelsAcrossDistributions", "interpretModel", "alwaysDisplayModelInformation",
                          "mixtureComponentPlot", "mixtureComponentPlotType", "mixtureComponentPlotKaplanMeier",
                          "predictionsConfidenceInterval", "predictionsConfidenceIntervalLevel",
                          "predictionsLifeTimeStepsType", "predictionsLifeTimeStepsNumber", "predictionsLifeTimeStepsFrom", "predictionsLifeTimeStepsSize",
                          "predictionsLifeTimeStepsTo", "predictionsLifeTimeCustom",
                          "colorPalette", "plotLegend", "plotTheme"
  )

  .sapSectionWrapper(
    jaspResults   = jaspResults,
    options       = options,
    fit           = fit,
    tableFunction = .sapmComponentPlotFun,
    name          = "mixtureComponentPlot",
    title         = gettext("Mixture Components"),
    dependencies  = outputDependencies,
    position      = 3.5
  )

  return()
}
.sapmFilterMixtures             <- function(fit, options) {

  # mixture specific output is shown only for models with multiple components
  if (!.saSurvivalReady(options))
    return(fit)

  return(Filter(function(x) attr(x, "components") > 1, fit))
}

.sapmComponentsTableFun         <- function(fit, options) {

  # create the table
  componentsTable <- createJaspTable()
  .sapAddColumnSubgroup(     componentsTable, options, output = "coefficientsCovarianceMatrix")
  .sapAddColumnDistribution( componentsTable, options, output = "coefficientsCovarianceMatrix")
  .sapAddColumnComponents(   componentsTable, options, output = "coefficientsCovarianceMatrix")
  .sapAddColumnModel(        componentsTable, options, output = "coefficientsCovarianceMatrix")
  componentsTable$addColumnInfo(name = "component", title = gettext("Component"),      type = "string")
  componentsTable$addColumnInfo(name = "quantity",  title = "",                        type = "string")
  componentsTable$addColumnInfo(name = "est",       title = gettext("Estimate"),       type = "number")
  componentsTable$addColumnInfo(name = "se",        title = gettext("Standard Error"), type = "number")
  if (options[["coefficientsConfidenceInterval"]]) {
    overtitleCi <- gettextf("%s%% CI", 100 * options[["coefficientsConfidenceIntervalLevel"]])
    componentsTable$addColumnInfo(name = "lower", title = gettext("Lower"), type = "number", overtitle = overtitleCi)
    componentsTable$addColumnInfo(name = "upper", title = gettext("Upper"), type = "number", overtitle = overtitleCi)
  }

  if (!.saSurvivalReady(options) || jaspBase::isTryError(fit))
    return(componentsTable)

  data <- .sapmComponentsTableData(fit, options[["coefficientsConfidenceIntervalLevel"]])

  data$subgroup        <- NA
  data$distribution    <- NA
  data$components      <- NA
  data$model           <- NA
  data$subgroup[1]     <- attr(fit, "subgroup")
  data$distribution[1] <- attr(fit, "distribution")
  data$components[1]   <- attr(fit, "components")
  data$model[1]        <- attr(fit, "modelTitle")

  # add footnotes
  if (!is.null(attr(fit, "label")) && attr(fit, "label") != "")
    componentsTable$addFootnote(attr(fit, "label"))
  componentsTable$addFootnote(gettext("Standard errors and confidence intervals are based on the delta method."))
  if (length(fit[["covpars"]]) > 0)
    componentsTable$addFootnote(gettext("The component parameters, means, and medians correspond to the reference level of factors and zero value of covariates."))
  for (message in .sapmFitMessages(fit, options))
    componentsTable$addFootnote(message, symbol = gettext("Warning:"))

  componentsTable$setData(data)
  componentsTable$showSpecifiedColumnsOnly <- TRUE

  return(componentsTable)
}
.sapmComponentsTableData        <- function(fit, level) {

  mixture    <- attr(fit, "mixture")
  family     <- .sapmFamily(mixture[["family"]])
  components <- mixture[["components"]]
  estimates  <- fit[["res.t"]][, "est"]

  # quantities of each component: mixing probability, parameters, mean, and median
  quantities <- function(estimates) {
    parameters <- .sapmBaseParameters(estimates, family, components)
    unlist(lapply(seq_len(components), function(k) {
      mean <- try(do.call(family[["mean"]], as.list(parameters[["base"]][[k]])), silent = TRUE)
      c(
        parameters[["probabilities"]][k],
        parameters[["base"]][[k]],
        if (jaspBase::isTryError(mean)) NA else mean,
        .sapmComponentMedian(family, parameters[["base"]][[k]])
      )
    }), use.names = FALSE)
  }

  # delta method with a numerical Jacobian (covariate effects do not affect the baseline quantities)
  estimate  <- quantities(estimates)
  jacobian  <- matrix(0, length(estimate), length(estimates))
  baseIndex <- setdiff(seq_along(estimates), fit[["covpars"]])
  for (i in baseIndex) {
    step          <- 1e-5 * max(abs(estimates[i]), 1)
    upper         <- lower <- estimates
    upper[i]      <- upper[i] + step
    lower[i]      <- lower[i] - step
    jacobian[, i] <- (quantities(upper) - quantities(lower)) / (2 * step)
  }
  standardError <- sqrt(pmax(diag(jacobian %*% fit[["cov"]] %*% t(jacobian)), 0))

  estimate[!is.finite(estimate)]           <- NA
  standardError[!is.finite(standardError)] <- NA

  # confidence intervals are computed on the link scale of each quantity (logit for probabilities, log for positive quantities)
  logParameter <- vapply(family[["transforms"]], function(f) identical(f, log), logical(1))
  links        <- rep(c("logit", ifelse(logParameter, "log", "identity"), "log", "log"), components)
  z            <- stats::qnorm((1 + level) / 2)
  lower        <- upper <- rep(NA_real_, length(estimate))
  logitLink    <- links == "logit"    & !is.na(standardError)
  logLink      <- links == "log"      & !is.na(standardError) & estimate > 0
  identityLink <- links == "identity" & !is.na(standardError)

  lower[logitLink]    <- stats::plogis(stats::qlogis(estimate[logitLink]) - z * standardError[logitLink] / (estimate[logitLink] * (1 - estimate[logitLink])))
  upper[logitLink]    <- stats::plogis(stats::qlogis(estimate[logitLink]) + z * standardError[logitLink] / (estimate[logitLink] * (1 - estimate[logitLink])))
  lower[logLink]      <- exp(log(estimate[logLink]) - z * standardError[logLink] / estimate[logLink])
  upper[logLink]      <- exp(log(estimate[logLink]) + z * standardError[logLink] / estimate[logLink])
  lower[identityLink] <- estimate[identityLink] - z * standardError[identityLink]
  upper[identityLink] <- estimate[identityLink] + z * standardError[identityLink]
  lower[!is.finite(lower)] <- NA
  upper[!is.finite(upper)] <- NA

  nQuantities <- length(family[["pars"]]) + 3
  component   <- rep(NA_character_, length(estimate))
  component[seq(1, length(estimate), by = nQuantities)] <- gettextf("Component %1$i", seq_len(components))

  return(data.frame(
    component = component,
    quantity  = rep(c(gettext("Mixing probability"), family[["pars"]], gettext("Mean"), gettext("Median")), components),
    est       = estimate,
    se        = standardError,
    lower     = lower,
    upper     = upper
  ))
}
.sapmClassificationTableFun     <- function(fit, options) {

  # create the table
  classificationTable <- createJaspTable()
  .sapAddColumnSubgroup(     classificationTable, options, output = "coefficientsCovarianceMatrix")
  .sapAddColumnDistribution( classificationTable, options, output = "coefficientsCovarianceMatrix")
  .sapAddColumnComponents(   classificationTable, options, output = "coefficientsCovarianceMatrix")
  .sapAddColumnModel(        classificationTable, options, output = "coefficientsCovarianceMatrix")
  classificationTable$addColumnInfo(name = "component",     title = gettext("Component"),                   type = "string")
  classificationTable$addColumnInfo(name = "probability",   title = gettext("Mixing probability"),          type = "number")
  classificationTable$addColumnInfo(name = "count",         title = gettext("Count"),                       type = "integer",  overtitle = gettext("Classified"))
  classificationTable$addColumnInfo(name = "proportion",    title = gettext("Proportion"),                  type = "number",   overtitle = gettext("Classified"))
  classificationTable$addColumnInfo(name = "meanPosterior", title = gettext("Mean posterior probability"),  type = "number",   overtitle = gettext("Classified"))

  if (!.saSurvivalReady(options) || jaspBase::isTryError(fit))
    return(classificationTable)

  mixture    <- attr(fit, "mixture")
  components <- mixture[["components"]]
  posterior  <- mixture[["posterior"]]
  dataset    <- attr(fit, "dataset")
  weights    <- if (options[["weights"]] != "") dataset[[options[["weights"]]]] else rep(1, nrow(posterior))

  # observations are classified to the component with the highest posterior probability
  assigned <- max.col(posterior, ties.method = "first")
  data     <- data.frame(
    component     = gettextf("Component %1$i", seq_len(components)),
    probability   = .sapmComponentEstimates(fit, .sapmFamily(mixture[["family"]]), components)[["probabilities"]],
    count         = vapply(seq_len(components), function(k) sum(weights[assigned == k]), numeric(1)),
    meanPosterior = vapply(seq_len(components), function(k) {
      if (!any(assigned == k)) return(NA_real_)
      return(stats::weighted.mean(posterior[assigned == k, k], weights[assigned == k]))
    }, numeric(1))
  )
  data$proportion <- data$count / sum(weights)

  data$subgroup        <- NA
  data$distribution    <- NA
  data$components      <- NA
  data$model           <- NA
  data$subgroup[1]     <- attr(fit, "subgroup")
  data$distribution[1] <- attr(fit, "distribution")
  data$components[1]   <- attr(fit, "components")
  data$model[1]        <- attr(fit, "modelTitle")

  # relative entropy (1 = perfectly separated components)
  entropy <- -sum(weights * rowSums(ifelse(posterior > 0, posterior * log(posterior), 0)))
  entropy <- 1 - entropy / (sum(weights) * log(components))

  # add footnotes
  if (!is.null(attr(fit, "label")) && attr(fit, "label") != "")
    classificationTable$addFootnote(attr(fit, "label"))
  classificationTable$addFootnote(gettextf("Observations are classified to the component with the highest posterior probability. The relative entropy of the classification is %1$.3f (values close to 1 indicate well-separated components).", entropy))

  classificationTable$setData(data)
  classificationTable$showSpecifiedColumnsOnly <- TRUE

  return(classificationTable)
}
.sapmComponentPlotFun           <- function(fit, options) {

  fit <- fit[[1]]

  estimateTitle <- switch(
    options[["mixtureComponentPlotType"]],
    "survival"           = gettext("Survival Probability"),
    "failureProbability" = gettext("Failure Probability"),
    "density"            = gettext("Density"),
    "hazard"             = gettext("Hazard")
  )

  if (!.saSurvivalReady(options) || jaspBase::isTryError(fit))
    return(createJaspPlot(title = estimateTitle))

  plotData <- try(.sapmComponentPlotData(fit, options))

  if (jaspBase::isTryError(plotData)) {
    tempPlot <- createJaspPlot(title = estimateTitle)
    tempPlot$setError(gettext("The model failed to produce predictions. Consider simplifying the model."))
    return(tempPlot)
  }

  hasLevel <- length(unique(plotData[["Level"]])) > 1
  options[["predictionsConfidenceInterval"]] <- options[["predictionsConfidenceInterval"]] && !all(is.na(plotData[["lCi"]]))

  # the mixture is displayed in black and the components follow the color palette
  colors <- c("black", jaspGraphs::JASPcolors(options[["colorPalette"]], asFunction = TRUE)(attr(fit, "components")))
  names(colors) <- levels(plotData[["Component"]])

  plot <- ggplot2::ggplot(data = plotData)

  if (options[["mixtureComponentPlotType"]] %in% c("survival", "failureProbability") && options[["mixtureComponentPlotKaplanMeier"]] && options[["censoringType"]] == "right") {
    kmTable <- .sapKaplanMeierStepData(attr(fit, "dataset"), options, failureProbability = options[["mixtureComponentPlotType"]] == "failureProbability")
    plot    <- plot + jaspGraphs::geom_line(mapping = ggplot2::aes(x = at, y = estimate), data = kmTable, color = "grey60")
  }

  if (options[["predictionsConfidenceInterval"]]) {
    aesCall <- list(
      x        = as.name("at"),
      ymin     = as.name("lCi"),
      ymax     = as.name("uCi"),
      group    = if (hasLevel) as.name("Level")
    )
    geomCall <- list(mapping = do.call(ggplot2::aes, aesCall[!sapply(aesCall, is.null)]), data = plotData[plotData[["Component"]] == levels(plotData[["Component"]])[1], ], fill = "grey60", alpha = 0.30)
    plot <- plot + do.call(ggplot2::geom_ribbon, geomCall)
  }

  aesCall <- list(
    x        = as.name("at"),
    y        = as.name("estimate"),
    color    = as.name("Component"),
    linetype = if (hasLevel) as.name("Level")
  )
  geomCall <- list(mapping = do.call(ggplot2::aes, aesCall[!sapply(aesCall, is.null)]))
  plot <- plot + do.call(jaspGraphs::geom_line, geomCall) +
    ggplot2::scale_color_manual(values = colors, name = gettext("Component"))

  xBreaks <- jaspGraphs::getPrettyAxisBreaks(range(plotData[["at"]], na.rm = TRUE))
  yBreaks <- jaspGraphs::getPrettyAxisBreaks(range(c(
    plotData[["estimate"]],
    if (options[["predictionsConfidenceInterval"]]) plotData[["lCi"]],
    if (options[["predictionsConfidenceInterval"]]) plotData[["uCi"]]), na.rm = TRUE))

  plot <- plot + jaspGraphs::scale_x_continuous(breaks = xBreaks, limits = range(xBreaks), oob = scales::oob_keep) +
    jaspGraphs::scale_y_continuous(breaks = yBreaks, limits = range(yBreaks), oob = scales::oob_keep) +
    ggplot2::ylab(estimateTitle) + ggplot2::xlab(gettext("Time"))

  # the detailed theme is available only for the survival probability plots
  if (options[["plotTheme"]] == "detailed")
    options[["plotTheme"]] <- "jasp"
  plot <- .sapPredictionPlotAddTheme(plot, options)

  tempPlot <- createJaspPlot(width = 550, height = 320)
  tempPlot$plotObject <- plot

  return(tempPlot)
}
.sapmComponentPlotData          <- function(fit, options) {

  mixture    <- attr(fit, "mixture")
  family     <- .sapmFamily(mixture[["family"]])
  components <- mixture[["components"]]
  type       <- options[["mixtureComponentPlotType"]]
  # the components might change rapidly, the time steps are not rounded for a smooth display
  options[["predictionsLifeTimeRoundSteps"]] <- FALSE
  times      <- .sapOptions2PredictionTime(options, fit, type = "mixtureComponents", plot = TRUE)
  ci         <- options[["predictionsConfidenceInterval"]]
  level      <- options[["predictionsConfidenceIntervalLevel"]]

  # the functions are evaluated with the covariate dependent parameters supplied by flexsurv
  mixtureDensity    <- function(t, start, ...) fit[["dfns"]][["d"]](t, ...)
  componentFunction <- function(k) {
    function(t, start, ...) {
      arguments  <- list(...)
      n          <- max(length(t), lengths(arguments))
      parameters <- lapply(stats::setNames(arguments[paste0(family[["pars"]], k)], family[["pars"]]), rep_len, length.out = n)
      t          <- rep_len(t, n)
      switch(
        type,
        "survival"           = do.call(family[["p"]], c(list(t), parameters, list(lower.tail = FALSE))),
        "failureProbability" = do.call(family[["p"]], c(list(t), parameters)),
        "density"            = .sapmStickBreaking(if (components > 1) do.call(cbind, lapply(arguments[paste0("v", seq_len(components - 1))], rep_len, length.out = n)), n)[, k] * do.call(family[["d"]], c(list(t), parameters)),
        "hazard"             = do.call(family[["h"]], c(list(t), parameters))
      )
    }
  }

  mixtureSummary <- switch(
    type,
    "survival"           = summary(fit, type = "survival", t = times, ci = ci, cl = level),
    "failureProbability" = summary(fit, type = "survival", t = times, ci = ci, cl = level),
    "density"            = summary(fit, fn = mixtureDensity, t = times, ci = ci, cl = level),
    "hazard"             = summary(fit, type = "hazard", t = times, ci = ci, cl = level)
  )
  componentSummaries <- lapply(seq_len(components), function(k) summary(fit, fn = componentFunction(k), t = times, ci = FALSE))

  componentLabels <- c(gettext("Mixture"), gettextf("Component %1$i", seq_len(components)))
  out <- list()
  for (j in seq_along(mixtureSummary)) {

    levelLabel <- if (length(mixtureSummary) > 1) decodeColNames(names(mixtureSummary)[j]) else NA

    mixtureData <- data.frame(
      at        = mixtureSummary[[j]][[1]],
      estimate  = mixtureSummary[[j]][["est"]],
      lCi       = if (ci) mixtureSummary[[j]][["lcl"]] else NA,
      uCi       = if (ci) mixtureSummary[[j]][["ucl"]] else NA,
      Component = componentLabels[1],
      Level     = levelLabel
    )
    if (type == "failureProbability") {
      mixtureData[["estimate"]] <- 1 - mixtureData[["estimate"]]
      mixtureData[c("lCi", "uCi")] <- 1 - mixtureData[c("uCi", "lCi")]
    }

    out[[length(out) + 1]] <- mixtureData
    for (k in seq_len(components)) {
      out[[length(out) + 1]] <- data.frame(
        at        = componentSummaries[[k]][[j]][["time"]],
        estimate  = componentSummaries[[k]][[j]][["est"]],
        lCi       = NA,
        uCi       = NA,
        Component = componentLabels[k + 1],
        Level     = levelLabel
      )
    }
  }

  out <- do.call(rbind, out)
  out[["Component"]] <- factor(out[["Component"]], levels = componentLabels)

  # set any Inf to NA
  out[["estimate"]][is.infinite(out[["estimate"]])] <- NA
  out[["lCi"]][is.infinite(out[["lCi"]])]           <- NA
  out[["uCi"]][is.infinite(out[["uCi"]])]           <- NA

  return(out)
}

# mixture messages
.sapmFitMessages                <- function(fit, options) {

  mixture  <- attr(fit, "mixture")
  messages <- NULL

  if (is.null(mixture))
    return(messages)

  if (!mixture[["emConverged"]])
    messages <- c(messages, gettextf("The EM algorithm did not converge within %1$i iterations.", mixture[["emIterations"]]))
  if (mixture[["emGap"]])
    messages <- c(messages, gettext("The EM algorithm had not converged; the reported fit is the polished solution."))
  if (length(mixture[["collapsed"]]) > 0)
    messages <- c(messages, sprintf(ngettext(
      length(mixture[["collapsed"]]),
      "Component %1$s has a negligible weight or diverging parameters; consider fewer components.",
      "Components %1$s have a negligible weight or diverging parameters; consider fewer components."
    ), paste(mixture[["collapsed"]], collapse = ", ")))
  if (mixture[["hessianWarning"]])
    messages <- c(messages, gettext("The Hessian of the likelihood is not positive definite; the standard errors might be unreliable."))

  return(messages)
}
.sapmSummaryMessages            <- function(fit, options) {

  messages <- NULL
  mixtures <- Filter(function(x) !jaspBase::isTryError(x) && !is.null(attr(x, "mixture")), fit)

  if (length(mixtures) == 0)
    return(messages)

  messages <- gettextf(
    "Mixture models were estimated with the EM algorithm (%1$s initialization, %2$i restarts) followed by a direct maximization of the likelihood.",
    switch(
      options[["mixtureInitialization"]],
      "kmeans"    = gettext("k-means"),
      "quantiles" = gettext("quantile-based"),
      "random"    = gettext("random")
    ),
    options[["mixtureRestarts"]]
  )

  for (i in seq_along(mixtures)) {
    fitMessages <- .sapmFitMessages(mixtures[[i]], options)
    if (length(fitMessages) > 0)
      messages <- c(messages, paste0(.sapmCellLabel(mixtures[[i]], options), ": ", fitMessages))
  }

  return(messages)
}
.sapmCellLabel                  <- function(fit, options) {
  return(gettextf(
    "%1$s model %2$s with %3$s%4$s",
    attr(fit, "distribution"),
    attr(fit, "modelTitle"),
    .sapComponentsLabel(attr(fit, "components")),
    if (options[["subgroup"]] != "") paste0(" (", attr(fit, "subgroupLabel"), ")") else ""
  ))
}
.sapmCoefficientsNames          <- function(coeffTable, fit) {

  mixture     <- attr(fit, "mixture")
  family      <- .sapmFamily(mixture[["family"]])
  components  <- mixture[["components"]]
  names       <- rownames(fit[["res"]])

  coeffTable[["mixtureComponent"]] <- NA_integer_

  for (k in seq_len(components)) {

    # component parameters
    for (par in family[["pars"]])
      coeffTable[["coefficient"]][names == paste0(par, k)] <- gettextf("%1$s (component %2$i)", par, k)

    # covariate effects on the location parameter of the component (the effects of the first component are not prefixed)
    index <- fit[["covpars"]][fit[["mx"]][[paste0(family[["location"]], k)]]]
    if (length(index) > 0) {
      coeffTable[["mixtureComponent"]][index] <- k
      if (k > 1) {
        prefix <- paste0(family[["location"]], k, "(")
        coeffTable[["coefficient"]][index] <- substr(names[index], nchar(prefix) + 1, nchar(names[index]) - 1)
      }
    }
  }

  # stick-breaking weights (the standard errors are not reported by flexsurv for the logit transformation)
  for (k in seq_len(components - 1)) {
    index <- names == paste0("v", k)
    coeffTable[["coefficient"]][index] <- if (components == 2) gettextf("Mixing probability (component %1$i)", k) else gettextf("Stick-breaking weight %1$i", k)
    coeffTable[["se"]][index]          <- fit[["res.t"]][index, "se"] * fit[["res"]][index, "est"] * (1 - fit[["res"]][index, "est"])
  }

  return(coeffTable)
}
