# Formula fields accept scalar expressions or lists of expressions. Keep parsing
# separate from the contextual range checks and messages at each call site.
.sapNumericExpression <- function(expression, multiple = FALSE, allowSpaces = FALSE) {

  if (!is.character(expression) || length(expression) != 1L)
    return(expression)

  text <- trimws(expression)
  if (multiple)
    text <- paste0("c(", gsub("[;|\r\n\t]+", ",", text), ")")
  value <- try(eval(parse(text = text), envir = baseenv()), silent = TRUE)

  # Coefficient vectors also accept whitespace separators; try them only after
  # normal parsing, so spaces inside an expression keep their usual meaning.
  if (multiple && allowSpaces && inherits(value, "try-error")) {
    text  <- gsub("[;|\r\n\t]+", ",", trimws(expression))
    text  <- paste(strsplit(text, "[[:space:]]+")[[1L]], collapse = ",")
    value <- try(eval(parse(text = paste0("c(", text, ")")), envir = baseenv()), silent = TRUE)
  }

  return(value)
}

.sapFiniteNumeric <- function(value, expectedLength = NULL) {

  valid <- is.numeric(value) && !is.complex(value) && length(value) > 0L && all(is.finite(value))
  if (!is.null(expectedLength))
    valid <- valid && length(value) == expectedLength

  return(valid)
}
