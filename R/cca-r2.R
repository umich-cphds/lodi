#' Single pollutant complete case analysis.
#'
#' lod.cca is a helper function that does complete case analysis for
#' single polluatant models. Its primary use is to compare it to clmi.
#' @param formula A R formula in the form outcome ~ exposure + covariates.
#' @param data A data.frame that contains the variables \code{formula}
#'  references.
#' @param type The type of regression to perform. Acceptable options are
#'   "linear" and "logistic".
#' @export
lod.cca <- function(formula, data, type)
{
  if (class(formula) != "formula")
      stop("formula must be a formula")
  if (!is.data.frame(data))
        stop("data must be a data.frame")
  type <- match.arg(type, c("linear", "logistic"))
  data <- stats::na.omit(data)

  mod.fr <- stats::model.frame(formula, data = data)
  if (type == "lm")
    model <- stats::lm(mod.fr)
  else
    model <- stats::glm(mod.fr, family = stats::binomial())
  list(model = model, formula = formula, data = data)
}

#' Single pollutant \code{sqrt(2)} imputation.
#'
#' lod.ro2 is a helper function that does \code{sqrt(2)} single imputation for
#' single polluatant models. It is pretty common to recode rehistered values
#' below the lod as \code{lod / sqrt(2)}. The function's primary purpose is to
#' compare it to clmi.
#' @param formula A R formula in the form outcome ~ exposure + covariates.
#' @param data A data.frame that contains the variables \code{formula}
#'  references.
#' @param exposure Name of the exposure variable
#' @param lod Name of the lod variable
#' @param type The type of regression to perform. Acceptable options are
#'   "linear" and "logistic".
#' @export
lod.root2 <- function(formula, data, exposure, lod, type)
{
  if (class(formula) != "formula")
    stop("formula must be a formula")
  if (!is.data.frame(data))
      stop("data must be a data.frame")
  if (!is.character(lod))
    stop("lod should refer to the name of the lod variable")
  if (is.null(data[[lod]]))
    stop(sprintf("%s not in data", lod))

  type <- match.arg(type, c("lm", "glm"))
  tmp <- data[[exposure]][i]
  for (i in 1:nrow(data))
    tmp[i] <- ifelse(is.na(tmp[i]), data[[lod]][i] / sqrt(2), tmp[i])
  data[[exposure]] <- tmp

  mod.fr <- stats::model.frame(formula, data = data)
  if (type == "lm")
    model <- stats::lm(mod.fr)
  else
    model <- stats::glm(mod.fr, family = stats::binomial())
  list(model = model, formula = formula, data = data)
}
