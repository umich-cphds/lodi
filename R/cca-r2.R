#' Single pollutant complete case analysis.
#'
#' lod_cca is a helper function that does complete case analysis for
#' single polluatant models. Its primary use is to compare it to clmi.
#' @param formula A R formula in the form outcome ~ exposure + covariates.
#' @param data A data.frame that contains the variables \code{formula}
#'  references.
#' @param type The type of regression to perform. Acceptable options are
#'   linear and logistic.
#' @examples
#' library(lodi)
#' # load lodi's toy data
#' data("toy-example")
#' x <- lod_cca(case_cntrl ~ poll + smoking + gender, toy.data, logistic)
#' # see the fit model
#' x$model
#' @export
lod_cca <- function(formula, data, type)
{
  if (class(formula) != "formula")
      stop("formula must be a formula")
  if (!is.data.frame(data))
        stop("data must be a data.frame")
  type <- as.character(substitute(type))
  type <- match.arg(type, c("linear", "logistic"))
  data <- stats::na.omit(data)

  mod.fr <- stats::model.frame(formula, data = data)
  if (type == "linear")
    model <- stats::lm(mod.fr)
  else
    model <- stats::glm(mod.fr, family = stats::binomial())
  list(model = model, formula = formula, data = data)
}
# lod_root2(case_cntrl ~ poll + smoking + gender, toy.data, poll, lod, logistic)

#' Single pollutant \code{sqrt(2)} imputation.
#'
#' lod_root2 is a helper function that does \code{sqrt(2)} single imputation for
#' single polluatant models. It is pretty common to recode observed values
#' below the limit of detection as \code{lod / sqrt(2)}. The function's primary
#' purpose is to compare it to clmi.
#' @param formula A R formula in the form outcome ~ exposure + covariates.
#' @param data A data.frame that contains the variables \code{formula}
#'  references.
#' @param exposure name of the exposure variable
#' @param lod name of the limit of detection variable
#' @param type The type of regression to perform. Acceptable options are
#'   "linear" and "logistic".
#' @examples
#' # load lodi's toy data
#' library(lodi)
#' data("toy_data")
#' lodi.out <- lod_root2(case_cntrl ~ poll + smoking + gender, toy_data, poll,
#'                         lod, logistic)
#' # see the fit model
#' lodi.out$model
#'
#' # we can log transform poll to make it normally distributed
#' lodi.out <- lod_root2(case_cntrl ~ log(poll) + smoking + gender, toy_data,
#'                         poll, lod, logistic)
#' lodi.out$model
#'
#' # transforming the exposure results in a new column being added to data,
#' # representing the transformed lod.
#' head(lodi.out$data)
#'
#' # You can even define your own transformation functions and use them
#' f <- function(x) exp(sqrt(x))
#' lodi.out <- lod_root2(case_cntrl ~ f(poll) + smoking + gender, toy_data,
#'                         poll, lod, logistic)
#' head(lodi.out$data)
#' @export
lod_root2 <- function(formula, data, exposure, lod, type)
{
  if (class(formula) != "formula")
    stop("formula must be a formula")
  if (!is.data.frame(data))
      stop("data must be a data.frame")

  # black magic
  transform <- str2lang(labels(terms(formula))[1])
  eval(substitute(assign(deparse(substitute(exposure)), quote(lod))))

  t.lod <- eval(substitute(substitute(transform)))
  exposure <- deparse(substitute(exposure))
  lod      <- deparse(substitute(lod))
  type     <- deparse(substitute(type))

  if (is.null(data[[lod]]))
    stop(sprintf("%s not in data", lod))
  type <- match.arg(type, c("linear", "logistic"))
  tmp <- data[[exposure]]
  for (i in 1:nrow(data))
    tmp[i] <- ifelse(is.na(tmp[i]), data[[lod]][i] / sqrt(2), tmp[i])
  data[[exposure]] <- tmp
  data[[deparse(t.lod)]] <- eval(t.lod, data)
  data[[deparse(transform)]] <- eval(transform, data)
  mod.fr <- stats::model.frame(formula, data = data)
  if (type == "lm")
    model <- stats::lm(mod.fr)
  else
    model <- stats::glm(mod.fr, family = stats::binomial())
  list(model = model, formula = formula, data = data)
}
