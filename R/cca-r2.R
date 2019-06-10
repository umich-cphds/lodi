#' Single pollutant complete case analysis.
#'
#' lod_cca is a helper function that does complete case analysis for
#' single polluatant models. Its primary use is to compare it to clmi.
#' @param formula A R formula in the form outcome ~ exposure + covariates.
#' @param df A data.frame that contains the variables \code{formula}
#'  references.
#' @param type The type of regression to perform. Acceptable options are
#'   linear and logistic.
#' @examples
#' library(lodi)
#' # load lodi's toy data
#' data("toy_data")
#' x <- lod_cca(case_cntrl ~ poll + smoking + gender, toy_data, logistic)
#' # see the fit model
#' x$model
#' @export
lod_cca <- function(formula, df, type)
{
  if (class(formula) != "formula")
      stop("formula must be a formula")
  if (!is.data.frame(df))
        stop("df must be a data.frame")
  type <- as.character(substitute(type))
  type <- match.arg(type, c("linear", "logistic"))
  df <- stats::na.omit(df)

  mod.fr <- stats::model.frame(formula, data = df)
  if (type == "linear")
    model <- stats::lm(mod.fr)
  else
    model <- stats::glm(mod.fr, family = stats::binomial())
  list(model = model, formula = formula, data = df)
}

#' Single pollutant \code{sqrt(2)} imputation.
#'
#' lod_root2 is a helper function that does \code{sqrt(2)} single imputation for
#' single polluatant models. It is pretty common to recode observed values
#' below the limit of detection as \code{lod / sqrt(2)}. The function's primary
#' purpose is to compare it to clmi.
#' @param formula A R formula in the form \code{outcome ~ exposure + covariates.}
#' @param df A data.frame that contains the variables \code{formula}
#'  references.
#' @param lod name of the limit of detection variable
#' @param type The type of regression to perform. Acceptable options are
#'   "linear" and "logistic".
#' @note
#'  Depending on the transformation used, a "Complicated transformation" error
#'  may occur. For example, the transformation \code{a * exposure} will cause an
#'  error.In thise case, define a transformation function as
#'  \code{f <- function(exposure) a * exposure} and use \code{f} in your
#'  formula. This techincal limitation is unavoidable at the moment.
#' @examples
#' # load lodi's toy data
#' library(lodi)
#' data("toy_data")
#' lodi.out <- lod_root2(case_cntrl ~ poll + smoking + gender, toy_data, lod,
#'                         logistic)
#' # see the fit model
#' lodi.out$model
#'
#' # we can log transform poll to make it normally distributed
#' lodi.out <- lod_root2(case_cntrl ~ log(poll) + smoking + gender, toy_data,
#'                         lod, logistic)
#' lodi.out$model
#'
#' # transforming the exposure results in a new column being added to data,
#' # representing the transformed lod.
#' head(lodi.out$data)
#'
#' # You can even define your own transformation functions and use them
#' f <- function(x) exp(sqrt(x))
#' lodi.out <- lod_root2(case_cntrl ~ f(poll) + smoking + gender, toy_data, lod,
#'                         logistic)
#' head(lodi.out$data)
#' @export
lod_root2 <- function(formula, df, lod, type)
{
  if (class(formula) != "formula")
    stop("formula must be a formula")
  if (!is.data.frame(df))
      stop("df must be a data.frame")

  # get the transformation on the exposure
  transform.init <- str2lang(labels(stats::terms(formula))[1])
  # get the exposure
  exposure  <- all.vars(transform.init)
  if (length(exposure) > 1)
    stop("Complicated transformation on exposure. See help for fix.")
  # Calculate the transformation function
  assign(exposure, quote(x))
  transform <- eval(substitute(substitute(transform.init)))
  t.function <- function(x) x
  body(t.function) <- transform
  environment(t.function) <- new.env()

  # calculate the column name of the transformed lod variable
  assign("x", substitute(lod))
  transform.lod <- deparse(eval(substitute(substitute(transform))))

  lod <- deparse(substitute(lod))
  if (is.null(df[[lod]]))
    stop(sprintf("%s not in data", lod))

  type <- deparse(substitute(type))
  type <- match.arg(type, c("linear", "logistic"))
  # impute NAs as lod / sqrt(2)
  tmp <- df[[exposure]]
  for (i in 1:nrow(df))
    tmp[i] <- ifelse(is.na(tmp[i]), df[[lod]][i] / sqrt(2), tmp[i])

  df[[exposure]] <- tmp
  df[[transform.lod]] <- t.function(df[[lod]])
  df[[deparse(transform.init)]] <- t.function(df[[exposure]])

  mod.fr <- stats::model.frame(formula, data = df)
  if (type == "lm")
    model <- stats::lm(mod.fr)
  else
    model <- stats::glm(mod.fr, family = stats::binomial())
  list(model = model, formula = formula, data = df, t.function = t.function)
}
