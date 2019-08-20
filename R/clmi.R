#' Censored Likelihood Multiple Imputation
#'
#' This function performs censored likelihood multiple imputation for
#' single-pollutant models where the pollutant of interest is subject to
#' varying detection limits across batches (this function will also work if
#' there is only one distinct detection limit). The function
#' outputs a list containing the imputed datasets and details regarding the
#' imputation procedure (i.e., number of imputed dataset, covariates used to
#' impute the non-detects, etc).
#'
#' \code{clmi} is somewhat picky regarding the \code{formula} parameter. It
#' tries to infer what transformation you'd like to apply to the exposure you
#' are imputing, what the exposure is, and what the outcome is. It attempts to
#' check to make sure that everything is working correctly, but it can fail.
#' Roughly, the rules are:
#' \itemize{
#'   \item The left hand side of formula should be the exposure you are trying
#'     to impute.
#'   \item the exposure may be optionally wrapped in a univariate transformation
#'     function. If the transformation function is not univariate, you ought to
#'     get an error about a "complicated" transformation.
#'   \item The first variable on the right hand side of \code{formula} should be
#'     your binary outcome of interest.
#'}
#' @param formula A formula in the form of \code{exposure ~ outcome + covariates}.
#' That is, the first variable on the right hand side of \code{formula} should
#' be the (binary) outcome of interest.
#' @param df A data.frame with \code{exposure}, \code{outcome} and
#'   \code{covariates}.
#' @param lod Name of limit of detection variable in \code{df}.
#' @param n.imps Number of datasets to impute. Default is 5.
#' @param seed For reproducability.
#' @param verbose If \code{TRUE} (default) print out useful debugging
#'   information while parsing \code{formula}.
#' @note \code{clmi} only supports categorical variables that are numeric,
#'   (i.e., not factors or characters). You can use the \code{model.matrix}
#'   function to convert a data frame with factors to a numeric design matrix
#'   and subsequently convert that matrix back into a data frame using
#'   \code{as.data.frame}.
#' @examples
#' library(lodi)
#'
#' # Note that the outcome of interest is the first variable on the right hand
#' # side of formula.
#' clmi.out <- clmi(poll ~ case_cntrl + smoking + gender, toy_data, lod, 1)
#'
#' # you can specify a transformation to the exposure in the formula
#' clmi.out <- clmi(log(poll) ~ case_cntrl + smoking + gender, toy_data, lod, 1)
#'
#' @references
#'   Boss J, Mukherjee B, Ferguson KK, et al. Estimating outcome-exposure
#'   associations when exposure biomarker detection limits vary across batches.
#'   Epidemiology. 2019;30(5):746-755.
#'   \href{https://doi.org/10.1097/EDE.0000000000001052}{10.1097/EDE.0000000000001052}
#' @export
clmi <- function(formula, df, lod, seed, n.imps = 5, verbose = TRUE)
{
  if (!rlang::is_formula(formula))
    stop("formula must be a formula")
  if (!is.data.frame(df))
    stop("df must be a data.frame.")

  if (verbose)
    print(paste("Formula:", rlang::expr_text(formula)))

  transform.init <- rlang::f_lhs(formula)
  exposure <- all.vars(transform.init)
  if (length(exposure) > 1)
    stop("Complicated transformation on exposure. See help for fix.")

  if (verbose)
    print(sprintf("Exposure variable: %s", exposure))

  vars <- all.vars(rlang::f_rhs(formula))
  outcome <- vars[1]

  if (verbose)
    print(sprintf("Outcome variable: %s", outcome))

  if (length(unique(df[[outcome]])) != 2)
    stop(paste(sprintf("Outcome (%s) variable is non binary.", outcome),
               "The right hand side of formula may be incorrectly ordered.")
    )

  # Calculate the transformation function
  assign(substitute(exposure), quote(x))
  transform <- eval(rlang::expr(substitute(!!transform.init)))
  t.function <- function(x) x
  rlang::fn_body(t.function) <- transform
  rlang::fn_env(t.function) <- new.env()

  if (verbose)
      print(sprintf("Transformation function: %s",
              gsub("\n", "", rlang::expr_text(t.function))
      ))

  lod <- deparse(substitute(lod))
  if (verbose)
    print(sprintf("LOD variable: %s", lod))

  if (is.null(df[[lod]]))
    stop(sprintf("%s not in df", lod))

  if (!is.numeric(df[[lod]]))
   stop(sprintf("%s must be numeric."))

  if (!is.numeric(seed))
    stop("seed must be a number.")
  set.seed(seed)
  if (!is.numeric(n.imps))
    stop("n.imps must be an integer.")
  if (n.imps < 1)
    stop("n.imps must be >= 1")

  vars <- c(exposure, vars)
  if (any(sapply(vars, function(x) !is.numeric(df[[x]]))))
    stop("clmi only supports floating point / integer variables.")
  if (any(sapply(vars[-1], function(x) any(is.na(df[[x]])))))
    stop("The covariates on the rhs of formula cannot contain missing values.")

  # columns to add back at the end
  leftovers <- df[, setdiff(names(df), c(vars,  lod))]

  # will be used to ensure the imputed column ordering is the same
  df.names <- names(df)
  df.rownames <- rownames(df)

  df <- df[, c(vars, lod)]
  t.imp.exp <- paste0(exposure, "_transform", "_imputed")

  df[[t.imp.exp]] <- df[[exposure]]

  obs.above.lod <- !is.na(df[[t.imp.exp]])
  # Subjects with concentration above LOD
  above.lod <- df[obs.above.lod,]
  # Subjects with concentration below LOD for each batch
  below.lod <- df[!obs.above.lod,]

  above.matrix <- as.matrix(above.lod[, vars[-(1:2)]])
  below.matrix <- as.matrix(below.lod[, vars[-(1:2)]])
  # Perform Multiple Imputation
  imp <- rep(list(below.lod), n.imps)
  for (j in 1:n.imps) {
    #Bootstrap data
    df.bs <- df[sample(nrow(df), nrow(df), T),]

    above.lod.bs <- df.bs[!is.na(df.bs[[t.imp.exp]]),]
    above.lod.bs[[t.imp.exp]] <- t.function(above.lod.bs[[t.imp.exp]])

    below.lod.bs <- df.bs[is.na(df.bs[[t.imp.exp]]),]
    below.lod.bs[[lod]] <- t.function(below.lod.bs[[lod]])

    above.matrix.bs <- as.matrix(above.lod.bs[, vars[-(1:2)]])
    below.matrix.bs <- as.matrix(below.lod.bs[, vars[-(1:2)]])

    # calculates the means of (1, exposure, covariates) given theta
    mu <- function(theta, outcome, covars)
    {
      theta[2] + theta[3] * outcome + covars %*% theta[-(1:3)]
    }

    # objective function for mle
    # is smooth and convex. unique global minimum
    objective <- function(theta)
    {
      mu.b <- mu(theta, below.lod.bs[[outcome]], below.matrix.bs)
      mu.a <- mu(theta, above.lod.bs[[outcome]], above.matrix.bs)
      -sum(log(stats::pnorm(below.lod.bs[[lod]], mu.b, sqrt(theta[1])))) +
        sum((above.lod.bs[[t.imp.exp]] - mu.a)^2) / (2 * theta[1]) +
        nrow(above.lod.bs) * 0.5 * log(2 * pi * theta[1])
    }

    # get MLE for Bootstrapped sample
    theta <- c(1, 0, 0, rep(0, length(vars[-(1:2)])))

    # set lower bounds for theta(1) (variance)
    lower <- c( 1e-12, -Inf, -Inf, rep(-Inf, length(vars[-(1:2)])))
    mle   <- stats::optim(theta, objective, method = "L-BFGS-B", lower = lower)

    # Impute missing values
    mus   <- mu(mle$par, below.lod[[outcome]], below.matrix)
    sigma <- sqrt(mle$par[1])

    normalize     <- function(x, mu, sd) (x - mu) / sd
    inv.normalize <- function(x, mu, sd) x * sd + mu

    probs <- stats::pnorm(normalize(t.function(below.lod[[lod]]), mus, sigma))
    zs    <- stats::qnorm(stats::runif(nrow(below.lod)) * probs)
    vals  <- inv.normalize(zs, mus, sigma)
    imp[[j]][[t.imp.exp]] <- as.vector(vals)
  }

  # Check MLE estimation for original dataset
  mu <- function(theta, outcome, covars)
  {
    theta[2] + theta[3] * outcome + covars %*% theta[-(1:3)]
  }

  objective <- function(theta)
  {
    mu.b <- mu(theta, below.lod[[outcome]], below.matrix)
    mu.a <- mu(theta, above.lod[[outcome]], above.matrix)
    -sum(log(stats::pnorm(t.function(below.lod[[lod]]), mu.b, sqrt(theta[1])))) +
      sum((t.function(above.lod[[t.imp.exp]]) - mu.a)^2) / (2 * theta[1]) +
      nrow(above.lod) * 0.5 * log(2 * pi * theta[1])
  }

  mle <- stats::optim(theta, objective, method = "L-BFGS-B", lower = lower,
                        hessian = T)

  param <- mle$par
  names(param) <- c("variance", "intercept", outcome, vars[-(1:2)])
  fisher.inf <- -solve(-mle$hessian)
  prop.sigma <- sqrt(diag(fisher.inf))
  #Aggregate imputed datasets with observed dataset
  imputed.dfs <- imp
  above.lod[[t.imp.exp]] <- t.function(above.lod[[t.imp.exp]])

  # add imputed values, then add back the rest of the columns to to imputed
  # dataframe and make sure it is in the same order as the original
  imputed.dfs <- lapply(imputed.dfs, function(df) {
      df <- rbind(df, above.lod)

      # reorder the rows so they match the original
      df <- df[df.rownames, ]
      df <- cbind(df, leftovers)

      # reorder the columns so they match the original
      df[, c(df.names, t.imp.exp)]
  })

  structure(
    list(formula = formula, nimp = n.imps, imputed.dfs = imputed.dfs,
           t.function = t.function, par.mle = param, fisher.inf = fisher.inf),
    class = "clmi.out")
}

#' Calculate pooled estimates from \code{clmi.out} objects using Rubin's rules
#' @param formula Formula to fit. Exposure variable should end in
#'   \code{_transform_imputed}.
#' @param clmi.out An object generated by clmi.
#' @param type Type of regression to pool. Valid types are
#'   logistic and linear.
#' @examples
#' # continue example from clmi
#' # fit model on imputed data and pool results
#' library(lodi)
#' data("toy_data")
#' clmi.out <- clmi(log(poll) ~ case_cntrl + smoking + gender, toy_data, lod, 1)
#' results <- pool.clmi(case_cntrl ~ poll_transform_imputed + smoking, clmi.out,
#'                        logistic)
#'
#' results$output
#' @export
pool.clmi <- function(formula, clmi.out, type)
{
  if (class(clmi.out) != "clmi.out")
    stop("clmi.out is not an clmi.out object.")

  type <- deparse(substitute(type))
  type <- match.arg(type, c("linear","logistic"))

  # Get names used in clmi() function
  imputed.dfs   <- clmi.out$imputed.dfs
  n.imps <- clmi.out$nimp
  nobs <- nrow(imputed.dfs[[1]])

  # Get number of coefficients to store
  # Get estimates and standard errors for each imputed dataset
  betas       <- c()
  varcov      <- vector("list", n.imps)
  regressions <- vector("list", n.imps)

  if (type == "logistic") {
    for(i in 1:n.imps) {
      data.fit <- imputed.dfs[[i]]
      logreg <- stats::glm(formula, family = "binomial", data = data.fit)
      betas <- rbind(betas, stats::coefficients(logreg))
      varcov[[i]]      <- summary(logreg)$cov.unscaled
      regressions[[i]] <- logreg
    }
  } else {
    for (i in 1:n.imps) {
      data.fit <- imputed.dfs[[i]]
      linreg <- stats::lm(formula, data = data.fit)
      betas <- rbind(betas, stats::coefficients(linreg))
      varcov[[i]]      <- stats::vcov(linreg)
      regressions[[i]] <- linreg
    }
  }

  # Pooled Inference for Multiply Imputed Datasets
  qbar <- apply(betas, 2, mean)
  ubar <- Reduce('+', varcov) / n.imps
  # Each column of this matrix is Qi - Qbar
  q_diff <- apply(betas, 1, function(x) x - qbar)
  ui <- matrix(0, ncol(betas), ncol(betas))
  for (m in 1:ncol(q_diff)) {
    ui <- ui + as.matrix(q_diff[,m]) %*% t(as.matrix(q_diff[,m]))
  }
  b <- ui / (n.imps-1)
  t <- ubar + (1 + (1 / n.imps)) * b
  gamma <- ((1 + (1 / n.imps)) * diag(b)) / diag(t)
  r <- ((1 + (1 / n.imps)) * diag(b)) / diag(ubar)
  v <- (n.imps - 1) * (1 + (1 / r))^2
  v0 <- nobs - ncol(betas)
  denom_adj_df <- ((1 - gamma) * v0 * (v0 + 1)) / (v0 + 3)
  vstar <- ((1 / v) + (1 / denom_adj_df))^(-1)

  #Get Confidence Intervals
  LCL <- qbar - stats::qt(0.975, vstar) * sqrt(diag(t))
  UCL <- qbar + stats::qt(0.975, vstar) * sqrt(diag(t))

  #Get p-values
  p_vals <- 2 * stats::pt(abs(qbar / sqrt(diag(t))), df = vstar, lower.tail = F)

  #Summary of each regression models on each imputed datasets
  regression.summaries <- lapply(regressions, summary)

  list(
    output = data.frame(est = qbar, se = sqrt(diag(t)), df = vstar,
      p.values = p_vals, LCL.95 = LCL, UCL.95 = UCL),
    pooled.vcov = t, regressions = regressions,
    regression.summaries = regression.summaries)
}
