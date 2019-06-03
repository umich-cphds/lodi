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
#' @param df a data.frame with contaminant concentration, batch number,
#'   covariates used in imputation, precision variables.
#' @param contaminant string corresponding to name of the contaminant variable.
#' @param batch string corresponding to name of batch indicator variable.
#' @param outcome string corresponding to name of binary health outcome
#' @param contam.covars covariates to use in imputing contaminant. Can be NULL.
#' @param lod.info a two column data.frame. The first column must be named
#'   "batch.info" corresponds to each batch. The second is named "lod" and
#'   contains the lods for each batch.
#' @param n.imps number of datasets to impute. Default is 5.
#' @param seed for reproducability
#' @param t.function transformation function to apply to contaminant to make
#'   it normally distributed. Default is no transformation
#' @note \code{clmi} only supports categorical variables that are numeric,
#'   (i.e., not factors or characters). You can use the  \code{model.matrix}
#'   function to convert data frame with factors to numeric design matrix and
#'   convert that matrix back into a data frame.
#' @examples
#' library(clmi)
#' data("toy-example")
#' covars <- c("smoking", "gender")
#' clmi.out <- clmi(toy.data, "poll", "batch1", "case_cntrl", covars,
#'                    lod.info = data.frame(batch.info = c("1","0"),
#'                      lod = c(0.8, 0.65)), 20, 12345,
#'                    t.function = function(x) log(x))
#'
#' @references Boss J, Mukherjee B, Ferguson KK, et al. Estimating
#'   outcome-exposure associations when exposure biomarker detection limits
#'   vary across batches. To appear in \emph{Epidemiology}. 2019.
#' @export
clmi <- function(df, contaminant, batch, outcome, contam.covars, lod.info,
                   n.imps = 5, seed, t.function = function(x) x)
{

  if (!is.data.frame(df))
    stop("df must be a data.frame.")

  if (!is.character(contaminant) || length(contaminant) > 1)
    stop("contaminant must be a length 1 character vector.")
  if (is.null(df[[contaminant]]))
    stop(sprintf("%s not in df.\n", contaminant))
  if (!any(is.na(df[[contaminant]])))
    stop("contaminant values below lod must be coded as NA")

  if (!is.character(batch) || length(batch) > 1)
    stop("batch must be a length 1 character vector.")
  if (is.null(df[[batch]]))
    stop(sprintf("%s not in df.\n", batch))

  if (!is.character(outcome) || length(outcome) > 1)
    stop("outcome must be a length 1 character vector.")
  if (is.null(df[[outcome]]))
    stop(sprintf("%s not in df.\n", outcome))
  if (any(is.na(df[[outcome]])))
    stop(sprintf("%s contains missing values.\n", outcome))

  if (!is.character(contam.covars))
    stop("contam.covars must be a character vector.")
  for (var in contam.covars) {
    if (is.null(df[[var]]))
      stop(sprintf("df does not contain %s.\n", var))
    if (any(is.na(df[[var]])))
      stop(sprintf("%s contains missing values.\n", var))
  }

  if (!is.data.frame(lod.info))
    stop("lod.info must be a data.frame.")
  if (is.null(lod.info$batch.info))
    stop("lod.info must contain a column named \"batch.info\".")

  if (is.null(lod.info$lod))
    stop("lod.info must contain a column named \"lod\".")
  if (!is.numeric(lod.info$lod))
    stop("lod.info$lod must be a numeric vector.")

  if (!is.numeric(seed))
    stop("seed must be a number.")

  if (!is.numeric(n.imps))
    stop("n.imps must be an integer.")
  if (n.imps < 1)
    stop("n.imps must be >= 1")

  if (!is.function(t.function))
    stop("t.function must be a function")

  vars <- c(contam.covars, outcome, contaminant)
  if (any(sapply(vars, function(x) !is.numeric(df[[x]]))))
    stop("clmi only supports floating point / integer variables.")

  set.seed(seed)
  # Organize Data into non-detects and detects
  if (!is.null(df$lod))
    stop("TBD")

  df$lod <- sapply(df[[batch]], function(b)
    lod.info[lod.info$batch.info == b, "lod"]
  )

  t.imp.contam <- paste0(contaminant, "_transform", "_imputed")
  df[[t.imp.contam]] <- df[[contaminant]]

  obs.above.lod <- !is.na(df[[t.imp.contam]])
  # Subjects with concentration above LOD
  above.lod <- df[obs.above.lod,]
  # Subjects with concentration below LOD for each batch
  below.lod <- df[!obs.above.lod,]

  above.matrix <- as.matrix(above.lod[, contam.covars])
  below.matrix <- as.matrix(below.lod[, contam.covars])

  # Perform Multiple Imputation
  imp <- rep(list(below.lod), n.imps)
  for (j in 1:n.imps) {
    #Bootstrap data
    df.bs <- df[sample(nrow(df), nrow(df), T),]

    above.lod.bs <- df.bs[!is.na(df.bs[[t.imp.contam]]),]
    above.lod.bs[[t.imp.contam]] <- t.function(above.lod.bs[[t.imp.contam]])

    below.lod.bs <- df.bs[is.na(df.bs[[t.imp.contam]]),]
    below.lod.bs$lod <- t.function(below.lod.bs$lod)

    above.matrix.bs <- as.matrix(above.lod.bs[, contam.covars])
    below.matrix.bs <- as.matrix(below.lod.bs[, contam.covars])

    # calculates the means of (1, contaminant, covariates) given theta
    mu <- function(theta, outcome, covars)
    {
      theta[1] + theta[2] * outcome + covars %*% theta[-(1:3)]
    }

    # objective function for mle
    # is smooth and convex. unique global minimum
    objective <- function(theta)
    {
      mu.b <- mu(theta, below.lod.bs[[outcome]], below.matrix.bs)
      mu.a <- mu(theta, above.lod.bs[[outcome]], above.matrix.bs)
      -sum(log(stats::pnorm(below.lod.bs$lod, mu.b, sqrt(theta[3])))) +
        sum((above.lod.bs[[t.imp.contam]] - mu.a)^2) / (2 * theta[3]) +
        nrow(above.lod.bs) * 0.5 * log(2 * pi * theta[3])
    }
    # get MLE for Bootstrapped sample
    theta <- c(0, 0, 1, rep(0, length(contam.covars)))
    # set lower bounds for theta(3) (variance)
    lower <- c(-Inf, -Inf, 1e-12, rep(-Inf, length(contam.covars)))
    mle   <- stats::optim(theta, objective, method = "L-BFGS-B", lower = lower)

    # Impute missing values
    mus   <- mu(mle$par, below.lod[[outcome]], below.matrix)
    sigma <- sqrt(mle$par[3])

    normalize     <- function(x, mu, sd) (x - mu) / sd
    inv.normalize <- function(x, mu, sd) x * sd + mu

    probs <- stats::pnorm(normalize(t.function(below.lod$lod), mus, sigma))
    zs    <- stats::qnorm(stats::runif(nrow(below.lod)) * probs)
    vals  <- inv.normalize(zs, mus, sigma)
    imp[[j]][[t.imp.contam]] <- as.vector(vals)
  }

  # Check MLE estimation for original dataset
  objective <- function(theta)
  {
    mu.b <- mu(theta, below.lod[[outcome]], below.matrix)
    mu.a <- mu(theta, above.lod[[outcome]], above.matrix)
    -sum(log(stats::pnorm(t.function(below.lod$lod), mu.b, sqrt(theta[3])))) +
      sum((t.function(above.lod[[t.imp.contam]]) - mu.a)^2) / (2 * theta[3]) +
      nrow(above.lod) * 0.5 * log(2 * pi * theta[3])
  }

  hessian <- stats::optimHess(theta, objective)

  fisher.inf <- -solve(-hessian)
  prop.sigma <- sqrt(diag(fisher.inf))
  #Aggregate imputed datasets with observed dataset
  imputed.dfs <- imp
  above.lod[[t.imp.contam]] <- t.function(above.lod[[t.imp.contam]])
  imputed.dfs <- lapply(imputed.dfs, function(df) rbind(df, above.lod))

  structure(
    list(nimp = n.imps, imputed.dfs = imputed.dfs, transform = t.function,
         mle = mle$par, fisher.inf = fisher.inf,
         vars = list(contaminant = contaminant, batch = batch,
           outcome = outcome, contam.covars = contam.covars, lod.info = lod.info,
           t.imp.contam = t.imp.contam)),
    class = "clmi.out")
}

#' Calculate pooled estimates from \code{clmi.out} objects
#'
#' @param clmi.out an object generated by clmi.
#' @param regression.type type of regression to pool. Valid types are
#'   "logistic" and "linear".
#' @param outcome.covars a character vector of variables associated
#'   with the outcome. Only include variable names for variables that were
#'   not used during imputation (variables that were used for imputation
#'   are automatically adjusted for).
#' @examples
#' # continue example from clmi
#' # fit model on imputed data and pool results
#' library(clmi)
#' data("toy-example")
#' covars <- c("smoking", "gender")
#' clmi.out <- clmi(toy.data, "poll", "batch1", "case_cntrl", covars,
#'                    lod.info = data.frame(batch.info = c("1","0"),
#'                      lod = c(0.8, 0.65)), 20, 12345,
#'                    t.function = function(x) log(x))
#' results <- pool.clmi(clmi.out, "logistic", NULL)
#'
#' results$output
#' @export
pool.clmi <- function(clmi.out, regression.type = c("logistic, linear"),
                        outcome.covars = NULL)
{
  if (class(clmi.out) != "clmi.out")
    stop("clmi.out is not an clmi.out object.")

  regression.type <- match.arg(regression.type, c("logistic, linear"))

  # Get names used in clmi() function
  imputed.dfs   <- clmi.out$imputed.dfs
  t.imp.contam  <- clmi.out$vars$t.imp.contam
  outcome       <- clmi.out$vars$outcome
  contam.covars <- clmi.out$vars$contam.covars
  n.imps <- clmi.out$nimp
  nobs <- nrow(imputed.dfs[[1]])

  for (covar in outcome.covars) {
    if (is.null(imputed.dfs[[1]][[covar]]))
      stop(sprintf("clmi.out does not contain %s.\n", covar))
    if (!is.numeric(imputed.dfs[[1]][[covar]]))
      stop(sprintf("%s must be a numeric variable.\n", covar))
  }

  # Get number of coefficients to store
  covars   <- unique(c(contam.covars, outcome.covars, t.imp.contam))
  data.fit <- imputed.dfs[[1]][, covars]
  if (is.vector(data.fit))
    num_continuous <- 1
  else
    num_continuous <- sum(sapply(data.fit, is.numeric))

  # Total number of parameters in logistic regression (+ 1 for the intercept)
  n.param <- 1 + num_continuous
  # Get estimates and standard errors for each imputed dataset
  betas       <- matrix(0, n.imps, n.param)
  varcov      <- vector("list", n.imps)
  regressions <- vector("list", n.imps)

  model <- paste0(outcome, " ~ .")
  if (regression.type == "logistic") {
    for(i in 1:n.imps) {
      data.fit <- imputed.dfs[[i]][, c(covars, outcome)]
      logreg <- stats::glm(model, family = "binomial", data = data.fit)

      betas[i,]        <- stats::coefficients(logreg)
      varcov[[i]]      <- summary(logreg)$cov.unscaled
      regressions[[i]] <- logreg
    }
  } else {
    for (i in 1:n.imps) {
      data.fit <- imputed.dfs[[i]][, c(covars, outcome)]
      linreg <- stats::lm(model, data = data.fit)
      betas[i,]        <- stats::coefficients(linreg)
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
