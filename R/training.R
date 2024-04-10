# Utilities used for model training, predictions and quality control.

# Model selection and training --------

#' Regularized linear model tuning and training.
#'
#' @description
#' Selects the optimal value of the shrinkage parameter \eqn{\lambda} by
#' minimizing model deviance/error or maximizing model accuracy
#' in cross-validation.
#' Subsequently, \code{\link[glmnet]{glmnet}} models are trained
#' (one for each column of the `y` response matrix) and returned along with the
#' optimal \eqn{\lambda} values, basic cross-validation stats and,
#' for Gaussian and Poisson family models, with additional metrics of goodness
#' of fit.
#'
#' @details
#' Model selection via cross-validation is performed by
#' \code{\link[glmnet]{cv.glmnet}}.
#' Independently of the modeling family, the function returns the tuning
#' statistic in out-of-fold predictions for the best \eqn{\lambda} values.
#' Specifically, for Gaussian family, the function calculates also
#' raw \eqn{R^2} with the formula \eqn{R^2 = 1 - \frac{MSE}{Var(y)}} for
#' out-of-fold predictions, where MSE stands for the out-of-fold
#' mean squared error.
#' For the Gaussian and Poisson family, Pearson's r and Spearman's
#' \eqn{\rho} for the correlation of predictions with the observed response
#' in the entire training data set are computed.
#' Additional evaluations for classification, i.e. binomial and multinomial
#' models are overall accuracy and Cohen's \eqn{\kappa} computed for comparison
#' of the observed and predicted class assiggnment in the training data set.
#'
#' The function will run in parallel provided a backend compatible with
#' \code{\link[furrr]{future_map}}
#' (e.g. declared by \code{\link[future]{plan}}).
#' Observations with missing response entries are silently removed.
#' Note: `glmnet` may have problems with fitting model to messy data
#' (i.e. a lot of NA) and may also fail if there are negative values in the `y`
#' matrix. If this is the case, you may try transform the responses accordingly.
#' The `train()` function is a S3 generic.
#'
#' @references
#' Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized
#' linear models via coordinate descent.
#' J Stat Softw (2010) 33:1–22. doi:10.18637/jss.v033.i01
#'
#' @return An object of class `modTrain` with the following components:
#'
#' * `models`: a list of \code{\link[glmnet]{cv.glmnet}} models used for
#' tuning of \eqn{\lambda}, named after the columns of the `y` matrix.
#'
#' * `stats`: the optimal \eqn{\lambda} parameter values, cross-validated
#' selection measures, and metrics of model performance in the
#' training data set.
#'
#' * `errors`: a list of error messages for responses for which tuning/training
#' failed.
#'
#' Note the `nobs()`, `components()`, `summary()`, `plot()`, and `predict`
#' methods defined for the class.
#'
#' @param x a numeric matrix with training explanatory factors, features are
#' specified in rows, observations are provided as columns. Alternatively, a
#' `modData` class object created e.g. with \code{\link{pre_process}}.
#'
#' @param y a numeric or factor matrix with the observed responses
#' (e.g. phenotype in the training data). The responses are provided as columns,
#' the observations are given in rows. The row names of `y` must match column
#' name of `x`.
#'
#' @param trans_fun a numeric function to transform responses for modeling.
#' A normality-improving transformation such as `log` may be provided here.
#' Used only for Gaussian family models and ignored elsewhere.
#'
#' @param standardize logical, should `x` values be standardized prior
#' to modeling?
#'
#' @param alpha numeric mixing parameter, which defines the modeling algorithm
#' (RIDGE regression: `alpha = 0`, LASSO: `alpha = 1`, Elastic Net otherwise).
#' See \code{\link[glmnet]{glmnet}} for details.
#'
#' @param family family of the `glmnet` models. Currently, only 'gaussian',
#' 'binomial', 'poisson', and 'multinomial' are implemented.
#' See: \code{\link[glmnet]{glmnet}} for details.
#'
#' @param type.measure model selection statistic (i.e. loss function). Defaults
#' to deviance. See \code{\link[glmnet]{cv.glmnet}} for details.
#'
#' @param safe logical, should modeling proceed in spite of errors for
#' some responses? .
#' If `TRUE`, any errors are logged in `errors` element of the output. This
#' option may be, however, considerably slower.
#'
#' @param ... additional arguments passed to \code{\link[glmnet]{cv.glmnet}} and
#' \code{\link[glmnet]{glmnet}}.
#'
#' @export train.default
#' @export

  train.default <- function(x,
                            y,
                            trans_fun = identity,
                            standardize = FALSE,
                            alpha = 1,
                            family = 'gaussian',
                            type.measure = c("default",
                                             "mse",
                                             "deviance",
                                             "class",
                                             "auc",
                                             "mae",
                                             "C"),
                            safe = FALSE, ...) {

    ## input control ----------

    ## validity of alpha, family and type.measure are checked by the
    ## downstream functions.

    start_time <- Sys.time()

    family <- match.arg(family[1],
                        c('gaussian', 'binomial', 'poisson', 'multinomial'))

    stopifnot(is.logical(safe))
    stopifnot(is.function(trans_fun))

    stopifnot(is.logical(standardize))

    type.measure <-
      match.arg(type.measure[1],
                c("default",
                  "mse",
                  "deviance",
                  "class",
                  "auc",
                  "mae",
                  "C"))

    inputs <- list(expl = x, pheno = y)

    if(any(!map_lgl(inputs, is.matrix))) {

      stop("'x' and 'y' have to be numeric matrices.", call. = FALSE)

    }

    if(any(!map_lgl(inputs, is.numeric))) {

      stop("'x' and 'y' have to be numeric matrices.", call. = FALSE)

    }

    if(sum(is.na(x)) > 0) {

      stop("There are NA in 'x': please remove them prior to modeling.",
           call. = FALSE)

    }

    common_observations <- base::intersect(colnames(x), rownames(y))

    if(length(common_observations) == 0) {

      stop("No matched observations in 'x' and 'y'. Wrong modeling matrices?",
           call. = FALSE)

    }

    ## benchmarking -------

    message(paste('\nTraining models for n =', ncol(y), 'responses'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    ## x matrix and response vectors -----------

    x <- x[, common_observations]
    y <- y[common_observations, ]

    x <- t(x)

    if(family == 'gaussian') y <- trans_fun(y)

    if(sum(y < 0, na.rm = TRUE) > 0) {

      warning(paste("There are negative values of in the response matrix 'y'.",
                    "Fitting of glmnet models may fail.",
                    "You may consider transforming the response to avoid",
                    "negative values and NA."),
              call. = FALSE)

    }

    y_vecs <- map(colnames(y), ~y[, .x])

    y_vecs <- set_names(y_vecs, colnames(y))

    y_vecs <- map(y_vecs, ~.x[!is.na(.x)])

    ## model tuning: specific for the modeling family  ------

    message('Model tuning and training')

    if(family %in% c('gaussian', 'poisson')) {

      fits <- fit_regression(y_vecs,
                             x,
                             safe = safe,
                             family = family,
                             type.measure = type.measure,
                             standardize = standardize,
                             alpha = alpha, ...)

    } else {

      fits <- fit_classification(y_vecs,
                                 x,
                                 safe = safe,
                                 family = family,
                                 type.measure = type.measure,
                                 standardize = standardize,
                                 alpha = alpha, ...)

    }

    cv_models <- fits$cv_models

    errors <- cv_models$errors

    if(length(cv_models) == 0) {

      warning('No models could be fitted. Returning errors.', call. = FALSE)

      return(errors)

    }

    y_vecs <- y_vecs[names(cv_models)]

    ## cross-validation stats and the best tunes ----------

    cvm <- NULL

    tunes <- map(cv_models, ~as_tibble(.x[c('lambda', 'cvm')]))

    if(type.measure %in% c('C', 'auc')) {

      best_tunes <-
        map(tunes,
            filter,
            cvm == max(cvm, na.rm = TRUE))

      loss_criterion <- 'max'

    } else {

      best_tunes <-
        map(tunes,
            filter,
            cvm == min(cvm, na.rm = TRUE))

      loss_criterion <- 'min'

    }

    best_tunes <- map(best_tunes, ~.x[1, ])

    response <- NULL
    lambda <- NULL
    loss_name <- NULL
    loss_crietrion <- NULL
    loss_value <- NULL
    n_features <- NULL

    cv_stats <-
      tibble(response = names(cv_models),
             n_observations = map_dbl(y_vecs, length),
             n_features = ncol(x),
             lambda = map_dbl(best_tunes, ~.x[['lambda']]),
             loss_name = type.measure,
             loss_criterion = loss_criterion,
             loss_value = map_dbl(best_tunes, ~.x[['cvm']]))

    if(family == 'gaussian' & type.measure %in% c('default', 'mse')) {

      variance <- map_dbl(y_vecs, var)

      rsq_oof <- NULL

      cv_stats <-
        mutate(cv_stats,
               variance = variance,
               rsq_oof = 1 - loss_value/variance)

    }

    ## evaluation --------

    message('Evalaution of the predictions')

    ## Gaussian and Poisson family:
    ## correlations with of the predicted outcome with the observed outcome
    ##
    ## Classification models: overall accuracy and kappa

    if(family %in% c('gaussian', 'poisson')) {

      ## Gaussian: correlation of the predicted and observed response
      ## Poisson: correlation of the observed response
      ## and linear predictor score

      train_evals <- evaluate_regression(model_lst = cv_models,
                                         lambda = cv_stats$lambda,
                                         y_lst = y_vecs,
                                         x = x)

      pearson <- NULL
      spearman <- NULL

      cv_stats <-
        mutate(cv_stats,
               pearson = train_evals$pearson,
               spearman = train_evals$spearman)

    } else {

      train_evals <- evaluate_classification(model_lst = cv_models,
                                             lambda = cv_stats$lambda,
                                             y_lst = y_vecs,
                                             x = x)

      accuracy <- NULL
      kappa <- NULL

      cv_stats <-
        mutate(cv_stats,
               accuracy = train_evals$accuracy,
               kappa = train_evals$kappa)

    }

    ## number of non-zero coefficients

    n_non_zero <- map2(cv_models,
                       cv_stats$lambda,
                       ~coef(.x$glmnet.fit,
                             s = .y))

    if(family != 'multinomial') {

      n_non_zero <- map_dbl(n_non_zero,
                            ~sum(as.matrix(.x) != 0))

    } else {

      n_non_zero <- map(n_non_zero,
                        map_dbl,
                        ~sum(as.matrix(.x) != 0))

      n_non_zero <- map_dbl(n_non_zero, sum)

    }

    n_coefficients <- NULL

    cv_stats <- mutate(cv_stats,
                       n_coefficients = n_non_zero,
                       .after = n_features)

    ## output -------

    mod_globals <- list(alpha = alpha,
                        family = family,
                        standardize = standardize,
                        type.measure = type.measure)

    dots <- list2(...)

    if(length(dots) > 1) mod_globals <- c(mod_globals, dots)

    return(modTrain(list(models = cv_models,
                         stats = cv_stats,
                         globals = mod_globals,
                         errors = errors)))

  }

#' @rdname train.default
#' @export train.modData
#' @export

  train.modData <- function(x,
                            y,
                            trans_fun = identity,
                            standardize = FALSE,
                            alpha = 1,
                            family = 'gaussian',
                            type.measure = c("default",
                                             "mse",
                                             "deviance",
                                             "class",
                                             "auc",
                                             "mae",
                                             "C"),
                            safe = FALSE, ...) {

    stopifnot(is_modData(x))

    train(x = components(x, 'train'),
          y = y,
          trans_fun = trans_fun,
          standardize = standardize,
          alpha = alpha,
          family = family,
          type.measure = type.measure,
          safe = safe, ...)

  }

# END ----
