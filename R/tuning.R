# Wrapper functions for tuning of GLMNET models by repeated CV

#' Tune `GLMNET` models by repeated cross-validation.
#'
#' @description
#' Finds the globally optimal value of the shrinkage parameter `lambda` for
#' regularized linear `GLMNET` models by minimizing or maximizing a cost
#' function in repeated cross-validation.
#'
#' @details
#' In principle, there is only one tuning parameter in regularized linear models
#' build by modeling tools of the `GLMNET` package: the shrinkage parameter
#' `lambda`, whose optimal value can be estimates in a repeated cross-validation
#' setting. Technically, this task is reduced to running
#' \code{\link[glmnet]{cv.glmnet}} for the requested number of repeats.
#' The function will run in parallel if a parallel backend is specified via
#' \code{\link[future]{plan}}.
#'
#' @return a \code{\link{cvTune}} cvlass object with the following components:
#'
#' * `stats`: model performance statistic, i.e. the cost function, for the
#' locally optimal `lambda` values in repeats of the cross-validation. Those
#' are the `lambda` values, one per repeat, which correspond to the optimum of
#' the cost function in out-of-fold predictions. In case there are two local
#' optima, the first one corresponding to a smaller `lambda` is returned.
#'
#' * `lambda`: the globally optimal value of `lambda`m, i.e. the optimum of
#' `lambdas` for particular repeats listed in `stats`. If multiple values of
#' cost function are available for the same value of lambda, then the cost
#' function is averaged per lambda value prior to finding the optimal lambda
#' value.
#'
#' * `model`: a `GLMNET` model fit by \code{\link[glmnet]{glmnet}} for the
#' globally optimal `lambda` value
#'
#' Note the methods \code{\link{nobs.cvTune}}, \code{\link{summary.cvTune}},
#' \code{\link{plot.cvTune}}, and \code{\link{predict.cvTune}} defined
#' for the class.
#'
#' @param x numeric matrix or a \code{\link[stats]{model.matrix}} object
#' with explanatory variables; observations in rows and variables in columns.
#' @param y response object (numeric vector, factor, a matrix or s survival
#' object); please consult \code{\link[glmnet]{glmnet}} for details.
#' @param fold_ids a list of numeric vectors with indexes (one per repeat),
#' which specify the assignment of observations to cross-validation folds.
#' See: \code{\link[glmnet]{cv.glmnet}} for details and
#' \code{\link[caret]{createFolds}} for a handy function returning the
#' fold-assortment vectors for each of the repeats.
#' @param type.measure cost function as specified for
#' \code{\link[glmnet]{cv.glmnet}}. For `type.measure = 'default'`, the model
#' deviance is minimized.
#' @param ... other arguments passed to \code{\link[glmnet]{glmnet}} and
#' \code{\link[glmnet]{cv.glmnet}}.
#'
#' @export

  tune_glmnet <- function(x,
                          y,
                          fold_ids,
                          type.measure = c("default",
                                           "mse", "deviance",
                                           "class", "auc",
                                           "mae", "C"),
                          ...) {

    ## input control -------

    if(!is.matrix(x)) stop("'x' has to be a matrix.", call. = FALSE)

    ## validation of y and a more detailed validation of `fold_ids`
    ## is done by `cv.glmnet()` function

    if(!is.list(fold_ids)) stop("'fold_ids' has to be a list.", call. = FALSE)

    fold_format <- map_lgl(fold_ids, is.numeric)

    if(any(!fold_format)) {

      stop("'x' has to be list of numeric vectors.", call. = FALSE)

    }

    if(length(fold_ids) < 2) {

      stop(paste("Not enough repeats to proceed; 'fold_ids'",
                 "must have at least two elements"),
           call. = FALSE)

    }

    if(is.null(names(fold_ids))) {

      names(fold_ids) <- paste('repeat', seq_along(fold_ids), sep = '_')

    }

    type.measure <-
      match.arg(type.measure[1],
                c("default", "mse", "deviance",
                  "class", "auc", "mae", "C"))

    opt_fun <-
      switch(type.measure,
             default = min,
             mse = min,
             deviance = min,
             class = min,
             auc = max,
             mae = min,
             C = max)

    ## tuning models and stats ------

    cvm <- NULL
    rep <- NULL
    lambda <- NULL

    tune_models <-
      future_map(fold_ids,
                 function(id) cv.glmnet(x = x,
                                        y = y,
                                        type.measure = type.measure,
                                        fold.id = id, ...),
                 .options = furrr_options(seed = TRUE))

    lambda_stats <-
      map(tune_models,
          ~as_tibble(.x[c('lambda', 'cvm', 'cvup', 'cvlo')]))

    lambda_stats <- map(lambda_stats,
                        filter,
                        cvm == opt_fun(cvm))

    lambda_stats <-
      map2_dfr(lambda_stats, names(lambda_stats),
               ~mutate(.x[1, ], rep = .y))

    ## finding the optimal stats ---------

    ## it is possible that during multiple runs of the algorithm
    ## land at the same optimal lambda with different cost functions
    ## hence, the cost function is averaged

    best_lambda_tbl <-
      group_by(lambda_stats[c('lambda', 'cvm')],
               lambda)

    best_lambda_tbl <- summarise(best_lambda_tbl,
                                 lambda = mean(lambda),
                                 cvm = mean(cvm))

    best_lambda_tbl <- ungroup(best_lambda_tbl)

    best_lambda <- filter(best_lambda_tbl, cvm == opt_fun(cvm))

    best_lambda <- best_lambda$lambda[[1]]

    ## GLMnet model for the optimal lambda -------

    glm_model <- glmnet(x = x,
                        y = y,
                        lambda = best_lambda, ...)

    ## ouput ------

    cvTune(list(stats = arrange(lambda_stats, cvm),
                lambda = best_lambda,
                model = glm_model))

  }

# END -------
