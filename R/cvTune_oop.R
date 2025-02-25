# S3 methods for `cvTune` objects that store results of cross-validated
# optimization of lambda

# class testing -------

#' Check for a `cvTune` object-
#'
#' @description
#' Tests if an object belong to \code{\link{cvTune}} class.
#'
#' @return a logical value.
#'
#' @param x an object.
#'
#' @export

 is_cvTune <- function(x) inherits(x, 'cvTune')

# Appearance of the class -------

#' Appearance of `cvTune` class.
#'
#' @description
#' `print()` display method for \code{\link{cvTune}} objects.
#'
#' @return nothing, called for its side effects.
#'
#' @param x an instance of `cvTune` class.
#' @param ... extra arguments passed to methods.
#'
#' @export

 print.cvTune <- function(x, ...) {

   stopifnot(is_cvTune(x))

   cat("GLMNET model tuning results. Best model:\n")
   print(x$model)

 }

# Observation number ---------

#' Number of observations, features and non-zero features in `cvTune` objects.
#'
#' @description
#' Retrieves the number of observations, all features, and features with
#' non-zero coefficients of the best performing `GLMNET` model.
#'
#' @return a data frame with number of observations,
#'
#' @param object an instance of the \code{\link{cvTune}} class.
#' @param ... extra arguments, currently none.
#'
#' @export

 nobs.cvTune <- function(object, ...) {

   observations <- NULL
   total_variables <- NULL
   non_zero_variables <- NULL

   tibble(observations = object$model$nobs,
          total_variables = object$model$dim[[1]],
          non_zero_variables = object$model$df)

 }

# Performance stats, summary --------

#' Performance statistics: summary of `cvTune` objects.
#'
#' @description
#' Retrieves locally optimal `lambda` values with cross-validated cost function
#' values from the \code{\link{cvTune}} model.
#'
#' @return
#' a sorted data frame with the globally optimal `lambda` value in the first
#' row and the following columns:
#'
#' * `rep`: identification of a repeat
#' * `cvm`: cross-validated cost function value for out-of-fold predictions in
#' particular repeats
#' * `cvup` and `cvlo`: upper bond of 95% confidence interval of `cvm`
#'
#' @param object an instance of \code{\link{cvTune}}.
#' @param ... extra arguments passed to methods; currently none.
#'
#' @export summary.cvTune
#' @export

  summary.cvTune <- function(object, ...) {

    stopifnot(is_cvTune(object))

    object$stats

  }

# Coefficients and predictions --------

#' Coeffcients and predictions by the tuned model in `cvTune` objects.
#'
#' @description
#' Methods `coef()` and `predict()` return coefficients and predictions of the
#' best performing `GLMNET` model.
#'
#' @return `coef()` returns a sparse one-column matrix with the model
#' coefficients as specified by \code{\link[glmnet]{coef.glmnet}}.
#' `predict()` returns linear predictor scores, predictions of the response or
#' class assignment as specified by \code{\link[glmnet]{predict.glmnet}}.
#'
#' @param object an instance of \code{\link{cvTune}}.
#' @param newx matrix of explanatory variables (a numeric matrix or an object
#' generated by \code{\link[stats]{model.matrix}}).
#' @param type of the prediction, as specified by
#' \code{\link[glmnet]{predict.glmnet}}.
#' @param ... additional arguments passed to \code{\link[glmnet]{coef.glmnet}}
#' or \code{\link[glmnet]{predict.glmnet}}.
#'
#' @export coef.cvTune
#' @export

  coef.cvTune <- function(object, ...) {

    stopifnot(is_cvTune(object))

    coef.glmnet(object = object$model,
                s = NULL, ...)

  }

#' @rdname coef.cvTune
#' @export predict.cvTune
#' @export

  predict.cvTune <- function(object,
                             newx,
                             type = c("link", "response",
                                      "coefficients", "nonzero", "class"),
                             ...) {

    stopifnot(is_cvTune(object))

    type <- match.arg(type[1],
                      c("link", "response",
                        "coefficients", "nonzero", "class"))

    predict(object = object$model,
            newx = newx,
            type = type,
            ...)

  }

# Diagnostic plots --------

#' Diagnostic plots for tuninig of `GLMNET` models.
#'
#' @description
#' Plots of the cost function values and `lambda` in repeats of the
#' cross-validated tuning of `GLMNET` models.
#'
#' @return a `ggplot` graphics, which can be customized by the user.
#'
#' @param x an instance of \code{\link{cvTune}} class.
#' @param palette a two element vector of names of colors for, respectively,
#' all data points and the best tune.
#' @param point_alpha data point alpha.
#' @param point_size size of the data points.
#' @param ... extra arguments, currently none.
#'
#' @export plot.cvTune
#' @export

  plot.cvTune <- function(x,
                          palette = c('steelblue', 'orangered3'),
                          point_alpha = 0.75,
                          point_size = 2, ...) {

    ## input control --------

    stopifnot(is_cvTune(x))

    if(length(palette) < 2) {

      stop("At least two colors provide as 'palette' are required.",
           call. = FALSE)

    }

    palette <- palette[1:2]

    stopifnot(is.numeric(point_alpha))
    stopifnot(is.numeric(point_size))

    ## plotting --------

    best <- NULL
    lambda <- NULL
    cvm <- NULL

    stat_data <-
      mutate(x$stats,
             best = ifelse(lambda == x$lambda,
                           'yes', 'no'))

    qc_plot <-
      ggplot(stat_data,
             aes(x = lambda,
                 y = cvm,
                 color = best)) +
      geom_point(shape = 16,
                 size = point_size,
                 alpha = point_alpha) +
      scale_color_manual(values = palette) +
      guides(color = 'none') +
      labs(title = 'Tuning of \u03BB',
           subtitle = paste('best \u03BB =', signif(x$lambda, 2)),
           x = '\u03BB',
           y = 'cost function, repeated CV')

    qc_plot

  }

# END ------
