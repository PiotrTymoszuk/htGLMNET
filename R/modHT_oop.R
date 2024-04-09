# S3 methods for the `modHT` class

# class inheritance -------

#' Test for an `modHT` class instance.
#'
#' @description
#' Tests is an object is an instance of the `modHT` class.
#'
#' @param x an object.
#'
#' @return a logical value.
#'
#' @export

  is_modHT <- function(x) inherits(x, 'modHT')

# Appearance ---------

#' Appearance of the `modHT` class.
#'
#' @description
#' `print()` method for the `modHT` class.
#'
#' @param x an instance of the `modHT` class.
#' @param ... extra arguments, currently none.
#'
#' @return none, called for its side effects.
#'
#' @export

  print.modHT <- function(x, ...) {

    stopifnot(is_modHT(x))

    cat(paste('modHT object with n =', nrow(x$model_stats), 'predictions'))

  }

# Numbers of observations, explanatory variables and models -------

#' Numbers of observations, explanatory variables and models.
#'
#' @description
#' Computes numbers of complete observations, explanatory variables
#' and models in a `modHT` object.
#'
#' @param object an instance of the `modHT` class.
#' @param ... extra arguments, currently none.
#'
#' @return a data frame with the response name,
#' numbers of complete observations, numbers of modeling features, and
#' numbers of non-zero coeffcients of the regularized models.
#'
#' @export nobs.modHT
#' @export

  nobs.modHT <- function(object, ...) {

    stopifnot(is_modHT(object))

    components(object, 'model_stats')[c('response',
                                        'n_observations',
                                        'n_features',
                                        'n_coefficients')]

  }

# Object components -----

#' Access components of a `modHT` object.
#'
#' @description
#' Accesses components of the `modHT` object.
#'
#' @details
#' The following are returned for the argument `type` specified as:
#'
#' * `predictions`: the prediction matrix or a list of prediction matrices.
#'
#' * `preprocess_stats`: feature distribution statistics after batch correction.
#'
#' * `model_stats`: model performance statistics.
#'
#' * `globals`: global parameter values used for pre-processing, modeling
#' and predictions.
#'
#' * `errors`: error messages for failed modeling tasks.
#'
#' @return see: `Details`.
#'
#' @param object a `modHT` object.
#' @param type the requested component of the object, see `Details`.
#' @param ... extra arguments, currently none.
#'
#' @export components.modHT
#' @export

  components.modHT <- function(object,
                                  type = c('predictions',
                                           'preprocess_stats',
                                           'model_stats',
                                           'globals',
                                           'errors'), ...) {

    stopifnot(is_modHT(object))

    type <- match.arg(type[1],
                      c('predictions',
                        'preprocess_stats',
                        'model_stats',
                        'globals',
                        'errors'))

    if(type == 'globals') return(object[c('preprocess_globals',
                                          'model_globals',
                                          'prediction_globals')])

    object[[type]]

  }

# Object summary -------

#' Object summary.
#'
#' @description
#' Performance statistics of the models used for predictions in a `modHT` object.
#'
#' @return a data frame with model performance statistic in cross-validation
#' (loss function, optionally \eqn{\R^2}) and in the training data set
#' (optional: correlation coefficients, classification statistics),
#' along with numbers of observations, modeling variables, and
#' model coefficients.
#'
#' @param object a `modHT` object.
#' @param ... extra argments, currently none.
#'
#' @export summary.modHT
#' @export

  summary.modHT <- function(object, ...) {

    components(object, 'model_stats')

  }

# Plotting -----------

#' @rdname plot.modTrain
#' @export plot.modHT
#' @export

  plot.modHT <- function(x,
                         shape_alpha = 0.5,
                         shape_color = 'steelblue',
                         point_alpha = 0.75,
                         point_color = shape_color,
                         point_size = 2,
                         point_wjitter = 0.0,
                         point_hjitter = 0.0,
                         show_txt = TRUE,
                         txt_size = 2.75, ...) {

    ## input control and plotting data -------

    stopifnot(is.numeric(shape_alpha))
    stopifnot(is.numeric(point_alpha))
    stopifnot(is.numeric(point_size))
    stopifnot(is.numeric(point_wjitter))
    stopifnot(is.numeric(point_hjitter))
    stopifnot(is.numeric(point_size))

    stopifnot(is.logical(show_txt))

    if(!show_txt) txt_color <- NA else txt_color <- 'black'

    ## other arguments are checked by downstream functions

    plot_stats(data = summary(x),
               family = x$model_globals$family,
               shape_alpha = shape_alpha,
               shape_color = shape_color,
               point_alpha = point_alpha,
               point_color = point_color,
               point_size = point_size,
               point_wjitter = point_wjitter,
               point_hjitter = point_hjitter,
               show_txt = show_txt,
               txt_size = txt_size,
               txt_color = txt_color)

  }

# END -------
