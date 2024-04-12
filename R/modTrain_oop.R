# S3 methods for the `modTrain` class

# Class testing ------

#' Test for an `modTrain` class instance.
#'
#' @description
#' Tests is an object is an instance of the `modTrain` class.
#'
#' @param x an object.
#'
#' @return a logical value.
#'
#' @export

  is_modTrain <- function(x) inherits(x, 'modTrain')

# Appearance --------

#' Appearance of the `modTrain` class.
#'
#' @description
#' `print()` method for the `modTrain` class.
#'
#' @param x an instance of the `modTrain` class.
#' @param ... extra arguments, currently none.
#'
#' @return none, called for its side effects.
#'
#' @export

  print.modTrain <- function(x, ...) {

    cat(paste('modTrain object with n =', nrow(x$stats), 'models'))

  }

# Numbers of observations, explanatory variables and models -------

#' Numbers of observations, explanatory variables and models.
#'
#' @description
#' Computes numbers of complete observations, explanatory variables
#' and models in a `modTrain` object.
#'
#' @param object an instance of the `modTrain` class.
#' @param ... extra arguments, currently none.
#'
#' @return a data frame with the response name,
#' numbers of complete observations, numbers of modeling features, and
#' numbers of non-zero coeffcients of the regularized models.
#'
#' @export nobs.modTrain
#' @export

  nobs.modTrain <- function(object, ...) {

    stopifnot(is_modTrain(object))

    components(object, 'stats')[c('response',
                                  'n_observations',
                                  'n_features',
                                  'n_coefficients')]

  }

# Object components -----

#' Access components of a `modTrain` object.
#'
#' @description
#' Accesses components of the `modTrain` object.
#'
#' @details
#' The following are returned for the argument `type` specified as:
#'
#' * `models`: `cv.glmnet` models.
#'
#' * `stats`: model performance statistics.
#'
#' * `globals`: global parameter values used for construction of `glmnet` models.
#'
#' * `errors`: error messages for failed modeling tasks.
#'
#' @return see: `Details`.
#'
#' @param object a `modTrain` object.
#' @param type the requested component of the object, see `Details`.
#' @param ... extra arguments, currently none.
#'
#' @export components.modTrain
#' @export

  components.modTrain <- function(object,
                                  type = c('models',
                                           'stats',
                                           'globals',
                                           'errors'), ...) {

    stopifnot(is_modTrain(object))

    type <- match.arg(type[1],
                      c('models',
                        'stats',
                        'globals',
                        'errors'))

    object[[type]]

  }

# Object summary -------

#' Object summary.
#'
#' @description
#' Performance statistics of the models stored in a `modTrain` object.
#'
#' @return a data frame with model performance statistic in cross-validation
#' (loss function, optionally \eqn{\R^2}) and in the training data set
#' (optional: correlation coefficients, classification statistics),
#' along with numbers of observations, modeling variables,
#' and model coefficients.
#'
#' @param object a `modTrain` object.
#' @param ... extra argments, currently none.
#'
#' @export summary.modTrain
#' @export

  summary.modTrain <- function(object, ...) {

    components(object, 'stats')

  }

# Diagnostic plots ---------

#' Diagnostic plots.
#'
#' @description
#' Diagnostic plots of model performance measures in cross-validation and in
#' the train data set.
#'
#' @details
#' The following plots are generated:
#'
#' * violin plot of the cross-validated performance statistic used for
#' selection of the optimal \eqn{\lambda}, \eqn{\R^2}, correlation coefficients,
#' and classification metrics.
#'
#' * optional scatter plots of cross-validated \eqn{\R^2} or cross-validated
#' loss function and correlation coefficients (Pearson and Spearman, training
#' data, Gaussian and Poisson family models), cross-validated loss function and
#' classification metrics (accuracy and Cohen's \eqn{\kappa}, training data,
#' classification models).
#'
#' @return a list of `ggplot` plots, appearance can be easily customized
#' by the user.
#'
#' @param x a `modTrain` or `modHT object` object.
#' @param shape_alpha alpha for violin plot.
#' @param shape_color color of the violin plot.
#' @param point_alpha data point alpha.
#' @param point_color color of the data points.
#' @param point_size size of the data points.
#' @param point_wjitter jittering width of the data points.
#' @param point_hjitter jittering  height of the data points.
#' @param show_txt logical, should response names be displayed in the plots?
#' @param txt_size size of the text.
#' @param ... extra arguments, currently none.
#'
#' @export plot.modTrain
#' @export

  plot.modTrain <- function(x,
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
               family = components(x, 'globals')$family,
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

# Prediction and model coefficients --------

#' Predict outcomes or compute coefficients.
#'
#' @description
#' Predicts outcome, logits, coefficients and similar for a `modTrain` object.
#'
#' @details
#' For Gaussian family models with
#' types `type = 'link'` or `type = 'response'`, a matrix with predicted values
#' on the same scale as response is returned with user-specified transformations
#' and, optionally, standardization.
#' For other families, `type = 'link'` returns
#' a matrix of linear predictor score, while `type = 'response'` returns a
#' matrix of probabilities (binomial: one matrix, and multinomial regression:
#' list of matrices, with each element corresponding to one response variable)
#' or a matrix of predicted means (Poisson).
#' `type = 'class'` applies only to classification models and produces a matrix
#' of classes (labels). For `type = 'coefficients'` a matrix of coefficients for
#' the optimal \eqn{\lambda} is returned.
#' is returned.
#' See: \code{\link[glmnet]{predict.glmnet}} for details.
#' The `coef()` method returns a matrix of model coefficients for the optimal
#' \eqn{\lambda}.
#'
#' @return specific for the prediction type, see: `Details`.
#' If a \code{\link{modData}} object with multiple test data
#' (see: \code{\link{multi_process}}) sets will be
#' provided as `newdata`, a list of matrices is returned.
#'
#' @param object an `modTrain` object.
#' @param newdata a numeric matrix of explanatory factors with features in rows
#' and observations in columns. Alternatively, a \code{\link{modData}} object
#' created e.g. with the \code{\link{pre_process}} function.
#' @param type type of prediction, see: `Details`.
#' @param ... extra arguments passed to \code{\link[glmnet]{predict.glmnet}}
#' or to code{\link[glmnet]{coef.glmnet}}.
#'
#' @export predict.modTrain
#' @export

  predict.modTrain <- function(object,
                               newdata = NULL,
                               type = c('link',
                                        'response',
                                        'class',
                                        'coefficients'),
                               ...) {

    # entry control -------

    stopifnot(is_modTrain(object))

    type <-match.arg(type[1],
                     c('link',
                       'response',
                       'class',
                       'coefficients'))

    # entry control ------

    if(type == 'coefficients') return(coef(object, ...))

    if(is_modData(newdata)) newdata <- newdata$test

    err_txt <-
      "'newdata' has to be a numeric matrix or a list of numeric matrices."

    if(!is.list(newdata)) {

      if(!is.matrix(newdata)) stop(err_txt, call. = FALSE)

      if(!is.numeric(newdata)) stop(err_txt, call. = FALSE)

      newdata <- t(newdata)

    } else {

      if(any(!map_lgl(newdata, is.matrix))) stop(err_txt, call. = FALSE)

      if(any(!map_lgl(newdata, is.numeric))) stop(err_txt, call. = FALSE)

      newdata <- map(newdata, t)

    }

    dots <- list2(...)

    mod_family <- object$globals$family

    new_offset <- dots$newoffset

    ## predictions -------

    if(!is.list(newdata)) {

      preds <- pred_(models = object$models,
                     lambdas = object$stats$lambda,
                     newx = newdata,
                     type = type,
                     family = mod_family,
                     newoffset = new_offset, ...)

    } else {

      preds <-
        map(newdata,
            ~pred_(models = object$models,
                   lambdas = object$stats$lambda,
                   newx = .x,
                   type = type,
                   family = mod_family,
                   newoffset = new_offset, ...))

    }

    preds

  }

#' @rdname predict.modTrain
#' @export coef.modTrain
#' @export

  coef.modTrain <- function(object, ...) {

    stopifnot(is_modTrain(object))

    coef_mats <-
      map2(object$models,
           object$stats$lambda,
           ~coef(.x$glmnet.fit,
                 s = .y))

    responses <- object$stats$response

    mod_family <- object$globals$family

    if(mod_family != 'multinomial') {

      for(i in seq_along(coef_mats)) {

        colnames(coef_mats[[i]]) <- responses[[i]]

      }

      coefs <- reduce(map(coef_mats, as.matrix), cbind)

      return(coefs)

    }

    ## for multinomial models, a list of matrices is returned
    ## each element corresponds to a class

    coef_mats <- transpose(coef_mats)

    for(i in names(coef_mats)) {

      for(j in seq_along(coef_mats[[i]])) {

        colnames(coef_mats[[i]][[j]]) <- responses[[j]]

      }

      coef_mats[[i]] <- reduce(map(coef_mats[[i]], as.matrix), cbind)

    }

    return(coef_mats)

  }

# END -----
