# Classes to be used by the package

# modData class --------

#' `modData` class.
#'
#' @description
#' Creates an object of `modData` class, which can be used for construction of
#' regularized regression models, e.g. with the \code{\link{pre_process}}
#' function.
#'
#' @details
#' Contains the following elements:
#'
#' * `train`: a batch-adjusted numeric training matrix of selected
#' explanatory factors, features in rows and observations in columns.
#'
#' * `test`: a batch-adjusted numeric test matrix of selected
#' explanatory factors, features in rows and observations in columns.
#'
#' * `stats`: distribution statistics of all features prior to the batch
#' adjustment: mean, median, Gini index,
#' and an indicator if the feature was selected as a modeling variable.
#'
#' * `mean_limits`: cutoffs of mean feature values used for selection of the
#' modeling variables.
#'
#' * `median_limits`: cutoffs of median feature values used for selection of
#' the modeling variables.
#'
#' * `gini_limits`: cutoffs of Gini coefficient used for selection of the
#' modeling variables.
#'
#' @param x a list with elements specified in `Details`.
#' @param ... extra arguments, currently none.
#'
#' @return An instance of `modData` class with the following methods specified:
#' \code{\link{nobs.modData}}, \code{\link{components.modData}},
#' \code{\link{summary.modData}}, and \code{\link{plot.modData}}.
#'
#' @export

  modData <- function(x, ...) {

    ## input control: x as a list --------

    stopifnot(is.list(x))

    if(is.null(names(x))) stop("'x' has to be a named list.", call. = FALSE)

    if(any(!c('train', 'test', 'stats', 'mean_limits', 'median_limits', 'gini_limits') %in% names(x))) {

      stop("'x' has to have 'train', 'test', 'stats', 'mean_limits', 'median_limits', and 'gini_limits' elements.",
           call. = FALSE)

    }

    ## input controls: numeric matrices --------

    if(!is.list(x[['test']])) {

      check_lst <- x[c('train', 'test')]

    } else {

      check_lst <- c(list(train = x[['train']]),
                     x[['test']])

    }

    matrix_test <- map_lgl(check_lst, is.matrix)

    if(any(!matrix_test)) {

      stop("'train' and 'test' elements have to be numeric matrices",
           call. = FALSE)

    }

    num_test <- map_lgl(check_lst, is.numeric)

    if(any(!num_test)) {

      stop("'train' and 'test' elements have to be numeric matrices",
           call. = FALSE)

    }

    na_test <- map_lgl(check_lst, ~any(is.na(.x)))

    if(any(na_test)) {

      warning(paste("NAs in 'train' and 'test' elements.",
                    'Wrong transformation function?'),
              call. = FALSE)

    }

    common_features <- map(check_lst, rownames)

    common_features <- reduce(common_features, base::intersect)

    if(length(common_features) != nrow(check_lst[[1]])) {

      stop("Row names of 'train' and 'test' matrices have to be the same.",
           call. = FALSE)

    }

    ## input control: stats and limits --------

    stopifnot(is.data.frame(x$stats))

    if(!is.null(x$mean_limits)) stopifnot(is.numeric(x$mean_limits))

    if(!is.null(x$median_limits)) stopifnot(is.numeric(x$median_limits))

    if(!is.null(x$gini_limits)) stopifnot(is.numeric(x$gini_limits))

    ## output -------

    if(is.list(x[['test']])) {

      classes <- c('modMData', 'modData')

    } else {

      classes <- 'modData'

    }

    structure(x, class = classes)

  }

# modTrain class --------

#' `modTrain` class.
#'
#' @description
#' Creates an object of the `modTrain` class, which stores tuned
#' `glmnet` models and their basic performance stats.
#' It can be obtained by calling the
#' \code{\link{train.default}} or \code{\link{train.modData}} functions.
#'
#' @details
#' Contains the following elements:
#'
#' * `models`: a list of \code{\link[glmnet]{cv.glmnet}} models used for
#' tuning of \eqn{\lambda}, named after the columns of the `y` matrix.
#'
#' * `stats`: the optimal \eqn{lambda} parameter values, cross-validated
#' selection measures, and statistics of model performance in the
#' training data set.
#'
#' * `globals`: globals used for construction of the GLMNET models such as
#' `alpha` or `family`.
#'
#' * `errors`: a list of error messages for responses for which tuning/training
#' failed.
#'
#' Note the \code{\link{nobs.modTrain}}, \code{\link{components.modTrain}},
#' \code{\link{summary.modTrain}}, \code{\link{plot.modTrain}},
#' \code{\link{coef.modTrain}}, and \code{\link{predict.modTrain}} methods
#' defined for the class.
#'
#' @return an instance of the `modTrain` class.
#'
#' @param x a list with the elements specified in `Details`.
#' @param ... extra arguments, currently none.
#'
#' @export

  modTrain <- function(x, ...) {

    ## input control ------

    err_txt <- paste("'x' has to be a list with the 'models', 'stats'",
                     "'globals', and 'errors' elements.")

    if(!is.list(x)) stop(err_txt, call. = FALSE)

    if(is.null(names(x))) stop(err_txt, call. = FALSE)

    if(any(!c('models', 'stats', 'globals', 'errors') %in% names(x))) {

      stop(err_txt, call. = FALSE)

    }

    stopifnot(is.data.frame(x$stats))

    if(any(!c('response', 'lambda', 'loss_name', 'loss_criterion') %in% names(x$stats))) {

      stop(paste("'stats' data frame must have 'variable', 'lambda',",
                 "'loss_name', and 'loss_criterion' columns."),
           call. = FALSE)

    }

    if(any(!c('alpha', 'family', 'standardize') %in% names(x$globals))) {

      stop("'globals' must have 'alpha', 'family', and 'standardize' elements",
           call. = FALSE)

    }

    ## output --------

    structure(x, class = 'modTrain')

  }

# modHT class ---------

#' `modHT` class.
#'
#' @description
#' Creates an object of the `modHT` class, which stores model predictions in
#' the training data set along with pre-processing and modeling performance
#' statistics. It can be created e.g. with the \code{\link{fitHT}} function.
#'
#' @details
#' The object has the following components:
#'
#' * `predictions`: model predictions for the test data as a matrix or a list
#' of matrices.
#'
#' * `preprocess_stats`: information of modeling variables after
#' the pre-processing step, like variable name, median and mean value,
#' Gini index.
#'
#' * `model_stats`: metrics of the regularized linear models, such as loss
#' function used for \eqn{\lambda} tuning, \eqn{R^2}, correlations of the
#' observed and predicted modeling response in the training data (regression
#' models), as well as overall accuracy and Cohen's \eqn{\kappa} for comparison
#' of the observed and predicted class assignment (classification models).
#'
#' * `preprocess_globals`: pre-processing parameters specified by the user like
#' modeling variable selection cutoffs or the function used to transform the
#' explanatory variables.
#'
#' * `model_globals`: modeling parameters specified by the user such as modeling
#' family, the function used to transform the response variables, loss function,
#' mixing parameter alpha.
#'
#' * `prediction_globals`: parameters of predictions such as prediction type.
#'
#' * `errors`: an optional vector of error messages.
#'
#' Note the methods \code{\link{nobs.modHT}}, \code{\link{components.modHT}},
#' \code{\link{summary.modHT}}, and \code{\link{plot.modHT}} defined
#' for the class.
#'
#' @return an instance of the `modHT` class with components
#' described in `Details`.
#'
#' @param x a named list with elements described in `Details`.
#' @param ... extra arguments, currently none.
#'
#' @export

  modHT <- function(x, ...) {

    ## input controls

    elements <- c('predictions',
                  'preprocess_stats',
                  'model_stats',
                  'preprocess_globals',
                  'model_globals',
                  'prediction_globals',
                  'errors')

    err_txt <- paste("'x' has to be a list with the 'predictions',",
                     "'preprocess_stats', 'model_stats', 'preprocess_globals',",
                     "'model_globals', 'prediction_globals', and 'errors'.")

    if(!is.list(x)) stop(err_txt, call. = FALSE)

    if(is.null(names(x))) stop(err_txt, call. = FALSE)

    if(any(!elements %in% names(x))) {

      stop(err_txt, call. = FALSE)

    }

    ## output ------

    structure(x, class = 'modHT')

  }

# END ---------
