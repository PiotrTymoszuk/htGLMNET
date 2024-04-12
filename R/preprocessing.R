# Functions used for pre-processing

# Pre-processing function --------

#' Pre-process training and test data sets of explanatory variables.
#'
#' @description
#' Functions for selection of common variable features in the training and
#' and test sets of explanatory variables (e.g. gene expression counts or
#' protein reads) and batch correction. These tools streamline the
#' pre-processing step of model development and offer few quality control
#' measures that allow the user to control the pre-processed data sets
#' prior to modeling.
#' The function `pre_process()` is designed to handle a pair of training
#' and test matrices and its output may be passed directly to the modeling tool
#' `train()`.
#' The function `multi_process()`, can handle any list of numeric matrices.
#'
#' @details
#' Technically, the functions identify common features in the `train` and
#' `test` numeric matrices of explanatory data.
#' Next, distribution metrics of such common features (median, mean, Gini index)
#' are computed, and modeling variables are selected based
#' on cutoffs of mean, median, or Gini index provided by the user.
#' Finally, the numeric matrices with the selected modeling features are
#' subjected to batch adjustment with the ComBat algorithm
#' (see: \code{\link[sva]{ComBat}} for details).
#' The function `multi_process()` will run in parallel, when a
#' `furrr`-compatible backend is declared.
#'
#' @references
#' Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD.
#' The sva package for removing batch effects and other unwanted variation
#' in high-throughput experiments.
#' Bioinformatics (2012) 28:882. doi:10.1093/BIOINFORMATICS/BTS034
#'
#' @param train a numeric matrix of training explanatory variables with
#' features (e.g. genes) in rows and observations in columns.
#' Note: it is recommended that it contains only positive values, if the
#' modeling features are selected with a Gini index cutoff.
#' @param test a numeric matrix of test explanatory variables
#' (function `pre_process()`) or a list of numeric matrices in the same format
#' as `train`.
#' @param mean_quantile a numeric that specifies the quantile cutoff
#' of mean value of the features to be selected. Default is `NULL`, which means
#' that no filtering by mean feature value will be applied.
#' @param median_quantile as above, a quantile of median used to select
#' the features. Defaults to `NULL`, which implicates that no median feature
#' filtering will be applied.
#' @param gini_quantile as above, a quantile of Gini coefficient used to select
#' the modeling features of sufficient variability. If `NULL` (default), no
#' filtering of modeling features by Gini coefficient will be applied.
#' @param trans_fun a function used to transform the output,
#' `identity` by default.
#' @param ... extra arguments passed to \code{\link[sva]{ComBat}}.
#'
#' @return
#' An instance of the `modData` class with the following components:
#'
#' * `train`: a batch-adjusted numeric training matrix of selected
#' explanatory factors, features in rows and observations in columns.
#'
#' * `test`: for `pre_process()`: a batch-adjusted numeric test matrix of
#' selected explanatory factors, features in rows and observations in columns.
#' For `multi_process()`: a list of batch-adjusted numeric test matrices of the
#' selected modeling explanatory features.
#'
#' * `stats`: distribution statistics of all features prior to the batch
#' correction: mean, median, Gini index,
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
#' Note the \code{\link{nobs.modData}}, \code{\link{components.modData}},
#' \code{\link{summary.modData}}, and \code{\link{plot.modData}} methods
#' defined for the class.
#'
#' @export

  pre_process <- function(train,
                          test,
                          mean_quantile = NULL,
                          median_quantile = NULL,
                          gini_quantile = NULL,
                          trans_fun = identity,
                          ...) {

    ## input control ---------

    input_data <- list(train = train,
                       test = test)

    if(any(!map_lgl(input_data, is.matrix))) {

      stop('Train and test have to be numeric matrices', call. = FALSE)

    }

    if(any(!map_lgl(input_data, is.numeric))) {

      stop('Train and test have to be numeric matrices', call. = FALSE)

    }

    if(any(map_lgl(input_data, ~is.null(rownames(.x))))) {

      stop('Train and test have to have row and column names', call. = FALSE)

    }

    if(any(map_lgl(input_data, ~is.null(colnames(.x))))) {

      stop('Train and test have to have row and column names', call. = FALSE)

    }

    mean_quantile <- mean_quantile[1]
    median_quantile <- median_quantile[1]
    gini_quantile <- gini_quantile[1]

    if(!is.null(mean_quantile)) {

      stopifnot(is.numeric(mean_quantile))

      if(mean_quantile < 0 | mean_quantile > 1) {

        stop("'mean_quantile' must be in [0, 1] range", call. = FALSE)

      }

    }

    if(!is.null(median_quantile)) {

      stopifnot(is.numeric(median_quantile))

      if(median_quantile < 0 | median_quantile > 1) {

        stop("'median_quantile' must be in [0, 1] range", call. = FALSE)

      }

    }

    if(!is.null(gini_quantile)) {

      stopifnot(is.numeric(gini_quantile))

      if(gini_quantile < 0 | gini_quantile > 1) {

        stop("'gini_quantile' must be in [0, 1] range", call. = FALSE)

      }

    }

    stopifnot(is.function(trans_fun))

    if(any(map_lgl(input_data, ~any(.x < 0)))) {

      if(!is.null(gini_quantile)) {

        warning(paste('The input train or test data contain negative values.',
                      'Selection of modeling features with a Gini index cutoff',
                      'is not recommended.'))

      }

    }

    ## common feature identifiers throw an error,
    ## not compatible with ComBat!

    common_observations <- map(input_data, colnames)

    common_observations <- reduce(common_observations, intersect)

    if(length(common_observations) > 0) {

      stop('Observation names must be unique.', call. = FALSE)

    }

    ## common features, feature selection and batch adjustment ------

    proc_res <- x_mats(x = input_data,
                       mean_quantile = mean_quantile,
                       median_quantile = median_quantile,
                       gini_quantile = gini_quantile,
                       trans_fun = trans_fun, ...)

    ## output ---------

    data_set <- NULL

    gene_stats <-
      map2_dfr(proc_res$stats, names(proc_res$stats),
               ~mutate(.x, data_set = .y))

    gene_stats <- relocate(gene_stats, data_set)

    modData(list(train = proc_res$x$train,
                 test = proc_res$x$test,
                 stats = gene_stats,
                 mean_limits = proc_res$mean_limits,
                 median_limits = proc_res$median_limits,
                 gini_limits = proc_res$gini_limits))

  }

#' @rdname pre_process
#' @export

  multi_process <- function(train,
                            test,
                            mean_quantile = NULL,
                            median_quantile = NULL,
                            gini_quantile = NULL,
                            trans_fun = identity,
                            ...) {

    ## input control ---------

    if(!is.list(test)) stop("'test' has to be a named list.", call. = FALSE)

    if(is.null(names(test))) stop("'test' has to be a named list.", call. = FALSE)

    input_data <- c(list(train = train), test)

    if(any(!map_lgl(input_data, is.matrix))) {

      stop('Train and test have to store numeric matrices', call. = FALSE)

    }

    if(any(!map_lgl(input_data, is.numeric))) {

      stop('Train and test have to store numeric matrices', call. = FALSE)

    }

    if(any(map_lgl(input_data, ~is.null(rownames(.x))))) {

      stop('Train and test matrices have to have row and column names',
           call. = FALSE)

    }

    if(any(map_lgl(input_data, ~is.null(colnames(.x))))) {

      stop('Train and test matrices have to have row and column names',
           call. = FALSE)

    }

    mean_quantile <- mean_quantile[1]
    median_quantile <- median_quantile[1]
    gini_quantile <- gini_quantile[1]

    if(!is.null(mean_quantile)) {

      stopifnot(is.numeric(mean_quantile))

      if(mean_quantile < 0 | mean_quantile > 1) {

        stop("'mean_quantile' must be in [0, 1] range", call. = FALSE)

      }

    }

    if(!is.null(median_quantile)) {

      stopifnot(is.numeric(median_quantile))

      if(median_quantile < 0 | median_quantile > 1) {

        stop("'median_quantile' must be in [0, 1] range", call. = FALSE)

      }

    }

    if(!is.null(gini_quantile)) {

      stopifnot(is.numeric(gini_quantile))

      if(gini_quantile < 0 | gini_quantile > 1) {

        stop("'gini_quantile' must be in [0, 1] range", call. = FALSE)

      }

    }

    stopifnot(is.function(trans_fun))

    if(any(map_lgl(input_data, ~any(.x < 0)))) {

      if(!is.null(gini_quantile)) {

        warning(paste('The input train or test data contain negative values.',
                      'Selection of modeling features with a Gini index cutoff',
                      'is not recommended.'))

      }

    }

    ## non-unique observation identifiers should return an error
    ## because of later problems with ComBat

    common_observations <- map(input_data, colnames)

    common_observations <- reduce(common_observations, intersect)

    if(length(common_observations) > 0) {

      stop('Observation names must be unique.', call. = FALSE)

    }

    ## common features, feature selection and batch adjustment ------

    proc_res <- x_mats(x = input_data,
                       mean_quantile = mean_quantile,
                       median_quantile = median_quantile,
                       gini_quantile = gini_quantile,
                       trans_fun = trans_fun, ...)

    ## output ---------

    data_set <- NULL

    gene_stats <-
      map2_dfr(proc_res$stats, names(proc_res$stats),
               ~mutate(.x, data_set = .y))

    gene_stats <- relocate(gene_stats, data_set)

    modData(list(train = proc_res$x$train,
                 test = proc_res$x[-1],
                 stats = gene_stats,
                 mean_limits = proc_res$mean_limits,
                 median_limits = proc_res$median_limits,
                 gini_limits = proc_res$gini_limits))

  }

# END -------
