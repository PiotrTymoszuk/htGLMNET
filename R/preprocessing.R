# Functions used for pre-processing

# Pre-processing function --------

#' Pre-process training and test data sets of explanatory variables.
#'
#' @description
#' A function for selection of common variable features in the training and
#' and test sets of explanatory variables (e.g. gene expression counts or
#' protein reads) and batch correction. This tool streamlines the pre-processing
#' step of model development and offers few quality control measures that
#' allow the user to control the pre-processed data sets prior to modeling.
#'
#' @details
#' Technically, the function identifies common features in the `train` and
#' `test` numeric matrices of explanatory data, computes their distribution
#' metrics (median, mean, Gini index), and selects the modeling variables based
#' on cutoffs of mean, median, or Gini index provided by the user.
#' Finally, the numeric matrices with the selected modeling features are
#' subjected to batch adjustment with the ComBat algorithm
#' (see: \code{\link[sva]{ComBat}} for details).
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
#' @param test a numeric matrix of test explanatory variables in the same format
#' as above.
#' @param mean_quantile a numeric that specifies the quantile cutoff
#' of mean value of the features to be selected. Default is `NULL`, which means
#' that no filtering by mean feature value will be applied.
#' @param median_quantile as above, a quantile of median used to select
#' the features. Defaults to `NULL`, which implicates that no median feature
#' filtering will be applied.
#' @param gini_quantile as above, a quantile of Gini coefficient used to select
#' the modeling features of sufficient variability. If `NULL` (default), no
#' filtering of modeling features by Gini coefficient will be applied.
#' @param trans_fun a function used to transform the output, `log` by default.
#' @param ... extra arguments passed to \code{\link[sva]{ComBat}}.
#'
#' @return
#' An instance of the `modData` class with the following components:
#'
#' * `train`: a batch-adjusted numeric training matrix of selected
#' explanatory factors, features in rows and observations in columns.
#'
#' * `test`: a batch-adjusted numeric test matrix of selected
#' explanatory factors, features in rows and observations in columns.
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

    ## common features and feature selection --------

    common_features <- map(input_data, rownames)

    common_features <- reduce(common_features, base::intersect)

    input_data <- map(input_data, ~.x[common_features, ])

    gene_stats <- map(input_data, row_stats)

    ## filtering out completely invariant genes

    zeroVar <- NULL
    variable <- NULL

    gene_stats <- map(gene_stats,
                      filter,
                      !zeroVar)

    gene_stats <- map(gene_stats,
                      ~.x[c('variable', 'median', 'mean', 'gini_coef')])

    gene_selection <- gene_stats

    mean_limits <- NULL
    median_limits <- NULL
    gini_limits <- NULL

    if(!is.null(mean_quantile)) {

      mean_limits <- map_dbl(gene_selection,
                             ~quantile(.x$mean,
                                       probs = mean_quantile,
                                       na.rm = TRUE))

      gene_selection <-
        map2(gene_selection, mean_limits,
             ~filter(.x, mean >= .y))

    }

    if(!is.null(median_quantile)) {

      median_limits <- map_dbl(gene_selection,
                               ~quantile(.x$median,
                                         probs = median_quantile,
                                         na.rm = TRUE))

      gene_selection <-
        map2(gene_selection, median_limits,
             ~filter(.x, median >= .y))

    }

    if(!is.null(gini_quantile)) {

      gini_limits <- map_dbl(gene_selection,
                             ~quantile(.x$gini_coef,
                                       probs = gini_quantile,
                                       na.rm = TRUE))

      gene_selection <-
        map2(gene_selection, gini_limits,
             ~filter(.x, gini_coef >= .y))

    }

    gene_selection <- map(gene_selection, ~.x$variable)

    gene_selection <- reduce(gene_selection, base::intersect)

    modeling_variable <- NULL

    gene_stats <- map(gene_stats,
                      mutate,
                      modeling_variable = ifelse(variable %in% gene_selection,
                                                 'yes', 'no'))

    input_data <- map(input_data, ~.x[gene_selection, ])

    ## batch adjustment ---------

    batch_ids <- map(input_data, colnames)

    batch_vct <- c(rep('train', ncol(input_data[[1]])),
                   rep('test', ncol(input_data[[2]])))

    cmm_data <- cbind(input_data[[1]], input_data[[2]])

    adj_data <- ComBat(dat = cmm_data,
                       batch = batch_vct, ...)

    adj_data <- trans_fun(adj_data)

    adj_data <- map(batch_ids, ~adj_data[, .x])

    ## output ---------

    data_set <- NULL

    gene_stats <-
      map2_dfr(gene_stats, names(gene_stats),
               ~mutate(.x, data_set = .y))

    gene_stats <- relocate(gene_stats, data_set)

    modData(list(train = adj_data$train,
                 test = adj_data$test,
                 stats = gene_stats,
                 mean_limits = mean_limits,
                 median_limits = median_limits,
                 gini_limits = gini_limits))

  }

# END -------
