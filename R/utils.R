# Non exported utilities

# Pre-processing -------

#' Pre-process a list of numeric matrices.
#'
#' @description
#' Pre-processing of a list of numeric matrices as described for
#' \code{\link{pre_process}}.
#'
#' @details
#' For internal use only!
#'
#' @return
#' A list with two elements:
#'
#' * `x`: a list of numeric matrices adjusted for bacth effects and containing
#' only the selected modeling features.
#'
#' * `stats`: feature distributions statistics (median, mean, Gini index) of
#' all features shared by the input data.
#'
#' * `mean_limits`, `median_limits`, `gini_limits`: numeric vectors of cutoffs
#' of mean, median, and Gini index of explanatory features used for selection
#' of the modeling variables.
#'
#' @param x a list of numeric matrices. Modeling features are provided in rows,
#' observations are in columns.#'
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

  x_mats <- function(x,
                     mean_quantile = NULL,
                     median_quantile = NULL,
                     gini_quantile = NULL,
                     trans_fun = identity, ...) {

    ## the input control is done by the upstream functions.

    ## common features and feature selection --------

    common_features <- map(x, rownames)

    common_features <- reduce(common_features, base::intersect)

    if(length(common_features) == 0) {

      stop('There are no shared modeling features.', call. = FALSE)

    }

    x <- map(x, ~.x[common_features, ])

    ## implementation of palatalization for larger inputs

    if(length(x) > 2) {

      gene_stats <- future_map(x,
                               row_stats,
                               .options = furrr_options(seed = TRUE))

    } else {

      gene_stats <- map(x, row_stats)

    }

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

    x <- map(x, ~.x[gene_selection, ])

    ## batch adjustment ---------

    batch_ids <- map(x, colnames)

    batch_vct <-
      map2(names(batch_ids),
           map_dbl(batch_ids, length),
           ~rep(.x, .y))

    batch_vct <- reduce(batch_vct, c)

    cmm_data <- reduce(x, cbind)

    adj_data <- ComBat(dat = cmm_data,
                       batch = batch_vct, ...)

    adj_data <- trans_fun(adj_data)

    adj_data <- map(batch_ids, ~adj_data[, .x])

    ## output ---------

    list(x = adj_data,
         stats = gene_stats,
         mean_limits = mean_limits,
         median_limits = median_limits,
         gini_limits = gini_limits)

  }

# Distribution stats of a list of numeric matrices ------

#' Distribution statistics for a list of numeric matrices.
#'
#' @description
#' An internal tool used by `summary()` methods, which computes
#' numbers of features, observations, missing data poins and distribution
#' statistics for a list of numeric matrices.
#'
#' @return a data frame with the statistics.
#'
#' @param x a list of numeric matrices.

  mat_stats <- function(x) {

    map2_dfr(x, names(x),
             ~tibble(data_set = .y,
                     n_features = nrow(.x),
                     n_observations = ncol(.x),
                     n_missing = sum(is.na(.x)),
                     mean = mean(.x, na.rm = TRUE),
                     sd = sd(.x, na.rm = TRUE),
                     median = median(.x, na.rm = TRUE),
                     q25 = quantile(.x, 0.25, na.rm = TRUE),
                     q75 = quantile(.x, 0.75, na.rm = TRUE),
                     q025 = quantile(.x, 0.025, na.rm = TRUE),
                     q975 = quantile(.x, 0.975, na.rm = TRUE)))

  }

# Fitting of Gaussian and Poisson family models ------

#' Serial fit of models.
#'
#' @description
#' Returns a list of `cv.glmnet` class models given a list of responses and
#' a common matrix of explanatory variables. `fit_regression()` is intended for
#' tuning and training regression models, `fit_classification()` is designed
#' for classification models (binomial and multinomial).
#'
#' @details
#' Intended for internal use. Will run in parallel, when a `furrr`-compatible
#' backend is declared.
#' Classification models are handled in a special way involving call evaluation,
#' because direct passing of additional arguments
#' to \code{\link[glmnet]{cv.glmnet}} and \code{\link[glmnet]{glmnet}} via
#' dots causes C++ errors. This is likely attributed to an unexpected behavior
#' of the `glmnet` `lognet()` function.
#'
#' @param y_lst a list of named vectors with the responses.
#' @param x a matrix with explanatory factors.
#' @param safe logical, if `TRUE`, the procedures is run with purrr's
#' `safely()`: correct evaluations are stored in `cv.models`,
#' errors are stored in `errors`.
#' @param ... other arguments passed to \code{\link[glmnet]{cv.glmnet}} and
#' \code{\link[glmnet]{glmnet}}, such as `alpha`, `family` or `type.measure`.
#'
#' @return A list with the following elements:
#'
#' * `cv_models`: a list `cv.glmnet` class models.
#'
#' * `errors`: a list of error messages.

  fit_regression <- function(y_lst, x, safe = FALSE, ...) {

    ## input control: done by the upstream function

    if(!safe) {

      cv_models <-
        future_map(y_lst,
                   ~cv.glmnet(x = x[names(.x), ],
                              y = .x, ...),
                   .options = furrr_options(seed = TRUE))

      errors <- NULL

    } else {

      cv_models <-
        future_map(y_lst,
                   ~safely(cv.glmnet)(x = x[names(.x), ],
                                      y = .x, ...),
                   .options = furrr_options(seed = TRUE))

      errors <- compact(map(cv_models, ~.x$error))

      cv_models <- compact(map(cv_models, ~.x$result))

    }

    list(cv_models = cv_models,
         errors = errors)
  }

#' @rdname fit_regression

  fit_classification <- function(y_lst, x, safe = FALSE, ...) {

    ## input control is performed by an upstream function

    dots <- list2(...)

    if(safe) fun <- safely(cv.glmnet) else fun <- cv.glmnet

    calls <-
      map(y_lst,
          ~call2(.fn = fun,
                 x = x[names(.x), ],
                 y = .x,
                 !!!dots))

    cv_models <-
      future_map(calls,
                 eval,
                 .options = furrr_options(seed = TRUE))

    if(!safe) {

      errors <- NULL

    } else {

      errors <- compact(map(cv_models, ~.x$error))

      cv_models <- compact(map(cv_models, ~.x$result))

    }

    list(cv_models = cv_models,
         errors = errors)

  }

# Model evaluation in the training data set --------

#' Evaluation of performance of `cv.glmnet` models in the training data set.
#'
#' @description
#' For regression models, i.e. Gaussian and Poisson family, the function
#' computes correlation coefficients between the predicted and observed outcome.
#' For classification models overall accuracy and Cohen's \eqn{\kappa} are
#' computed for comparison of the actual and predicted class assignment.
#'
#' @details
#' Overall accuracy and Cohen's \eqn{\kappa} are computed with
#' \code{\link[caret]{multiClassSummary}}.
#'
#' @return `evaluate_regression()` returns a list of two numeric vectors, with
#' Pearson's and Spearman's correlation coefficients, respectively.
#' `evaluate_classification()` returns a list of two numeric vectors ,
#' with overall accuracy and Cohen's \eqn{\kappa}, respectively.
#'
#' @param model_lst a list of `cv.glmnet` models.
#' @param lambda a vector with lambda values, for which the predictions will
#' be made.
#' @param y_lst a list of named vectors with the observed responses.
#' @param x a numeric matrix of explanatory factors.

  evaluate_regression <- function(model_lst, lambda, y_lst, x) {

    ## input control is done by the main function

    train_preds <- map2(model_lst,
                        lambda,
                        ~predict(.x,
                                 newx = x,
                                 type = 'link',
                                 s = .y,
                                 newoffset = .x$glmnet.fit$offset))

    train_preds <- map2(train_preds, y_lst,
                        ~cbind(predicted = as.numeric(.x[names(.y), 1]),
                               observed = .y))

    pearson <- map(train_preds,
                   cor,
                   method = 'pearson',
                   use = 'complete.obs')

    spearman <-  map(train_preds,
                     cor,
                     method = 'spearman',
                     use = 'complete.obs')

    list(pearson = map_dbl(pearson, ~.x[1, 2]),
         spearman = map_dbl(spearman, ~.x[1, 2]))

  }

#' @rdname evaluate_regression

  evaluate_classification <- function(model_lst, lambda, y_lst, x) {

    # input control is done by the main function

    train_preds <- map2(model_lst,
                        lambda,
                        ~predict(.x,
                                 newx = x,
                                 type = 'class',
                                 s = .y,
                                 newoffset = .x$glmnet.fit$offset))

    train_preds <- map(train_preds, ~set_names(as.numeric(.x), rownames(x)))

    train_preds <-
      map2(train_preds,
           y_lst,
           ~data.frame(pred = .x[names(.y)],
                       obs = .y))

    obs_levs <- map(train_preds, ~levels(factor(.x$obs)))

    train_preds <-
      map2(train_preds,
           obs_levs,
           ~mutate(.x,
                   pred = factor(pred, levels = .y),
                   obs = factor(obs, levels = .y)))

    train_preds <- map(train_preds,
                       ~filter(.x, complete.cases(.x)))

    train_stats <-
      map2(train_preds,
           obs_levs,
           ~multiClassSummary(data = .x, lev = .y))

    list(accuracy = map_dbl(train_stats, ~.x[['Accuracy']]),
         kappa = map_dbl(train_stats, ~.x[['Kappa']]))

  }

# Diagnostic plots ---------

#' Model diagnostic plots.
#'
#' @description
#' A function intended for internal use that takes a data frame with plotting
#' stats and a modeling family name and returns family-specific violin and
#' acatter plots.
#'
#' @return a list of `ggplot` objects.
#'
#' @param data a data frame with model performance metrics.
#' @param family name of the modeling family.
#' @param shape_alpha alpha for violin plot.
#' @param shape_color color of the violin plot.
#' @param point_alpha data point alpha.
#' @param point_color color of the data points.
#' @param point_size size of the data points.
#' @param point_wjitter jittering width of the data points.
#' @param point_hjitter jittering  height of the data points.
#' @param show_txt logical, should response names be displayed in the plots?
#' @param txt_size size of the text.
#' @param txt_color color of the text.

  plot_stats <- function(data,
                         family,
                         shape_alpha = 0.5,
                         shape_color = 'steelblue',
                         point_alpha = 0.75,
                         point_color = shape_color,
                         point_size = 2,
                         point_wjitter = 0.0,
                         point_hjitter = 0.0,
                         show_txt = TRUE,
                         txt_size = 2.75,
                         txt_color = 'black') {

    ## entry control is done by upstream functions

    ## plot axis labels ---------

    ax_title <- data$loss_name[1]

    ax_title <-
      switch(ax_title,
             default = 'default loss function',
             mse = 'MSE',
             deviance = 'deviance',
             class = 'classification error',
             auc = 'AUC',
             mae = 'MAE',
             C = 'C-index')

    if(family == 'gaussian' & ax_title == 'default') {

      ax_title <- 'MSE'

    }

    ax_title <- paste('loss function,', ax_title)

    ## violin plots -------

    loss_value <- NULL
    response <- NULL

    if(family %in% c('gaussian', 'poisson')) {

      violin_list <-
        list(x_var = c('CV', 'CV', 'training', 'traininig'),
             y_var = c('loss_value', 'rsq_oof', 'pearson', 'spearman'),
             plot_title = c('Loss function, cross-validation',
                            'Explained variance, cross-validation',
                            'Correlation with the outcome, training data',
                            'Correlation with the outcome, training data'),
             y_lab = c(ax_title,
                       'R\u00B2',
                       "Pearson's r",
                       "Spearman's \u03C1"))

      violin_names <- c('loss_function', 'rsq', 'pearson', 'spearman')

      if(family == 'poisson') {

        violin_names <- violin_names[-2]

        violin_list <- map(violin_list, ~.x[-2])

      }

    } else {

      violin_list <-
        list(x_var = c('CV', 'training', 'traininig'),
             y_var = c('loss_value', 'accuracy', 'kappa'),
             plot_title = c('Loss function, cross-validation',
                            'Overall accuracy, training data',
                            'Overall accuracy, training data'),
             y_lab = c(ax_title,
                       'overall accuracy',
                       "Cohen's \u03BA"))

      violin_names <- c('loss_function', 'accuracy', 'kappa')

    }

    violin_plots <-
      pmap(violin_list,
           function(x_var, y_var, plot_title, y_lab) ggplot(data,
                                                            aes(x = x_var,
                                                                y = .data[[y_var]])) +
             geom_violin(alpha = shape_alpha,
                         fill = shape_color) +
             geom_point(shape = 21,
                        fill = point_color,
                        size = point_size,
                        alpha = point_alpha,
                        position = position_jitter(width = 0.1,
                                                   height = 0)) +
             geom_text_repel(aes(label = response),
                             size = txt_size,
                             color = txt_color) +
             theme_bw() +
             theme(axis.title.x = element_blank()) +
             labs(title = plot_title,
                  y = y_lab))

    violin_plots <- set_names(violin_plots, violin_names)

    ## scatter plots, family specific ------

    if(family %in% c('gaussian', 'poisson')) {

      scatter_list <-
        list(x_var = c(rep('loss_value', 2), rep('rsq_oof', 2)),
             y_var = rep(c('pearson', 'spearman'), 2),
             x_lab = c(rep(paste0(ax_title, ', cross-validation'), 2),
                       rep('R\u00B2, cross-validation', 2)),
             y_lab = rep(c('r, training', '\u03C1, training'), 2))

      scatter_names <- c('loss_pearson', 'loss_spearman',
                         'rsq_pearson', 'rsq_spearman')

      if(family == 'poisson') {

        scatter_list <- map(scatter_list, ~.x[1:2])

        scatter_names <- scatter_names[1:2]

      }

    } else {

      scatter_list <-
        list(x_var = rep('loss_value', 2),
             y_var = c('accuracy', 'kappa'),
             x_lab = rep(paste0(ax_title, ', cross-validation'), 2),
             y_lab = c('overall accuracy, training',
                       "Cohen's \u03BA, training"))

      scatter_names <- c('loss_accuracy', 'loss_kapppa')

    }

    scatter_plots <-
      pmap(scatter_list,
           function(x_var, y_var, x_lab, y_lab) ggplot(data,
                                                       aes(x = .data[[x_var]],
                                                           y = .data[[y_var]])) +
             geom_point(shape = 21,
                        fill = point_color,
                        size = point_size,
                        alpha = point_alpha,
                        position = position_jitter(width = point_wjitter,
                                                   height = point_hjitter)) +
             geom_text_repel(aes(label = response),
                             size = txt_size,
                             color = txt_color) +
             theme_bw() +
             labs(title = 'Model performance',
                  x = x_lab,
                  y = y_lab))

    scatter_plots <- set_names(scatter_plots, scatter_names)

    c(violin_plots, scatter_plots)


  }

# Predictions --------

#' Predictions for a list `cv.glmnet` models.
#'
#' @description
#' Predicts responses, linear predictor scores, coefficients and similar,
#' for a list of `cv.glmnet` objects. Intended for the internal use.
#'
#' @return As described for \code{\link{predict.modTrain}}.
#'
#' @param models a list of `cv.glmnet` models.
#' @param lambdas a numeric vector of lambdas to make predictions for.
#' @param newx a numeric matrix of explanatory factors.
#' @param type of prediction, see: \code{\link{predict.modTrain}}.
#' @param family modeling family.
#' @param newoffset modeling offset.
#' @param ... extra arguments passed to \code{\link[glmnet]{predict.glmnet}}.

  pred_ <- function(models,
                    lambdas,
                    newx,
                    type,
                    family,
                    newoffset = NULL,
                    ...) {

    ## entry control is accomplished by upstream functions

    responses <- names(models)

    if(!is.null(newoffset)) {

      preds <- map2(models,
                    lambdas,
                    ~predict(.x,
                             newx = newx,
                             s = .y,
                             type = type,
                             newoffset = newoffset))

    } else {

      preds <- map2(models,
                    lambdas,
                    ~predict(.x,
                             newx = newx,
                             s = .y,
                             type = type,
                             newoffset = .x$glmnet.fit$offset))

    }

    if(type == 'class' | family %in% c('gaussian', 'poisson', 'binomial')) {

      for(i in seq_along(preds)) {

        colnames(preds[[i]]) <- responses[[i]]

      }

      pred_mat <- reduce(map(preds, as.matrix), cbind)

      ## for multinomial model class predictions, row names are not
      ## provided by default and will be appended

      if(family == 'multinomial') {

        rownames(pred_mat) <- rownames(newx)

      }

      return(pred_mat)

    }

    ## for predictions of multinomial models of type 'response' and 'class'
    ## lists of matrices with elements corresponding to the response variables
    ## are produced

    preds <- map(preds, ~.x[, , 1])

    return(preds)

  }

# END -----
