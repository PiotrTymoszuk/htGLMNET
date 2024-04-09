# S3 method for the `modData` class.

# Class testing -------

#' Test for an `modData` class instance.
#'
#' @description
#' Tests is an object is an instance of the `modData` class.
#'
#' @param x an object.
#'
#' @return a logical value.
#'
#' @export

  is_modData <- function(x) inherits(x, 'modData')

# Appearance -------

#' Appearance of the `modData` class.
#'
#' @description
#' `print()` method for the `modData` class.
#'
#' @param x an instance of the `modData` class.
#' @param ... extra arguments, currently none.
#'
#' @return none, called for its side effects.
#'
#' @export

  print.modData <- function(x, ...) {

    stopifnot(is_modData(x))

    cat(paste('modData object with', nrow(x$train), 'features'))

  }

# Numbers of observations and variables -----

#' Numbers of variables and observations.
#'
#' @description
#' Numbers of observations and variables in a `modData` object.
#'
#' @param object an instance of the `modData` class.
#' @param ... extra arguments, currently none.
#'
#' @return a numeric vector.
#'
#' @export nobs.modData
#' @export

  nobs.modData <- function(object, ...) {

    stopifnot(is_modData(object))

    n_features <- nrow(object$train)

    n_observations <- map(object[c('train', 'test')], ncol)

    c(n_training = n_observations[[1]],
      n_test = n_observations[[2]],
      n_variables = n_features)

  }

# Accessing the object's components ------

#' Access components of a `modData` object.
#'
#' @description
#' Accesses components of the `modData` object.
#'
#' @details
#' The following are returned for the argument `type` specified as:
#'
#' * `train`: the training numeric matrix with the selected modeling variables
#' and adjusted for the batch effects by ComBat.
#'
#' * `test`: the test numeric matrix with the selected modeling variables
#' and adjusted for the batch effects by ComBat.
#'
#' * `limits`: cutoffs of median, mean and Gini index used for selection of the
#' modeling variables.
#'
#' * `stats`: feature distribution statistics (median, mean, Gini index) prior
#' to the batch selection (all variables) and following the batch correction via
#' ComBat (modeling variables only).
#'
#' @return see: `Details`.
#'
#' @param object a `modData` object.
#' @param type the requested component of the object, see `Details`.
#' @param ... extra arguments, currently none.
#'
#' @export components.modData
#' @export

  components.modData <- function(object,
                                  type = c('train', 'test', 'limits', 'stats'),
                                  ...) {

    stopifnot(is_modData(object))

    type <- match.arg(type[1],
                      c('train', 'test', 'limits', 'stats'))

    if(type %in% c('train', 'test')) return(object[[type]])

    if(type == 'limits') {

      statistic <- NULL

      limits <-
        set_names(object[c('mean_limits', 'median_limits', 'gini_limits')],
                  c('mean', 'median', 'Gini index'))

      limits <- compact(limits)

      limits <- map2_dfr(limits, names(limits),
                         ~tibble(statistic = .y,
                                 train = .x[[1]],
                                 test = .x[[2]]))

      return(limits)

    } else {

      before_stats <- object[['stats']]

      after_stats <- map(object[c('train', 'test')],
                         row_stats)

      data_set <- NULL

      after_stats <-
        map2_dfr(after_stats, names(after_stats),
                 ~mutate(.x[c('variable', 'median', 'mean', 'gini_coef')],
                         data_set = .y))

      after_stats <- relocate(after_stats, data_set)

      return(list(before_adjustment = before_stats,
                  after_adjustment = after_stats))

    }

  }

# Summary --------

#' Summary of a `modData` object.
#'
#' @description
#' Computes distribution statistics (mean with SD, median with interquartile
#' range and 95% quantile range) for the entire training and test data sets
#' after the Combat adjustment.
#'
#' @return A data frame with whole-matrix distribution statistics.
#'
#' @param object a `modData` object.
#' @param ... extra arguments, currently none.
#'
#' @export summary.modData
#' @export

  summary.modData <- function(object, ...) {

    stopifnot(is_modData(object))

    map2_dfr(object[c('train', 'test')],
             c('train', 'test'),
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

# Plotting -------

#' Diagnostic plots for a `modData` object.
#'
#' @description
#' Generates a series of diagnostic plots for a `modData` object.
#'
#' @details
#' The following plots are returned:
#'
#' * violin plots of feature medians, means and Gini indexes prior to
#' the adjustment with the user-specified cutoffs denoted as dashed lines.
#'
#' * violin plots of feature medians and means after the batch effect
#' adjustment.
#'
#' @param x a `modData` object.
#' @param palette colors for the training and test data sets.
#' @param alpha alpha of the violins.
#' @param stat_color color of the points and whiskers representing
#' medians and interquartile ranges.
#' @param ... extra arguments passed to \code{\link[ggplot2]{geom_violin}}.
#'
#' @return a list of `ggplot` objects.
#'
#' @export plot.modData
#' @export

  plot.modData <- function(x,
                           palette = c(train = 'indianred',
                                       test = 'steelblue'),
                           alpha = 0.5,
                           stat_color = 'orangered2', ...) {

    stopifnot(is_modData(x))
    stopifnot(is.numeric(alpha))

    palette <- palette[1:2]

    ## plotting data -------

    plot_data <- components(x, 'stats')

    ft_numbers <- map(plot_data, ~table(.x[['data_set']]))

    train <- NULL
    test <- NULL

    axis_labs <- map(ft_numbers,
                     ~c(train = paste('train\nn =', .x['train']),
                        test = paste('test\nn =', .x['test'])))

    selection_limits <-
      set_names(x[c('mean_limits', 'median_limits', 'gini_limits')],
                c('mean', 'median', 'gini_coef'))

    ## medians and interquartile ranges to be shown in the plots -------

    median_tbl <- list()

    data_set <- NULL
    median_val <- NULL
    q25 <- NULL
    q75 <- NULL

    for(i in names(plot_data)) {

      group_data <- group_by(plot_data[[i]], data_set)

      median_tbl[[i]] <-
        map(c('mean', 'median', 'gini_coef'),
            ~summarise(group_data,
                       median_val = median(.data[[.x]], na.rm = TRUE),
                       q25 = quantile(.data[[.x]], 0.25, na.rm = TRUE),
                       q75 = quantile(.data[[.x]], 0.75, na.rm = TRUE)))


    }

    ## plotting --------

    diag_plots <- list()

    for(i in c('before_adjustment', 'after_adjustment')) {

      diag_plots[[i]] <-
        pmap(list(x = c('mean', 'median', 'gini_coef'),
                  y = c('Mean feature values',
                        'Median feature values',
                        'Gini coefficients'),
                  z = c('mean', 'median', 'Gini coefficient'),
                  w = median_tbl[[i]]),
             function(x, y, z, w) ggplot(plot_data[[i]],
                                         aes(x = data_set,
                                             y = .data[[x]],
                                             fill = data_set)) +
               geom_violin(alpha = alpha,
                           ...) +
               geom_errorbar(data = w,
                             aes(y = median_val,
                                 ymin = q25,
                                 ymax = q75),
                             width = 0,
                             color = stat_color) +
               geom_point(data = w,
                          aes(y = median_val),
                          shape = 18,
                          size = 3,
                          color = stat_color) +
               scale_x_discrete(limits = c('train', 'test'),
                                labels = axis_labs[[i]]) +
               scale_fill_manual(values = palette) +
               guides(fill = 'none') +
               theme(axis.title.x = element_blank()) +
               labs(title = y,
                    y = z))

      diag_plots[[i]] <- set_names(diag_plots[[i]],
                                   c('mean', 'median', 'gini_coef'))

    }

    ## appending the plots of distribution prior to the adjustment
    ## with the selection cutoffs

    for(i in names(diag_plots[[1]])) {

      if(is.null(selection_limits[[i]])) next

      diag_plots[[1]][[i]] <- diag_plots[[1]][[i]] +
        geom_hline(yintercept = selection_limits[[i]][[1]],
                   color = palette[[1]],
                   linetype = 'dashed') +
        geom_hline(yintercept = selection_limits[[i]][[2]],
                   color = palette[[2]],
                   linetype = 'dashed')

    }

    map(diag_plots,
        map, ~.x + theme_bw())

  }

# END ------
