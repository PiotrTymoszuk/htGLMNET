# Tools for prediction of drug sensitivity

# packages ------

  library(tidyverse)
  library(trafo)
  library(rlang)

  library(sva)
  library(caret)
  library(caretExtra)
  library(glmnet)

  library(microViz)

# model training and predictions -------

  train_drugs <- function(train_expression,
                          train_phenotype,
                          family = 'gaussian',
                          alpha = 0.5, ...) {

    ## trains a series of GLMNET models of the drug's IC50
    ## lambda is found by CV tuning. Returns a list of CV stats, best lambdas
    ## and the optimized models

    ## tuning --------

    x_mat <- t(train_expression)

    y_mats <- map(colnames(train_phenotype),
                  ~train_phenotype[, .x])

    y_mats <- set_names(y_mats, colnames(train_phenotype))

    y_mats <- map(y_mats, ~.x[!is.na(.x)])

    cv_models <-
      future_map(y_mats,
                 ~cv.glmnet(x = x_mat[names(.x), ],
                            y = log(.x),
                            family = family,
                            alpha = alpha, ...),
                 .options = furrr_options(seed = TRUE))

    ## cross-validated stats and best tunes -------

    tunes <- map(cv_models, ~as_tibble(.x[c('lambda', 'cvm')]))

    best_tunes <-
      map(tunes,
          filter,
          cvm == min(cvm, na.rm = TRUE))

    best_tunes <- map(best_tunes, ~.x[1, ])

    cv_stats <-
      tibble(variable = names(y_mats),
             lambda = map_dbl(best_tunes, ~.x[['lambda']]),
             variance = map_dbl(y_mats, ~var(log(.x))),
             mse = map_dbl(best_tunes, ~.x[['cvm']]))

    cv_stats <- mutate(cv_stats,
                       rsq = 1 - mse/variance)

    ## glmet models ----------

    train_models <-
      map2(y_mats, cv_stats$lambda,
           ~glmnet(x = x_mat[names(.x), ],
                   y = log(.x),
                   lambda = .y,
                   family = family,
                   alpha = alpha, ...))

    list(cv_stats = cv_stats,
         glmnet_models = train_models)

  }

  calculate_drug_scores <- function(glmnet_models,
                                    expression) {
    scores <- map(glmnet_models,
                  ~predict.glmnet(.x,
                                  newx = t(expression),
                                  newoffset = .x$offset))

    scores <- map2(scores, names(scores),
                   set_colnames)

    scores <- as.data.frame(reduce(scores, cbind))

    as_tibble(rownames_to_column(scores, 'sample_id'))

  }

# A common function for calculating the drug sensitivity scores -----
