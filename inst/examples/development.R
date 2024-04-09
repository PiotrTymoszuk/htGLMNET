# Development and testing of the package components

# tools and data -------

  load('./inst/examples/pck_devel.RData')

  library(htGLMNET)
  library(glmnet)

  library(tidyverse)
  library(microViz)
  library(trafo)
  library(furrr)

# Explanatory variables ----------

  pre_pro_cptac <-
    pre_process(pck_devel$train_expression,
                pck_devel$test_expression$cptac,
                trans_fun = function(x) x,
                mean_quantile = NULL,
                median_quantile = 0.2,
                gini_quantile = 0.2)

  pre_pro_tcga <-
    pre_process(pck_devel$train_expression,
                pck_devel$test_expression$tcga,
                trans_fun = function(x) x,
                mean_quantile = NULL,
                median_quantile = 0.2,
                gini_quantile = 0.2)

# Responses --------

  ## numeric responses, Gaussian and Poisson family

  gaussian_responses <- pck_devel$train_pheno[, 1:50] * 10^6

  poisson_responses <- as.data.frame(pck_devel$train_pheno[, 1:50])

  poisson_responses <- poisson_responses %>%
    map_dfc(rank) %>%
    as.matrix %>%
    set_rownames(rownames(poisson_responses))

  ## binomial responses: IC50 below median is deemed sensitive

  binomial_responses <- as.data.frame(pck_devel$train_pheno[, 1:10])

  binomial_cuts <- colMedians(binomial_responses)

  binomial_responses <-
    map2_dfc(binomial_responses, binomial_cuts,
         ~cut(.x, c(-Inf, .y, Inf))) %>%
    map_dfc(~as.numeric(.x) - 1) %>%
    as.matrix %>%
    set_rownames(rownames(binomial_responses))

  ## multinomial responses

  multi_responses <- as.data.frame(pck_devel$train_pheno[, 1:10])

  multi_responses <- multi_responses %>%
    map_dfc(~cut(.x,
                 c(-Inf, quantile(.x, c(0.25, 0.5, 0.75), na.rm = TRUE), Inf))) %>%
    map_dfc(~as.numeric(.x) - 1) %>%
    as.matrix %>%
    set_rownames(rownames(multi_responses))

# Tuning and training ---------

  plan('multisession')

  train_gaussian <- train(x = pre_pro_cptac,
                          y = gaussian_responses,
                          trans_fun = log,
                          alpha = 0.5,
                          standardize = FALSE,
                          safe = TRUE)

  train_poisson <- train(x = pre_pro_cptac,
                         y = poisson_responses,
                         family = 'poisson',
                         alpha = 0.5,
                         standardize = FALSE,
                         safe = TRUE)

  train_binomial <- train(x = pre_pro_cptac,
                          y = binomial_responses,
                          family = 'binomial',
                          alpha = 0.5,
                          standardize = FALSE,
                          safe = TRUE)

  train_multinomial <- train(x = pre_pro_cptac,
                             y = multi_responses,
                             family = 'multinomial',
                             alpha = 0.5,
                             standardize = FALSE,
                             safe = TRUE)

  plan('sequential')

# Predictions --------

  test_gaussian <- predict(train_gaussian,
                           newdata = pre_pro_cptac,
                           type = 'link')

  test_poisson <- predict(train_poisson,
                          newdata = pre_pro_cptac,
                          type = 'response')

  test_binomial <- predict(train_binomial,
                           newdata = pre_pro_cptac,
                           type = 'response')

  test_multinomial <- predict(train_multinomial,
                              newdata = pre_pro_cptac,
                              type = 'class')

# High throughput modeling ------

  plan('multisession')

  ht_gaussian <- fitHT.default(x_train = pck_devel$train_expression,
                               x_test = pck_devel$test_expression$tcga,
                               y = gaussian_responses,
                               trans_fun = log,
                               alpha = 0.5,
                               standardize = FALSE,
                               safe = TRUE)

  plan('sequential')

  plot(ht_gaussian)

  plot(train_gaussian)

# END ------
