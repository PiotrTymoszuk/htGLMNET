# Example code for the GitHub repository of the package

# tools and data -------

  library(tidyverse)
  library(htGLMNET)
  library(rstatix)

  library(future)

  library(patchwork)

  ## training data: GDSC experiment:
  ## whole-genome expression in untreated cancer cell lines,
  ## log transformed to improve normality

  training_x <- gdsc_expression

  training_x[1:5, 1:5]

  ## test data: whole-tumor (bulk) gene expression in pancreatic cancers
  ## of the CPTAC and TCGA cohorts
  ## already log-transformed, to improve normality

  test_x <- list(CPTAC = paad_cptac,
                 TCGA = paad_tcga)

  test_x$CPTAC[1:5, 1:5]

  ## IC50 concentrations in µM of ten anti-cancer drugs
  ## log transformation greatly improves normality!
  ## Yet, `glmnet` may have problems with negative values after
  ## the log transformation, to avoid that, we're transforming the IC50 to pM

  training_ic50 <- gdsc_ic50

  training_ic50[1:5, 1:5]

  training_ic50 %>%
    as.data.frame %>%
    map_dfr(shapiro_test) %>%
    mutate(variable = colnames(training_ic50))

  training_ic50 %>%
    log %>%
    as.data.frame %>%
    map_dfr(shapiro_test) %>%
    mutate(variable = colnames(training_ic50))

  ## picomole conversion and log transformation

  training_ic50 <- log(1e6 * training_ic50)

# Pre-processing ---------

  ## pre-processing involves the following steps:
  ##
  ## 1) selection of common genes in the GDSC cell lines and pancreatic cancers
  ## such genes must be fairly abundantly expressed and variable,
  ## hence using 0.2 quantiles of median expression and Gini coefficients
  ## as modeling variable selection criteria
  ##
  ## 2) batch adjustment via ComBat

  pre_pro_data <- multi_process(train = training_x,
                                test = test_x,
                                median_quantile = 0.2,
                                gini_quantile = 0.2)

  ## the effects of batch adjustment and common gene numbers
  ## can be inspected by `summary()` and `plot()`

  summary(pre_pro_data)

  pre_pro_plots <- plot(pre_pro_data)

  ## easy customization of the diagnostic plots

  pre_pro_plots$before_adjustment$median +
    labs(title = 'Median, before ComBat',
         x = 'Data set') +
    scale_y_continuous(limits = c(0, 15)) +
    pre_pro_plots$after_adjustment$median +
    labs(title = 'Median, after ComBat',
         x = 'Data set') +
    scale_y_continuous(limits = c(0, 15))

  pre_pro_plots$before_adjustment$gini_coef +
    labs(title = 'Gini index, before ComBat',
         x = 'Data set') +
    scale_y_continuous(limits = c(0, 1)) +
    pre_pro_plots$after_adjustment$gini_coef +
    labs(title = 'Gini index, after ComBat',
         x = 'Data set') +
    scale_y_continuous(limits = c(0, 1))

# Tuning, training and evaluation ---------

  ## Gaussian family Elastic Net models are fitted
  ## (family = 'gaussian', alpha = 0.5)
  ##
  ## the optimal lambda is found by minimizing mean squared error
  ## in cross-validation (type.measure = 'mse')
  ##
  ## no standardization, because the gene expression values are expected
  ## to be on the same scale.
  ##
  ## The function runs in parallel upon defining a multisession backend

  plan('multisession')

  set.seed(12345)

  train_object <- train(x = pre_pro_data,
                        y = training_ic50,
                        standardize = FALSE,
                        family = 'gaussian',
                        alpha = 0.5,
                        type.measure = 'mse')

  plan('sequential')

  ## model evaluation stats can be explored by `summary()` and `plot()`

  summary(train_object)

  ## again, diagnostic plots in `ggplot` format can be easily
  ## customized. In this cases, we're adding lines of cutoffs of
  ## substantial correlation (r >= 0.5) and explained variance (R^2 >= 0.26)

  train_plots <- plot(train_object)

  train_plots$rsq_pearson +
    labs(title = 'Molde performance: CV and training') +
    geom_hline(yintercept = 0.5,
               linetype = 'dashed') +
    geom_vline(xintercept = 0.26,
               linetype = 'dashed')

# Gene signatures of drug response -------

  ## gene signatures of drug response may be obtained by extracting
  ## non-zero coefficients of the regularized model.
  ## By this way, we're selecting only the features that contribute to
  ## prediction of drug sensitivity
  ##
  ## Coefficients of the models are obtained via `coef()`, which returns
  ## a matrix for all investigated genes.

  train_coefficients <- coef(train_object)

  train_coefficients[1:5, 1:5]

  ## removal of the intercepts and selection of non-zero model coefficients

  non_zero <- train_coefficients[-1, ]

  non_zero <- colnames(non_zero) %>%
    map(~non_zero[, .x]) %>%
    map(as.data.frame) %>%
    map(set_names, 'beta')

  non_zero <- non_zero %>%
    map(filter, beta != 0) %>%
    map(rownames_to_column, 'gene_symbol') %>%
    map(arrange, -beta) %>%
    map(as_tibble) %>%
    set_names(colnames(train_coefficients))

  ## plotting coefficients of the signature members for cis-platin

  non_zero$Cisplatin_1005 %>%
    ggplot(aes(x = beta,
               y = reorder(gene_symbol, beta),
               fill = factor(sign(beta)))) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = c('-1' = 'steelblue',
                                 '1' = 'indianred3')) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(face = 'italic'),
          legend.position = 'none') +
    labs(title = 'Cis-platin response signature',
         x = expression(beta * ', association with IC50'))

# Prediction of drug sensitivity for the bulk cancer data -------

  ## IC50 samples for individuals cancers may be computed
  ## by calling `predict(type = 'response')`.
  ## Of convenience, as `newdata`, the result of `pre_process()` may be
  ## supplied. The function picks the test data set automatically!

  test_ic50 <- predict(train_object,
                       newdata = pre_pro_data$test,
                       type = 'response')

  ## a list of matrices is returned
  ## transforming the predicted sensitivities back to µM

  test_ic50 <- test_ic50 %>%
    map(~exp(.x) * 1e-6)

  test_ic50$CPTAC[1:5, 1:5]

  ## let's visualize predictions for
  ## the models with the best performance
  ## (cross-validated R^2 >= 0.26, Pearson's r >= 0.5)

  best_models <- train_object %>%
    summary %>%
    filter(rsq_oof >= 0.26,
           pearson > 0.5) %>%
    .$response

  ## plotting data in long format

  test_plot_data <- test_ic50 %>%
    map(~.x[, best_models]) %>%
    map(as.data.frame) %>%
    map(as_tibble)

  test_plot_data <-
    map2_dfr(test_plot_data, names(test_plot_data),
             ~mutate(.x, cohort = .y))

  test_plot_data <- test_plot_data %>%
    pivot_longer(cols = all_of(best_models),
                 names_to = 'compound',
                 values_to = 'ic50')

  test_plot_data %>%
    ggplot(aes(x = reorder(compound, ic50),
               y = ic50,
               fill = cohort)) +
    geom_violin(scale = 'width',
                position = position_dodge(0.9)) +
    geom_point(shape = 16,
               size = 1,
               position = position_jitterdodge(jitter.width = 0.05,
                                               jitter.height = 0,
                                               dodge.width = 0.9),
               alpha = 0.25,
               show.legend = FALSE) +
    scale_y_continuous(trans = 'log10') +
    guides(x = guide_axis(angle = 45)) +
    theme_bw() +
    labs(title = 'Predicted IC50',
         subtitle = 'Bulk cancer samples',
         x = 'Compound',
         y = 'IC50, µM')

# END -----
