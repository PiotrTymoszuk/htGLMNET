# Preparing example data sets to be delivered with the package.
#
# The training data sets for the GDSC2 drug screening were
# obtained from https://osf.io/c6tfx/:
#
# a matrix of log-transformed gene expression values (RMA normalized) for
# untreated cell lines.
#
# a matrix of IC50 concentrations of anti-cancer drugs in the cell lines.
#
# The test data are log-transformed gene expression counts  for pancreatic
# cancer samples of he CPTAC cohort. The training data were fetched from
# https://www.cbioportal.org/study/summary?id=paad_cptac_2021

# processing and saving --------

  load('./inst/examples/pck_devel.RData')

  ## test set expression: CPTAC

  paad_cptac <- pck_devel$test_expression$cptac

  save(paad_cptac, file = './data/paad_cptac.RData', compress = 'xz')

  ## test set expression: TCGA

  paad_tcga <- pck_devel$test_expression$tcga

  save(paad_tcga, file = './data/paad_tcga.RData', compress = 'xz')

  ## training data set phenotype: the first 10 drugs
  ## and randomly selected 200 cell lines with the complete data

  gdsc_ic50 <- pck_devel$train_pheno[, 1:10]

  gdsc_ic50 <- gdsc_ic50[complete.cases(gdsc_ic50), ]

  set.seed(12345)

  gdsc_cell_lines <- sample(rownames(gdsc_ic50), size = 200, replace = FALSE)

  gdsc_ic50 <- gdsc_ic50[gdsc_cell_lines, ]

  save(gdsc_ic50, file = './data/gdsc_ic50.RData', compress = 'xz')

  gdsc_expression <- pck_devel$train_expression

  gdsc_expression <- gdsc_expression[, gdsc_cell_lines]

  save(gdsc_expression, file = './data/gdsc_expression.RData', compress = 'xz')

# END -----
