# Descriptions of the data sets:
#
# The general idea of prediction of drug sensitivity with regularized regression
# comes from two papers:
#
# Geeleher P, Cox NJ, Huang RS. Clinical drug response can be predicted
# using baseline gene expression levels and in vitro drug sensitivity in
# cell lines. Genome Biol (2014) 15: doi:10.1186/gb-2014-15-3-r47
#
# Geeleher P, Cox N, Stephanie Huang R. pRRophetic: An R Package for Prediction
# of Clinical Chemotherapeutic Response from Tumor Gene Expression Levels.
# PLoS One (2014) 9:e107468. doi:10.1371/JOURNAL.PONE.0107468
#
# Maeser D, Gruener RF, Huang RS. oncoPredict: an R package for predicting
# in vivo or cancer patient drug response and biomarkers from cell line
# screening data. Brief Bioinform (2021) 22:1–7. doi:10.1093/BIB/BBAB260
#
# The training data sets for the GDSC2 drug screening:
#
# Yang W, Soares J, Greninger P, Edelman EJ, Lightfoot H, Forbes S,
# Bindal N, Beare D, Smith JA, Thompson IR, et al. Genomics of Drug
# Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker
# discovery in cancer cells. Nucleic Acids Res (2013) 41:D955–D961.
# doi:10.1093/NAR/GKS1111
#
# were obtained from https://osf.io/c6tfx/, the repository by Measer et al.:
#
# a matrix of log-transformed gene expression values (RMA normalized) for
# untreated cell lines.
#
# a matrix of IC50 concentrations of anti-cancer drugs in the cell lines.
#
# The test data are log-transformed gene expression counts  for pancreatic
# cancer samples of the CPTAC cohort:
#
# Cao L, Huang C, Cui Zhou D, Hu Y, Lih TM, Savage SR, Krug K, Clark DJ,
# Schnaubelt M, Chen L, et al. Proteogenomic characterization of pancreatic
# ductal adenocarcinoma. Cell (2021) 184:5031-5052.e26.
# doi:10.1016/J.CELL.2021.08.023
#
# The training data were fetched from:
# https://www.cbioportal.org/study/summary?id=paad_cptac_2021

#' Baseline gene expression of cancer cell lines, GDSC2 study.
#'
#' @description
#' RMA-normalized and log-transformed gene expression values in 200 randomly
#' selected cancer cell lines from the GDSC2 drug treatment screening without
#' any treatment. The data set may be used as an example training set for
#' exploring tools of the package. The full data set is available
#' in the repository by Maeser and colleagues at https://osf.io/c6tfx/.
#'
#' @references
#' Yang W, Soares J, Greninger P, Edelman EJ, Lightfoot H, Forbes S,
#' Bindal N, Beare D, Smith JA, Thompson IR, et al. Genomics of Drug
#' Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker
#' discovery in cancer cells. Nucleic Acids Res (2013) 41:D955–D961.
#' doi:10.1093/NAR/GKS1111
#'
#' @references
#' Maeser D, Gruener RF, Huang RS. oncoPredict: an R package for predicting
#' in vivo or cancer patient drug response and biomarkers from cell line
#' screening data. Brief Bioinform (2021) 22:1–7. doi:10.1093/BIB/BBAB260
#'
#' @references
#' Geeleher P, Cox NJ, Huang RS. Clinical drug response can be predicted
#' using baseline gene expression levels and in vitro drug sensitivity in
#' cell lines. Genome Biol (2014) 15: doi:10.1186/gb-2014-15-3-r47
#'
#' @references
#' Geeleher P, Cox N, Stephanie Huang R. pRRophetic: An R Package for Prediction
#' of Clinical Chemotherapeutic Response from Tumor Gene Expression Levels.
#' PLoS One (2014) 9:e107468. doi:10.1371/JOURNAL.PONE.0107468
#'
#' @format
#' a numeric matrix of RMA normalized and log-transformed gene
#' expression levels, genes in rows and cell line identifiers in columns.
#' The matrix has 17419 rows and 200 columns.
#'
#' @source https://osf.io/c6tfx/
#'
#' @docType data
#'
#' @name gdsc_expression
#'
#' @usage data(gdsc_expression)

  NULL

#' Drug sensitivity. GDSC2 study.
#'
#' @description
#' Untransformed IC50 concentrations for 200 randomly selected cell lines
#' treated with 10 cancer drugs. The data set may be used for exploring tools
#' of the package. Its full version is provided in a repository by Maeser at al.
#' (https://osf.io/c6tfx/).
#'
#' @references
#' Yang W, Soares J, Greninger P, Edelman EJ, Lightfoot H, Forbes S,
#' Bindal N, Beare D, Smith JA, Thompson IR, et al. Genomics of Drug
#' Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker
#' discovery in cancer cells. Nucleic Acids Res (2013) 41:D955–D961.
#' doi:10.1093/NAR/GKS1111
#'
#' @references
#' Maeser D, Gruener RF, Huang RS. oncoPredict: an R package for predicting
#' in vivo or cancer patient drug response and biomarkers from cell line
#' screening data. Brief Bioinform (2021) 22:1–7. doi:10.1093/BIB/BBAB260
#'
#' @references
#' Geeleher P, Cox NJ, Huang RS. Clinical drug response can be predicted
#' using baseline gene expression levels and in vitro drug sensitivity in
#' cell lines. Genome Biol (2014) 15: doi:10.1186/gb-2014-15-3-r47
#'
#' @references
#' Geeleher P, Cox N, Stephanie Huang R. pRRophetic: An R Package for Prediction
#' of Clinical Chemotherapeutic Response from Tumor Gene Expression Levels.
#' PLoS One (2014) 9:e107468. doi:10.1371/JOURNAL.PONE.0107468
#'
#' @format
#' a numeric matrix of non-transformed IC50 concentrations (µM) for 200 randomly
#' selected cell lines. Cell line identifiers are provided in rows, drugs are
#' placed in columns. The matrix has 200 rows and 10 columns.
#'
#' @source https://osf.io/c6tfx/
#'
#' @docType data
#'
#' @name gdsc_ic50
#'
#' @usage data(gdsc_ic50)

  NULL

#' Bulk cancer RNA sequencing, PAAD CPTAC study.
#'
#' @description
#' Bulk cancer expression data (RSEM, log transformed) for 140 pancreatic
#' ductal cancer samples from the PAAD CPTAC study.
#'
#' @references
#' Cao L, Huang C, Cui Zhou D, Hu Y, Lih TM, Savage SR, Krug K, Clark DJ,
#' Schnaubelt M, Chen L, et al. Proteogenomic characterization of pancreatic
#' ductal adenocarcinoma. Cell (2021) 184:5031-5052.e26.
#' doi:10.1016/J.CELL.2021.08.023
#'
#' @format
#' a numeric matrix with RNA sequencing gene expression levels normalized by
#' RSEM and log-transformed. The matrix contains 2364 rows and 140 columns.
#' Genes (symbols) are placed in rows, cancer sample identifiers are in columns.
#'
#' @source https://www.cbioportal.org/study/summary?id=paad_cptac_2021
#'
#' @docType data
#'
#' @name paad_cptac
#'
#' @usage data(paad_cptac )

  NULL

# END ------
