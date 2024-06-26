% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/high_throughput.R
\name{fitHT}
\alias{fitHT}
\alias{fitHT.default}
\alias{fitHT.modData}
\title{High-throughput regularized linear modeling.}
\usage{
fitHT(...)

\method{fitHT}{default}(
  x_train,
  x_test,
  y,
  trans_fun = identity,
  standardize = FALSE,
  alpha = 1,
  family = "gaussian",
  type.measure = c("default", "mse", "deviance", "class", "auc", "mae", "C"),
  type.predict = c("link", "response", "class"),
  safe = FALSE,
  preprocess.args = list(),
  predict.args = list(),
  ...
)

\method{fitHT}{modData}(
  x,
  y,
  trans_fun = identity,
  standardize = FALSE,
  alpha = 1,
  family = "gaussian",
  type.measure = c("default", "mse", "deviance", "class", "auc", "mae", "C"),
  type.predict = c("link", "response", "class"),
  safe = FALSE,
  predict.args = list(),
  ...
)
}
\arguments{
\item{...}{additional arguments passed to \code{\link[glmnet]{cv.glmnet}} and
\code{\link[glmnet]{glmnet}}.}

\item{x_train}{a numeric matrix with training explanatory variables.
Features are coded by rows, observations are provided in columns.}

\item{x_test}{a numeric matrix with test explanatory variables.
Features are coded by rows, observations are provided in columns.}

\item{y}{a numeric or factor matrix with the observed responses
(e.g. phenotype in the training data). The responses are provided as columns,
the observations are given in rows. The row names of \code{y} must match column
name of \code{x}.}

\item{trans_fun}{a numeric function to transform responses for modeling.
A normality-improving transformation such as \code{log} may be provided here.
Used only for Gaussian family models and ignored elsewhere.}

\item{standardize}{logical, should \code{x} values be standardized prior
to modeling?}

\item{alpha}{numeric mixing parameter, which defines the modeling algorithm
(RIDGE regression: \code{alpha = 0}, LASSO: \code{alpha = 1}, Elastic Net otherwise).
See \code{\link[glmnet]{glmnet}} for details.}

\item{family}{family of the \code{glmnet} models. Currently, only 'gaussian',
'binomial', 'poisson', and 'multinomial' are implemented.
See: \code{\link[glmnet]{glmnet}} for details.}

\item{type.measure}{model selection statistic (i.e. loss function). Defaults
to deviance. See \code{\link[glmnet]{cv.glmnet}} for details.}

\item{type.predict}{prediction type. See: \code{\link{predict.modTrain}} for
details.}

\item{safe}{logical, should modeling proceed in spite of errors for
some responses? .
If \code{TRUE}, any errors are logged in \code{errors} element of the output. This
option may be, however, considerably slower.}

\item{preprocess.args}{a named list with arguments passed to
\code{\link{pre_process}}, such as variable selection cutoffs. Not available
if a \code{modData} object is used as \code{x}.}

\item{predict.args}{a named list with arguments passed to
\code{\link{predict.modTrain}}, such as variable offset.}

\item{x}{an object of the \code{modData} class created e.g. with the
\code{\link{pre_process}} function.}
}
\value{
an instance of the \code{modHT} class with the following components:
\itemize{
\item \code{predictions}: model predictions for the test data as a matrix or a list
of matrices.
\item \code{preprocess_stats}: information of modeling variables after
the pre-processing step, like variable name, median and mean value,
Gini index.
\item \code{model_stats}: metrics of the regularized linear models, such as loss
function used for \eqn{\lambda} tuning, \eqn{R^2}, correlations of the
observed and predicted modeling response in the training data (regression
models), as well as overall accuracy and Cohen's \eqn{\kappa} for comparison
of the observed and predicted class assignment (classification models).
\item \code{preprocess_globals}: pre-processing parameters specified by the user like
modeling variable selection cutoffs or the function used to transform the
explanatory variables.
\item \code{model_globals}: modeling parameters specified by the user such as modeling
family, the function used to transform the response variables, loss function,
mixing parameter alpha.
\item \code{prediction_globals}: parameters of predictions such as prediction type.
\item \code{errors}: an optional vector of error messages.
}

Note the methods \code{\link{nobs.modHT}}, \code{\link{components.modHT}},
\code{\link{summary.modHT}}, and \code{\link{plot.modHT}} defined
for the class.
}
\description{
Functions for light-weight, high-throughput modeling of multiple responses
with common large set of explanatory data.
}
\details{
The \code{fitHT()} functions bundle pre-processing (modeling
variable selection, batch-adjustment), tuning and training (choice of
\eqn{\lambda}), model evaluation, and predictions.
Technically, the pre-processing step is accomplished by the
\code{\link{pre_process}} function and is skipped when the user provides
already pre-processed data as an \code{\link{modData}} object.
For the tuning, training and evaluation steps, the function
\code{\link{train.default}} is employed.
Finally, predictions in a format requested by the user and specified by the
\code{type} argument, are generated by the \code{\link{predict.modTrain}} method.
The tools will run in parallel, is a furrr-compatible backend is declared,
e.g. with \code{\link[future]{plan}}.
}
\references{
Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD. The sva package for
removing batch effects and other unwanted variation in high-throughput
experiments. Bioinformatics (2012) 28:882. doi:10.1093/BIOINFORMATICS/BTS034

Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized
linear models via coordinate descent.
J Stat Softw (2010) 33:1–22. doi:10.18637/jss.v033.i01
}
