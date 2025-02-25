% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\name{modHT}
\alias{modHT}
\title{\code{modHT} class.}
\usage{
modHT(x, ...)
}
\arguments{
\item{x}{a named list with elements described in \code{Details}.}

\item{...}{extra arguments, currently none.}
}
\value{
an instance of the \code{modHT} class with components
described in \code{Details}.
}
\description{
Creates an object of the \code{modHT} class, which stores model predictions in
the training data set along with pre-processing and modeling performance
statistics. It can be created e.g. with the \code{\link{fitHT}} function.
}
\details{
The object has the following components:
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
