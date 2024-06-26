% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modHT_oop.R
\name{components.modHT}
\alias{components.modHT}
\title{Access components of a \code{modHT} object.}
\usage{
\method{components}{modHT}(
  object,
  type = c("predictions", "preprocess_stats", "model_stats", "globals", "errors"),
  ...
)
}
\arguments{
\item{object}{a \code{modHT} object.}

\item{type}{the requested component of the object, see \code{Details}.}

\item{...}{extra arguments, currently none.}
}
\value{
see: \code{Details}.
}
\description{
Accesses components of the \code{modHT} object.
}
\details{
The following are returned for the argument \code{type} specified as:
\itemize{
\item \code{predictions}: the prediction matrix or a list of prediction matrices.
\item \code{preprocess_stats}: feature distribution statistics after batch correction.
\item \code{model_stats}: model performance statistics.
\item \code{globals}: global parameter values used for pre-processing, modeling
and predictions.
\item \code{errors}: error messages for failed modeling tasks.
}
}
