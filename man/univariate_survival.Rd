% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predsurv.R
\name{univariate_survival}
\alias{univariate_survival}
\title{Fit univariate survival
Get the cost and coefficients for univariate survival analysis.}
\usage{
univariate_survival(d_train, x_vars, d_test = NULL, weights = NULL,
  ncores = 1L, time_var = "time", status_var = "status",
  cost = "ll", mod = c("cox", "reg"), save_file)
}
\arguments{
\item{x_vars}{vars to analysis}

\item{ncores}{for parallel.}

\item{time_var}{defults to time}

\item{status_var}{defaults to status}

\item{cost}{cost measure}

\item{save_file}{should results be saved on file}

\item{dat}{dataframe}
}
\description{
Fit univariate survival
Get the cost and coefficients for univariate survival analysis.
}
