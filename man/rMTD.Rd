% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simu.R
\name{rMTD}
\alias{rMTD}
\title{Simulation from an MTD model}
\usage{
rMTD(mtdorder, weight, family, param, size, trun = FALSE, init_data = NULL)
}
\arguments{
\item{mtdorder}{the model order. Note mtdorder > 1.}

\item{weight}{a numeric vector of weights of the MTD model.}

\item{family}{a quoted keyword that specifies the stationary marginal distribution.
Supported keywords are: "gaussian", "poisson, negative binomial and lomax".}

\item{param}{a list of pre-specified parameter values of the model.}

\item{size}{the number of data to be generated conditional on the initial values.}

\item{trun}{a logical value. If TRUE, the first mtdorder elements need to be specified.}

\item{init_data}{a numeric vector of the initial values of the data to be generated. 
The length of the vector is mtdorder.}
}
\value{
The return object is a list comprising:
\describe{
  \item{y}{a numeric vector of simulated data.}
  \item{label}{a numeric vector of configuration variables that identify which component of the mtd was used to generate the data.}
  \item{lags}{a numeric vector that contains the specific lags were used to generate the data.}
  \item{u}{if there is, a numeric vector of the auxiliary variables used to generate the data.}
}
}
\description{
This function generates data from an MTD model with a given family of stationary marginal distribution.
}
\references{
Zheng, X., Kottas, A., and Sansó, B. (2021),
"On construction and estimation of stationary mixture transition distribution models,"
\emph{arXiv:2010.12696}.
}
