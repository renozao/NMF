\name{setNMFClass}
\alias{setNMFClass}
\title{Factory Method for NMF Model Classes}
\usage{
  setNMFClass(Class, ..., where = toppackage(),
    contains = "NMFstd", VERBOSE = TRUE)
}
\description{
  Defines two S4 classes for representing NMF models: one
  to hold data from the actual model, the other one to hold
  fitting data for model estimated with the function
  \code{\link{nmf}}.
}

