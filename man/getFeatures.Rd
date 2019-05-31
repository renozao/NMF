\name{getFeatures}
\alias{getFeatures}
\title{It retrieves the top most contributing feature associated with a metagene estimated by \code{\link{nmf}}}
\description{It applies the \code{\link{extractFeatures}} function with the method of choice (default is 'max')
to identify the top most contributing features of the metagenes.}
\usage{
getFeatures(nmf.result, method = "max", original.data)
}

\arguments{
  \item{nmf.result}{Object returned by \code{\link{nmf}} function, run with a single rank and a single algorithm.}
  \item{method}{Scoring or selection method ("kim" or "max").}
  \item{original.data}{The original data matrix submiited to the \code{\link{nmf}} analysis.}
}
\value{list with the selected features.}

\author{
Cintia C Palu
}
