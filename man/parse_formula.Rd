\name{parse_formula}
\alias{parse_formula}
\title{Simple Parsing of Formula}
\usage{
  parse_formula(x)
}
\arguments{
  \item{x}{formula to parse}
}
\value{
  a list with the following elements: \item{response}{
  logical that indicates if the formula has a response
  term.} \item{y}{ name of the response variable.}
  \item{x}{ list of regressor variable names.} \item{n}{
  number of regressor variables.}
}
\description{
  Formula parser for formula-based NMF models.
}
\keyword{internal}

