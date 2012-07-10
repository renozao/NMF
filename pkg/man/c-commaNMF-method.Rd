\docType{methods}
\name{c,NMF-method}
\alias{c,NMF-method}
\title{Concatenating NMF Models}
\usage{
  \S4method{c}{NMF}(x, ..., margin = 3, recursive = FALSE)
}
\arguments{
  \item{x}{an NMF model}

  \item{...}{other objects to concatenate. Currently only
  two objects at a time can be concatenated (i.e. \code{x}
  and \code{..1}).}

  \item{margin}{integer that indicates the margin along
  which to concatenate (only used when \code{..1} is a
  matrix): \describe{ \item{1L}{} \item{2L}{} \item{3L}{}
  \item{4L}{} } If missing the margin is heuristically
  determined by looking at common dimensions between the
  objects.}

  \item{recursive}{logical.  If \code{recursive = TRUE},
  the function recursively descends through lists (and
  pairlists) combining all their elements into a vector.}
}
\description{
  Binds compatible matrices and NMF models together.
}
\keyword{internal}
\keyword{methods}

