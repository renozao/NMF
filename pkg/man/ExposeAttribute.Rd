\name{ExposeAttribute}
\alias{EAmode}
\alias{EAmode<-}
\alias{ExposeAttribute}
\title{Exposing Object Attributes}
\usage{
  ExposeAttribute(object, mode = "rw")

  EAmode(x)

  EAmode(x) <- value
}
\arguments{
  \item{object}{any R object whose attributes need to be
  exposed}

  \item{mode}{access mode: \describe{ \item{\dQuote{r}:}{
  (read-only) only method \code{$} is defined}
  \item{\dQuote{w}:}{ (write-only) only method \code{$<-}
  is defined} \item{\dQuote{rw}:}{ (read-write) both
  methods \code{$} and \code{$<-} are defined} }}

  \item{x}{an \code{ExposeAttribute} object}

  \item{value}{replacement value for mode.}
}
\description{
  The function \code{ExposeAttribute} creates an S3 object
  that exposes all attributes of any R object, by making
  them accessible via methods \code{\link{$}} and/or
  \code{\link{$<-}}.

  \code{EAmode} and \code{EAmode<-} get and sets the access
  mode of \code{ExposeAttribute} objects.
}

