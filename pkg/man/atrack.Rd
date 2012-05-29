\docType{methods}
\name{.atrack}
\alias{adata}
\alias{alength}
\alias{amargin}
\alias{anames}
\alias{annotationTrack}
\alias{atrack}
\alias{.atrack}
\alias{.atrack,ANY-method}
\alias{.atrack,data.frame-method}
\alias{.atrack,matrix-method}
\alias{.atrack-methods}
\alias{is.atrack}
\title{Annotation Tracks}
\usage{
  .atrack(object, ...)

  is.atrack(x)

  adata(x, value, ...)

  amargin(x, value)

  anames(x, default.margin)

  alength(x, default.margin)

  \S4method{.atrack}{ANY}(object, data = NULL, ...)

  atrack(..., order = NULL, enforceNames = FALSE,
    .DATA = NULL)

  annotationTrack(x)
}
\arguments{
  \item{object}{an object from which is extracted
  annotation tracks}

  \item{...}{extra arguments to allow extensions and passed
  to the next method call. For \code{atrack}, arguments in
  \code{...} are concatenated into a single
  \code{annotationTrack} object.}

  \item{x}{an R object}

  \item{value}{replacement value for the complete
  annotation data list}

  \item{default.margin}{margin to use if no margin data is
  stored in the \code{x}.}

  \item{data}{object used to extend the annotation track
  within a given data context. It is typically a
  matrix-like object, against which annotation
  specifications are matched using
  \code{\link{match_atrack}}.}

  \item{order}{an integer vector that indicates the order
  of the annotation tracks in the result list}

  \item{enforceNames}{logical that indicates if missing
  track names should be generated as \code{X<i>}}

  \item{.DATA}{data used to match and extend annotation
  specifications. It is passed to argument \code{data} of
  the \code{.atrack} methods, which in turn use pass it to
  \code{\link{match_atrack}}.}
}
\description{
  \code{.atrack} is an S4 generic method that converts an
  object into an annotation track object. It provides a
  general and flexible annotation framework that is used by
  \code{\link{aheatmap}} to annotates heatmap rows and
  columns.

  \code{is.atrack} tests if an object is an
  \code{annotationTrack} object.

  \code{adata} get/sets the annotation parameters on an
  object

  \code{amargin} get/sets the annotation margin, i.e. along
  which dimension of the data the annotations are to be
  considered.

  \code{anames} returns the reference margin names for
  annotation tracks, from their embedded annotation data
  object.

  \code{alength} returns the reference length for
  annotation tracks, from their embedded annotation data
  object

  \code{atrack} creates/concatenates \code{annotationTrack}
  objects

  \code{annotationTrack} is constructor function for
  \code{annotationTrack} object
}
\details{
  Methods for \code{.atrack} exist for common type of
  objects, which should provide enough options for new
  methods to define how annotation track are extracted from
  more complex objects, by coercing/filtering them into a
  supported type.
}
\section{Methods}{
  \describe{

  \item{.atrack}{\code{signature(object = "ANY")}: The
  default method converts character or integer vectors into
  factors. Numeric vectors, factors, a single NA or
  \code{annotationTrack} objects are returned unchanged
  (except from reordering by argument \code{order}). Data
  frames are not changed either, but class
  'annotationTrack' is appended to their original class
  set. }

  }
}
\keyword{internal}
\keyword{methods}

