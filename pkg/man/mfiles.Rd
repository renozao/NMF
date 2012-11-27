\name{mfiles}
\alias{mfiles}
\title{M Files}
\usage{
  mfiles(..., pattern = "mfile_", dir = tempdir())
}
\arguments{
  \item{...}{specification of a .m files as character
  arguments. The elements of the vector can be either file
  paths or plain Octave/Matlab code, which are then written
  to disk in -- temporary -- .m files. Note that the paths
  do not need to correspond to existing files.}

  \item{dir}{existing directory where to write the .m files
  generated from the plain code elements of \var{x}.}

  \item{pattern}{a non-empty character vector giving the
  initial part of the name.}
}
\description{
  \code{mfiles} convert mfile specifications in to real
  paths to .m files that can be sourced with
  \code{\link[RcppOctave]{o_source}}.
}

