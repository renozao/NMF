\name{nmf.options}
\alias{nmf.getOption}
\alias{nmf.options}
\alias{nmf.printOptions}
\alias{nmf.resetOptions}
\title{NMF Package Specific Options}
\usage{
  nmf.options(...)

  nmf.getOption(x, default = NULL)

  nmf.resetOptions(..., ALL = FALSE)

  nmf.printOptions()
}
\arguments{
  \item{...}{option specifications. For \code{nmf.options}
  this can be named arguments or a single unnamed argument
  that is a named list (see \code{\link{options}}.

  For \code{nmf.resetOptions}, this must be the names of
  the options to reset. Note that \pkg{pkgmaker} version >=
  0.9.1 is required for this to work correctly, when
  options other than the default ones have been set after
  the package is loaded.}

  \item{ALL}{logical that indicates if options that are not
  part of the default set of options should be removed.
  Note that in \pkg{pkgmaker <= 0.9} this argument is only
  taken into account when no other argument is present.
  This is fixed in version 0.9.1.}

  \item{x}{a character string holding an option name.}

  \item{default}{if the specified option is not set in the
  options list, this value is returned.  This facilitates
  retrieving an option and checking whether it is set and
  setting it separately if not.}
}
\description{
  \code{nmf.options} sets/get single or multiple options,
  that are specific to the NMF package. It behaves in the
  same way as \code{\link[base]{options}}.

  \code{nmf.getOption} returns the value of a single
  option, that is specific to the NMF package. It behaves
  in the same way as \code{\link[base]{getOption}}.

  \code{nmf.resetOptions} reset all NMF specific options to
  their default values.

  \code{nmf.printOptions} prints all NMF specific options
  along with their default values, in a relatively compact
  way.
}
\examples{
# show all NMF specific options
nmf.printOptions()

# get some options
nmf.getOption('verbose')
nmf.getOption('backend')
# set new values
nmf.options(verbose=TRUE)
nmf.options(backend='mc', default.algorithm='lee')
nmf.printOptions()

# reset to default
nmf.resetOptions()
nmf.printOptions()
}

