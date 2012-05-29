\name{nmf.options}
\alias{nmf.getOption}
\alias{nmf.options}
\title{Package Specific Options}
\usage{
  nmf.options(...)

  nmf.getOption(x, default = NULL)
}
\arguments{
  \item{...}{any options can be defined, using \code{name =
  value} or by passing a list of such tagged values.
  However, only the ones below are used in base \R.
  Further, \code{options('name') == options()['name']}, see
  the example.  }

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
}

