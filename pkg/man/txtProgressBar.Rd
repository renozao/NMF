\name{txtProgressBar}
\alias{txtProgressBar}
\title{Simple Progress Bar}
\usage{
  txtProgressBar(min = 0, max = 1, initial = 0, char = "=",
    width = NA, title = if (style == 3) " ", label,
    style = 1, file = "", shared = NULL)
}
\arguments{
  \item{shared}{specification of a shared directory to use
  when the progress bar is to be used by multiple
  processes.}

  \item{min}{(finite) numeric values for the extremes of
  the progress bar. Must have \code{min < max}.}

  \item{max}{(finite) numeric values for the extremes of
  the progress bar. Must have \code{min < max}.}

  \item{initial}{initial or new value for the progress bar.
  See \sQuote{Details} for what happens with invalid
  values.}

  \item{char}{the character (or character string) to form
  the progress bar.}

  \item{width}{the width of the progress bar, as a multiple
  of the width of \code{char}.  If \code{NA}, the default,
  the number of characters is that which fits into
  \code{getOption("width")}.}

  \item{title}{ignored, for compatibility with other
  progress bars.}

  \item{label}{ignored, for compatibility with other
  progress bars.}

  \item{style}{the \sQuote{style} of the bar -- see
  \sQuote{Details}.}

  \item{file}{an open connection object or \code{""} which
  indicates the console: \code{\link{stderr}()} might be
  useful here.}
}
\description{
  Creates a simple progress bar with title. This function
  is identical to \code{utils::txtProgressBar} but allow
  adding a title to the progress bar, and can be shared by
  multiple processes, e.g., in multicore or multi-hosts
  computations.
}
\author{
  R Core Team
}
\keyword{internal}

