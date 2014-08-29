\docType{methods}
\name{plot,NMFfit,missing-method}
\alias{plot,NMFfit,missing-method}
\title{Plots the residual track computed at regular interval during the fit of
the NMF model \code{x}.}
\usage{
  \S4method{plot}{NMFfit,missing}(x, y, skip = -1, ...)
}
\arguments{
  \item{skip}{an integer that indicates the number of points to skip/remove from the beginning
of the curve.
If \code{skip=1L} (default) only the initial residual -- that is computed before any iteration, is
skipped, if present in the track (it associated with iteration 0).}

  \item{x}{the coordinates of points in the plot. Alternatively, a
      single plotting structure, function or \emph{any \R object with a
        \code{plot} method} can be provided.}

  \item{y}{the y coordinates of points in the plot, \emph{optional}
    if \code{x} is an appropriate structure.}

  \item{...}{Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
    Many methods will accept the following arguments:
    \describe{
    \item{\code{type}}{what type of plot should be drawn.  Possible types are
      \itemize{
        \item \code{"p"} for \bold{p}oints,
        \item \code{"l"} for \bold{l}ines,
        \item \code{"b"} for \bold{b}oth,
        \item \code{"c"} for the lines part alone of \code{"b"},
        \item \code{"o"} for both \sQuote{\bold{o}verplotted},
        \item \code{"h"} for \sQuote{\bold{h}istogram} like (or
        \sQuote{high-density}) vertical lines,
        \item \code{"s"} for stair \bold{s}teps,
        \item \code{"S"} for other \bold{s}teps, see \sQuote{Details} below,
        \item \code{"n"} for no plotting.
      }
      All other \code{type}s give a warning or an error; using, e.g.,
      \code{type = "punkte"} being equivalent to \code{type = "p"} for S
      compatibility.  Note that some methods,
      e.g. \code{\link{plot.factor}}, do not accept this.
    }
    \item{\code{main}}{an overall title for the plot: see \code{\link{title}}.}
    \item{\code{sub}}{a sub title for the plot: see \code{\link{title}}.}
    \item{\code{xlab}}{a title for the x axis: see \code{\link{title}}.}
    \item{\code{ylab}}{a title for the y axis: see \code{\link{title}}.}
    \item{\code{asp}}{the \eqn{y/x} aspect ratio,
      see \code{\link{plot.window}}.}
    }
  }
}
\description{
  Plots the residual track computed at regular interval
  during the fit of the NMF model \code{x}.
}
\keyword{methods}

