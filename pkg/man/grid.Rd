\name{grid.Call.graphics}
\alias{current.vpPath}
\alias{grid.Call}
\alias{grid.Call.graphics}
\alias{tryViewport}
\title{Grid Internal Functions}
\usage{
  grid.Call.graphics(fnname, ...)

  grid.Call(fnname, ...)

  current.vpPath()

  tryViewport(name, verbose = FALSE)
}
\arguments{
  \item{fnname}{native \code{grid} function name.}

  \item{name}{viewport name}

  \item{verbose}{toggle verbosity}
}
\description{
  These functions are generally redefinitions of their
  counter parts in the \code{grid} package, except that
  they do not call \code{'L_gridDirty'}. This is to enable
  mixing base and grid graphics in \code{\link{aheatmap}}.

  \code{tryViewport}It uses \code{\link[grid]{grid.ls}} and
  not \code{\link{seekViewport}} as the latter would reset
  the graphic device.
}
\keyword{internal}

