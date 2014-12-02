\name{tryViewport}
\alias{current.vpPath_patched}
\alias{tryViewport}
\alias{.use.grid.patch}
\title{Internal Grid Extension}
\usage{
  tryViewport(name, verbose = FALSE)

  current.vpPath_patched()

  .use.grid.patch()
}
\arguments{
  \item{name}{viewport name}

  \item{verbose}{toggle verbosity}
}
\description{
  These functions enable mixing base and grid graphics in \code{\link{aheatmap}},
by avoiding calls to the grid internal function \code{'L_gridDirty'}.
They are not exported (i.e. not tampering core functions) and are only meant for internal
use within the \pkg{NMF} package.

\code{tryViewport} tries to go down to a viewport in the current tree,
given its name.

\code{current.vpPath_patched} aims at substituting \code{\link[grid]{current.vpPath}},
so that the graphic engine is not reset.
This is essentially to prevent outputting a blank page at the beginning of PDF
graphic engines.

\code{.use.grid.patch} tells if the user enabled patching grid.
}
\details{
  \code{tryViewport} uses \code{\link[grid]{grid.ls}} and not
\code{\link{seekViewport}} as the latter would reset the graphic device
and break the mix grid/base graphic capability.
}
\keyword{internal}

