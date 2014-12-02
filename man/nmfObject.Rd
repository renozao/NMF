\name{nmfObject}
\alias{nmfObject}
\title{Updating NMF Objects}
\usage{
  nmfObject(object, verbose = FALSE)
}
\arguments{
  \item{object}{an R object created by the NMF package, e.g., an object of class
\code{\linkS4class{NMF}} or \code{\linkS4class{NMFfit}}.}

  \item{verbose}{logical to toggle verbose messages.}
}
\description{
  This function serves to update an objects created with previous versions of the
NMF package, which would otherwise be incompatible with the current version,
due to changes in their S4 class definition.
}
\details{
  This function makes use of heuristics to automatically update object slots,
which have been borrowed from the BiocGenerics package, the function
\code{updateObjectFromSlots} in particular.
}

