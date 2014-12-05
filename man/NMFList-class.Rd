\docType{class}
\name{NMFList-class}
\alias{NMFList-class}
\title{Class for Storing Heterogeneous NMF fits}
\description{
  This class wraps a list of NMF fit objects, which may come from different
runs of the function \code{\link{nmf}}, using different parameters, methods, etc..
These can be either from a single run (NMFfit) or multiple runs (NMFfitX).

Note that its definition/interface is very likely to change in the future.
}
\section{Methods}{
\describe{

\item{algorithm}{\code{signature(object = "NMFList")}: Returns the method names used to compute the NMF fits in the list.
It returns \code{NULL} if the list is empty.
}


\item{runtime}{\code{signature(object = "NMFList")}: Returns the CPU time required to compute all NMF fits in the list.
It returns \code{NULL} if the list is empty.
If no timing data are available, the sequential time is returned.
}


\item{seqtime}{\code{signature(object = "NMFList")}: Returns the CPU time that would be required to sequentially compute all NMF
fits stored in \code{object}.

This method calls the function \code{runtime} on each fit and sum up the
results.
It returns \code{NULL} on an empty object.
}


\item{show}{\code{signature(object = "NMFList")}: Show method for objects of class \code{NMFList}
}

}

}

