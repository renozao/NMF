\name{NMFList-class}

\docType{class}
\alias{NMFList-class}

\alias{algorithm,NMFList-method}
\alias{runtime,NMFList-method}
\alias{show,NMFList-method}

\title{Class "NMFList" to Handle the Comparison of NMF Results}
\description{
Class "NMFList" is used to wrap into a list the results of different NMF runs
with the objective to compare them.

While it handles indifferently any kind of NMF result, it is usually used to 
compare NMF results from different algorithms.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{as.NMFList(...)}.
	  
}
\section{Slots}{
  \describe{    
    \item{\code{runtime}:}{Object of class \code{"proc_time"} that stores the 
    CPU time used to compute all the results in the list.}
    \item{\code{.Data}:}{Object of class \code{"list"}, inherited from 
    class \code{list}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\link[Biobase]{AssayData-class}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{algorithm}{: returns the character vector of the names of the algorithms 
     of the results in the list.}
    \item{plot}{ plots the error track of each result in the list onto a single plot.
    To work properly, the results must be computed with argument \code{.options='t'} or 
    \code{.options=list(track=TRUE, ...)}.
    See \code{\link{nmf-compare}} for more details.
    }
    \item{runtime}{: returns the CPU time used to compute all the results in the list, 
    as stored in slot \code{runtime}. If no time data is available then it returns \code{NULL}.}
    \item{show}{: show method for class \code{NMFList}.}
    \item{summary}{ returns a \code{data.frame} whose lines are the summary measures of each result
     in the list.}
	 }
}

\author{ Renaud Gaujoux  }

\seealso{
	\code{\link{nmf-compare}}, \code{\link{nmf-multiple}}, \code{\link{nmf}}
}
\examples{
showClass("NMFList")
}
\keyword{classes}