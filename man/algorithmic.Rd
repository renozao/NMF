\docType{methods}
\name{algorithmic-NMF}
\alias{algorithm}
\alias{algorithm<-}
\alias{algorithmic-NMF}
\alias{algorithm<--methods}
\alias{algorithm-methods}
\alias{algorithm<-,NMFfit,ANY-method}
\alias{algorithm,NMFfit-method}
\alias{algorithm,NMFfitXn-method}
\alias{algorithm<-,NMFSeed,function-method}
\alias{algorithm,NMFSeed-method}
\alias{algorithm<-,NMFStrategyFunction,function-method}
\alias{algorithm,NMFStrategyFunction-method}
\alias{algorithm<-,NMFStrategyOctave,character-method}
\alias{compare}
\alias{compare-methods}
\alias{compare,NMFfitXn-method}
\alias{logs}
\alias{logs,ANY-method}
\alias{logs-methods}
\alias{modelname}
\alias{modelname,ANY-method}
\alias{modelname-methods}
\alias{modelname,NMFfit-method}
\alias{modelname,NMFfitXn-method}
\alias{modelname,NMFStrategy-method}
\alias{niter}
\alias{niter<-}
\alias{niter<--methods}
\alias{niter-methods}
\alias{niter,NMFfit-method}
\alias{niter<-,NMFfit,numeric-method}
\alias{nrun}
\alias{nrun,ANY-method}
\alias{nrun-methods}
\alias{nrun,NMFfit-method}
\alias{nrun,NMFfitX1-method}
\alias{nrun,NMFfitX-method}
\alias{nrun,NMFfitXn-method}
\alias{objective}
\alias{objective<-}
\alias{objective<--methods}
\alias{objective-methods}
\alias{objective<-,NMFfit,ANY-method}
\alias{run}
\alias{run-methods}
\alias{runtime}
\alias{runtime.all}
\alias{runtime.all-methods}
\alias{runtime.all,NMFfit-method}
\alias{runtime.all,NMFfitX-method}
\alias{runtime-methods}
\alias{runtime,NMFfit-method}
\alias{seeding}
\alias{seeding<-}
\alias{seeding<--methods}
\alias{seeding-methods}
\alias{seeding<-,NMFfit-method}
\alias{seeding,NMFfit-method}
\alias{seeding,NMFfitXn-method}
\alias{seqtime}
\alias{seqtime-methods}
\alias{seqtime,NMFfitXn-method}
\alias{seqtime,NMFList-method}
\title{Generic Interface for Algorithms}
\usage{
  algorithm(object, ...)

  algorithm(object, ...)<-value

  seeding(object, ...)

  seeding(object, ...)<-value

  niter(object, ...)

  niter(object, ...)<-value

  nrun(object, ...)

  objective(object, ...)

  objective(object, ...)<-value

  runtime(object, ...)

  runtime.all(object, ...)

  seqtime(object, ...)

  modelname(object, ...)

  run(object, y, x, ...)

  logs(object, ...)

  compare(object, ...)
}
\arguments{
  \item{object}{an object computed using some algorithm, or that describes an algorithm
itself.}

  \item{value}{replacement value}

  \item{...}{extra arguments to allow extension}

  \item{y}{data object, e.g. a target matrix}

  \item{x}{a model object used as a starting point by the algorithm,
e.g. a non-empty NMF model.}
}
\description{
  The functions documented here are S4 generics that define
  an general interface for -- optimisation -- algorithms.

  This interface builds upon the broad definition of an
  algorithm as a workhorse function to which is associated
  auxiliary objects such as an underlying model or an
  objective function that measures the adequation of the
  model with observed data. It aims at complementing the
  interface provided by the \code{\link{stats}} package.
}
\details{
  \code{algorithm} and \code{algorithm<-} get/set an object
  that describes the algorithm used to compute another
  object, or with which it is associated. It may be a
  simple character string that gives the algorithm's names,
  or an object that includes the algorithm's definition
  itself (e.g. an \code{\link{NMFStrategy}} object).

  \code{seeding} get/set the seeding method used to
  initialise the computation of an object, i.e. usually the
  function that sets the starting point of an algorithm.

  \code{niter} and \code{niter<-} get/set the number of
  iterations performed to compute an object. The function
  \code{niter<-} would usually be called just before
  returning the result of an algorithm, when putting
  together data about the fit.

  \code{nrun} returns the number of times the algorithm has
  been run to compute an object. Usually this will be 1,
  but may be be more if the algorithm involves multiple
  starting points.

  \code{objective} and \code{objective<-} get/set the
  objective function associated with an object. Some
  methods for \code{objective} may also compute the
  objective value with respect to some target/observed
  data.

  \code{runtime} returns the CPU time required to compute
  an object. This would generally be an object of class
  \code{\link[=proc.time]{proc_time}}.

  \code{runtime.all} returns the CPU time required to
  compute a collection of objects, e.g. a sequence of
  independent fits.

  \code{seqtime} returns the sequential CPU time -- that
  would be -- required to compute a collection of objects.
  It would differ from \code{runtime.all} if the
  computations were performed in parallel.

  \code{modelname} returns a the type of model associated
  with an object.

  \code{run} calls the workhorse function that actually
  implements a strategy/algorithm, and run it on some data
  object.

  \code{logs} returns the log messages output during the
  computation of an object.

  \code{compare} compares objects obtained from running
  separate algorithms.
}
\section{Methods}{
\describe{

\item{algorithm}{\code{signature(object = "NMFfit")}: Returns the name of the algorithm that fitted the NMF model \code{object}.
}

\item{algorithm}{\code{signature(object = "NMFList")}: Returns the method names used to compute the NMF fits in the list.
It returns \code{NULL} if the list is empty.

See \code{\link{algorithm,NMFList-method}} for more details.
}

\item{algorithm}{\code{signature(object = "NMFfitXn")}: Returns the name of the common NMF algorithm used to compute all fits
stored in \code{object}

Since all fits are computed with the same algorithm, this method returns the
name of algorithm that computed the first fit.
It returns \code{NULL} if the object is empty.
}

\item{algorithm}{\code{signature(object = "NMFSeed")}: Returns the workhorse function of the seeding method described by \code{object}.
}

\item{algorithm}{\code{signature(object = "NMFStrategyFunction")}: Returns the single R function that implements the NMF algorithm -- as stored in
slot \code{algorithm}.
}

\item{algorithm}{\code{signature(object = "NMFStrategyOctave")}: Returns the name of the Octave/Matlab function that implements the NMF algorithm -- as stored in
slot \code{algorithm}.

See \code{\link{algorithm,NMFStrategyOctave-method}} for more details.
}


\item{algorithm<-}{\code{signature(object = "NMFSeed", value = "function")}: Sets the workhorse function of the seeding method described by \code{object}.
}

\item{algorithm<-}{\code{signature(object = "NMFStrategyFunction", value = "function")}: Sets the function that implements the NMF algorithm, stored in slot \code{algorithm}.
}

\item{algorithm<-}{\code{signature(object = "NMFStrategyOctave", value = "character")}: Sets the name of the Octave/Matlab function that implements the NMF algorithm.
It is stored in slot \code{algorithm}.
}


\item{compare}{\code{signature(object = "NMFfitXn")}: Compares the fits obtained by separate runs of NMF, in a single
call to \code{\link{nmf}}.
}


\item{logs}{\code{signature(object = "ANY")}: Default method that returns the value of attribute/slot \code{'logs'} or, if this latter
does not exists, the value of element \code{'logs'} if \code{object} is a \code{list}.
It returns \code{NULL} if no logging data was found.
}


\item{modelname}{\code{signature(object = "ANY")}: Default method which returns the class name(s) of \code{object}.
This should work for objects representing models on their own.

For NMF objects, this is the type of NMF model, that corresponds to the
name of the S4 sub-class of \code{\linkS4class{NMF}}, inherited by \code{object}.
}

\item{modelname}{\code{signature(object = "NMFfit")}: Returns the type of a fitted NMF model.
It is a shortcut for \code{modelname(fit(object)}.
}

\item{modelname}{\code{signature(object = "NMFfitXn")}: Returns the common type NMF model of all fits stored in \code{object}

Since all fits are from the same NMF model, this method returns the
model type of the first fit.
It returns \code{NULL} if the object is empty.
}

\item{modelname}{\code{signature(object = "NMFStrategy")}: Returns the model(s) that an NMF algorithm can fit.
}


\item{niter}{\code{signature(object = "NMFfit")}: Returns the number of iteration performed to fit an NMF model, typically
with function \code{\link{nmf}}.

Currently this data is stored in slot \code{'extra'}, but this might change
in the future.
}


\item{niter<-}{\code{signature(object = "NMFfit", value = "numeric")}: Sets the number of iteration performed to fit an NMF model.

This function is used internally by the function \code{\link{nmf}}.
It is not meant to be called by the user, except when developing
new NMF algorithms implemented as single function, to set the number
of iterations performed by the algorithm on the seed, before returning it
(see \code{\linkS4class{NMFStrategyFunction}}).
}


\item{nrun}{\code{signature(object = "ANY")}: Default method that returns the value of attribute \sQuote{nrun}.

Such an attribute my be attached to objects to keep track of data about
the parent fit object (e.g. by method \code{\link{consensus}}), which
can be used by subsequent function calls such as plot functions
(e.g. see \code{\link{consensusmap}}).
This method returns \code{NULL} if no suitable data was found.
}

\item{nrun}{\code{signature(object = "NMFfitX")}: Returns the number of NMF runs performed to create \code{object}.

It is a pure virtual method defined to ensure \code{nrun} is defined
for sub-classes of \code{NMFfitX}, which throws an error if called.

Note that because the \code{\link{nmf}} function allows to run the NMF
computation keeping only the best fit, \code{nrun} may return a value
greater than one, while only the result of the best run is stored in
the object (cf. option \code{'k'} in method \code{\link{nmf}}).
}

\item{nrun}{\code{signature(object = "NMFfit")}: This method always returns 1, since an \code{NMFfit} object is obtained
from a single NMF run.
}

\item{nrun}{\code{signature(object = "NMFfitX1")}: Returns the number of NMF runs performed, amongst which \code{object} was
selected as the best fit.
}

\item{nrun}{\code{signature(object = "NMFfitXn")}: Returns the number of runs performed to compute the fits stored in the list
(i.e. the length of the list itself).
}


\item{objective}{\code{signature(object = "NMFfit")}: Returns the objective function associated with the algorithm that computed the
fitted NMF model \code{object}, or the objective value with respect to a given
target matrix \code{y} if it is supplied.

See \code{\link{objective,NMFfit-method}} for more details.
}


\item{runtime}{\code{signature(object = "NMFfit")}: Returns the CPU time required to compute a single NMF fit.
}

\item{runtime}{\code{signature(object = "NMFList")}: Returns the CPU time required to compute all NMF fits in the list.
It returns \code{NULL} if the list is empty.
If no timing data are available, the sequential time is returned.

See \code{\link{runtime,NMFList-method}} for more details.
}


\item{runtime.all}{\code{signature(object = "NMFfit")}: Identical to \code{runtime}, since their is a single fit.
}

\item{runtime.all}{\code{signature(object = "NMFfitX")}: Returns the CPU time required to compute all the NMF runs.
It returns \code{NULL} if no CPU data is available.
}

\item{runtime.all}{\code{signature(object = "NMFfitXn")}: If no time data is available from in slot \sQuote{runtime.all} and argument
\code{null=TRUE}, then the sequential time as computed by
\code{\link{seqtime}} is returned, and a warning is thrown unless \code{warning=FALSE}.

See \code{\link{runtime.all,NMFfitXn-method}} for more details.
}


\item{seeding}{\code{signature(object = "NMFfit")}: Returns the name of the seeding method that generated the starting point
for the NMF algorithm that fitted the NMF model \code{object}.
}

\item{seeding}{\code{signature(object = "NMFfitXn")}: Returns the name of the common seeding method used the computation of all fits
stored in \code{object}

Since all fits are seeded using the same method, this method returns the
name of the seeding method used for the first fit.
It returns \code{NULL} if the object is empty.
}


\item{seqtime}{\code{signature(object = "NMFList")}: Returns the CPU time that would be required to sequentially compute all NMF
fits stored in \code{object}.

This method calls the function \code{runtime} on each fit and sum up the
results.
It returns \code{NULL} on an empty object.
}

\item{seqtime}{\code{signature(object = "NMFfitXn")}: Returns the CPU time that would be required to sequentially compute all NMF
fits stored in \code{object}.

This method calls the function \code{runtime} on each fit and sum up the
results.
It returns \code{NULL} on an empty object.
}

}

}
\section{Interface fo NMF algorithms}{

This interface is implemented for NMF algorithms by the classes \code{\link{NMFfit}},
\code{\link{NMFfitX}} and \code{\link{NMFStrategy}}, and their respective sub-classes.
The examples given in this documentation page are mainly based on this implementation.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

#----------
# modelname,ANY-method
#----------
# get the type of an NMF model
modelname(nmfModel(3))
modelname(nmfModel(3, model='NMFns'))
modelname(nmfModel(3, model='NMFOffset'))

#----------
# modelname,NMFStrategy-method
#----------
# get the type of model(s) associated with an NMF algorithm
modelname( nmfAlgorithm('brunet') )
modelname( nmfAlgorithm('nsNMF') )
modelname( nmfAlgorithm('offset') )
}
\keyword{methods}

