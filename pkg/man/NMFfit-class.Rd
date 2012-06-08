\docType{class}
\name{NMFfit-class}
\alias{NMFfit}
\alias{NMFfit-class}
\title{Base Class for to store Nonnegative Matrix Factorisation results}
\usage{
  NMFfit(fit = nmfModel(), ..., rng = NULL)
}
\arguments{
  \item{fit}{an NMF model}

  \item{...}{extra argument used to initialise slots in the
  instantiating \code{NMFfit} object.}

  \item{rng}{RNG settings specification (typically a
  suitable value for \code{\link{.Random.seed}}).}
}
\description{
  Base class to handle the results of general
  \strong{Nonnegative Matrix Factorisation} algorithms
  (NMF).

  The function \code{NMFfit} is a factory method for NMFfit
  objects, that should not need to be called by the user.
  It is used internally by the functions \code{\link{nmf}}
  and \code{seed} to instantiate the starting point of NMF
  algorithms.
}
\details{
  It provides a general structure and generic functions to
  manage the results of NMF algorithms.  It contains a slot
  with the fitted NMF model (see slot \code{fit}) as well
  as data about the methods and parameters used to compute
  the factorization.

  The purpose of this class is to handle in a generic way
  the results of NMF algorithms. Its slot \code{fit}
  contains the fitted NMF model as an object of class
  \code{\linkS4class{NMF}}.

  Other slots contains data about how the factorization has
  been computed, such as the algorithm and seeding method,
  the computation time, the final residuals, etc\dots{}

  Class \code{NMFfit} acts as a wrapper class for its slot
  \code{fit}.  It inherits from interface class
  \code{\linkS4class{NMF}} defined for generic NMF models.
  Therefore, all the methods defined by this interface can
  be called directly on objects of class \code{NMFfit}. The
  calls are simply dispatched on slot \code{fit}, i.e.  the
  results are the same as if calling the methods directly
  on slot \code{fit}.
}
\section{Slots}{
  \describe{

  \item{fit}{An object that inherits from class
  \code{\linkS4class{NMF}}, and contains the fitted NMF
  model.

  NB: class \code{NMF} is a virtual class. The default
  class for this slot is \code{NMFstd}, that implements the
  standard NMF model.}

  \item{residuals}{A \code{numeric} vector that contains
  the final residuals or the residuals track between the
  target matrix and its NMF estimate(s).  Default value is
  \code{numeric()}.

  See method \code{\link{residuals}} for details on
  accessor methods and main interface \code{\link{nmf}} for
  details on how to compute NMF with residuals tracking.}

  \item{method}{a single \code{character} string that
  contains the name of the algorithm used to fit the model.
  Default value is \code{''}.}

  \item{seed}{a single \code{character} string that
  contains the name of the seeding method used to seed the
  algorithm that fitted the NMF model. Default value is
  \code{''}.  See \code{\link{nmf}} for more details.}

  \item{rng}{an object that contains the RNG settings used
  for the fit. Currently the settings are stored as an
  integer vector, the value of \code{\link{.Random.seed}}
  at the time the object is created. It is initialized by
  the \code{initialized} method. See \code{\link{getRNG}}
  for more details.}

  \item{distance}{either a single \code{"character"} string
  that contains the name of the built-in objective
  function, or a \code{function} that measures the
  residuals between the target matrix and its NMF estimate.
  See \code{\link{objective}} and
  \code{\link{deviance,NMF-method}}.}

  \item{parameters}{a \code{list} that contains the extra
  parameters -- usually specific to the algorithm -- that
  were used to fit the model.}

  \item{runtime}{object of class \code{"proc_time"} that
  contains various measures of the time spent to fit the
  model. See \code{\link[base]{system.time}}}

  \item{options}{a \code{list} that contains the options
  used to compute the object.}

  \item{extra}{a \code{list} that contains extra
  miscellaneous data for internal usage only. For example
  it can be used to store extra parameters or temporary
  data, without the need to explicitly extend the
  \code{NMFfit} class. Currently built-in algorithms only
  use this slot to store the number of iterations performed
  to fit the object.

  Data that need to be easily accessible by the end-user
  should rather be set using the methods \code{$<-} that
  sets elements in the \code{list} slot \code{misc} -- that
  is inherited from class \code{\linkS4class{NMF}}.}

  }
}
\section{Methods}{
  \describe{

  \item{algorithm}{\code{signature(object = "NMFfit")}:
  Returns the name of the algorithm that fitted the NMF
  model \code{object}. }

  \item{basis}{\code{signature(object = "NMFfit")}: Returns
  the basis matrix from an NMF model fitted with function
  \code{\link{nmf}}.

  It is a shortcut for \code{basis(fit(object), ...)},
  dispatching the call to the \code{basis} method of the
  actual NMF model. }

  \item{basis<-}{\code{signature(object = "NMFfit", value =
  "matrix")}: Sets the the basis matrix of an NMF model
  fitted with function \code{\link{nmf}}.

  It is a shortcut for \code{basis(fit(object)) <- value},
  dispatching the call to the \code{basis<-} method of the
  actual NMF model. It is not meant to be used by the user,
  except when developing new NMF algorithms, to update the
  basis matrix of the seed object before returning it. }

  \item{coef}{\code{signature(object = "NMFfit")}: Returns
  the the coefficient matrix from an NMF model fitted with
  function \code{\link{nmf}}.

  It is a shortcut for \code{coef(fit(object), ...)},
  dispatching the call to the \code{coef} method of the
  actual NMF model. }

  \item{coef<-}{\code{signature(object = "NMFfit", value =
  "matrix")}: Sets the the coefficient matrix of an NMF
  model fitted with function \code{\link{nmf}}.

  It is a shortcut for \code{coef(fit(object)) <- value},
  dispatching the call to the \code{coef<-} method of the
  actual NMF model. It is not meant to be used by the user,
  except when developing new NMF algorithms, to update the
  coefficient matrix in the seed object before returning
  it. }

  \item{compare}{\code{signature(object = "NMFfit")}:
  Compare multiple NMF fits passed as arguments. }

  \item{deviance}{\code{signature(object = "NMFfit")}:
  Returns the deviance of a fitted NMF model.

  This method returns the final residual value if the
  target matrix \code{y} is not supplied, or the
  approximation error between the fitted NMF model stored
  in \code{object} and \code{y}. In this case, the
  computation is performed using the objective function
  \code{method} if not missing, or the objective of the
  algorithm that fitted the model (stored in slot
  \code{'distance'}).

  See \code{\link{deviance,NMFfit-method}} for more
  details. }

  \item{fit}{\code{signature(object = "NMFfit")}: Returns
  the NMF model object stored in slot \code{'fit'}. }

  \item{fit<-}{\code{signature(object = "NMFfit", value =
  "NMF")}: Updates the NMF model object stored in slot
  \code{'fit'} with a new value. }

  \item{fitted}{\code{signature(object = "NMFfit")}:
  Computes and return the estimated target matrix from an
  NMF model fitted with function \code{\link{nmf}}.

  It is a shortcut for \code{fitted(fit(object), ...)},
  dispatching the call to the \code{fitted} method of the
  actual NMF model. }

  \item{minfit}{\code{signature(object = "NMFfit")}:
  Returns the object its self, since there it is the result
  of a single NMF run. }

  \item{modelname}{\code{signature(object = "NMFfit")}:
  Returns the type of a fitted NMF model. It is a shortcut
  for \code{modelname(fit(object)}. }

  \item{niter}{\code{signature(object = "NMFfit")}: Returns
  the number of iteration performed to fit an NMF model,
  typically with function \code{\link{nmf}}.

  Currently this data is stored in slot \code{'extra'}, but
  this might change in the future. }

  \item{niter<-}{\code{signature(object = "NMFfit", value =
  "numeric")}: Sets the number of iteration performed to
  fit an NMF model.

  This function is used internally by the function
  \code{\link{nmf}}. It is not meant to be called by the
  user, except when developing new NMF algorithms
  implemented as single function, to set the number of
  iterations performed by the algorithm on the seed, before
  returning it (see
  \code{\linkS4class{NMFStrategyFunction}}). }

  \item{nmf.equal}{\code{signature(x = "NMFfit", y =
  "NMF")}: Compares two NMF models when at least one comes
  from a NMFfit object, i.e. an object returned by a single
  run of \code{\link{nmf}}. }

  \item{nmf.equal}{\code{signature(x = "NMFfit", y =
  "NMFfit")}: Compares two fitted NMF models, i.e. objects
  returned by single runs of \code{\link{nmf}}. }

  \item{NMFfitX}{\code{signature(object = "NMFfit")}:
  Creates an \code{NMFfitX1} object from a single fit. This
  is used in \code{\link{nmf}} when only the best fit is
  kept in memory or on disk. }

  \item{nrun}{\code{signature(object = "NMFfit")}: This
  method always returns 1, since an \code{NMFfit} object is
  obtained from a single NMF run. }

  \item{objective}{\code{signature(object = "NMFfit")}:
  Returns the objective function associated with the
  algorithm that computed the fitted NMF model
  \code{object}, or the objective value with respect to a
  given target matrix \code{y} if it is supplied. }

  \item{plot}{\code{signature(x = "NMFfit", y =
  "missing")}: Plots the residual track computed at regular
  interval during the fit of the NMF model \code{x}. }

  \item{residuals}{\code{signature(object = "NMFfit")}:
  Returns the residuals -- track -- between the target
  matrix and the NMF fit \code{object}. }

  \item{runtime}{\code{signature(object = "NMFfit")}:
  Returns the CPU time required to compute a single NMF
  fit. }

  \item{runtime.all}{\code{signature(object = "NMFfit")}:
  Identical to \code{runtime}, since their is a single fit.
  }

  \item{seeding}{\code{signature(object = "NMFfit")}:
  Returns the name of the seeding method that generated the
  starting point for the NMF algorithm that fitted the NMF
  model \code{object}. }

  \item{show}{\code{signature(object = "NMFfit")}: Show
  method for objects of class \code{NMFfit} }

  \item{summary}{\code{signature(object = "NMFfit")}:
  Computes summary measures for a single fit from
  \code{\link{nmf}}.

  This method adds the folowing measures to the measures
  computed by the method \code{summary,NMF}:

  See \code{\link{summary,NMFfit-method}} for more details.
  }

  }
}
\examples{
# run default NMF algorithm on a random matrix
n <- 50; r <- 3; p <- 20
V <- rmatrix(n, p)
res <- nmf(V, r)

# result class is NMFfit
class(res)
isNMFfit(res)

# show result
res

# compute summary measures
summary(res, target=V)
}

