\docType{methods}
\name{nmf}
\alias{nmf}
\alias{nmf,data.frame,ANY,ANY-method}
\alias{nmf,formula,ANY,ANY-method}
\alias{nmf,matrix,data.frame,ANY-method}
\alias{nmf,matrix,matrix,ANY-method}
\alias{nmf,matrix,missing,ANY-method}
\alias{nmf,matrix,NMF,ANY-method}
\alias{nmf,matrix,NULL,ANY-method}
\alias{nmf,matrix,numeric,character-method}
\alias{nmf,matrix,numeric,function-method}
\alias{nmf,matrix,numeric,list-method}
\alias{nmf,matrix,numeric,missing-method}
\alias{nmf,matrix,numeric,NMFStrategy-method}
\alias{nmf,matrix,numeric,NULL-method}
\alias{nmf-methods}
\title{Running NMF algorithms}
\usage{
  nmf(x, rank, method, ...)

  \S4method{nmf}{matrix,numeric,NULL}(x, rank, method,
    seed = NULL, model = NULL, ...)

  \S4method{nmf}{matrix,numeric,function}(x, rank, method,
    seed, model = "NMFstd", ..., name,
    objective = "euclidean", mixed = FALSE)

  \S4method{nmf}{matrix,NMF,ANY}(x, rank, method, seed,
    ...)

  \S4method{nmf}{matrix,NULL,ANY}(x, rank, method, seed,
    ...)

  \S4method{nmf}{matrix,matrix,ANY}(x, rank, method, seed,
    model = list(), ...)

  \S4method{nmf}{formula,ANY,ANY}(x, rank, method, ...,
    model = NULL)

  \S4method{nmf}{matrix,numeric,NMFStrategy}(x, rank,
    method, seed = nmf.getOption("default.seed"),
    rng = NULL, nrun = if (length(rank) > 1) 30 else 1,
    model = NULL, .options = list(),
    .pbackend = nmf.getOption("pbackend"),
    .callback = NULL, ...)
}
\arguments{
  \item{x}{target data to fit, i.e. a matrix-like object}

  \item{rank}{specification of the factorization rank. It
  is usually a single numeric value, but other type of
  values are possible (e.g. matrix), for which specific
  methods are implemented. See for example methods
  \code{nmf,matrix,matrix,ANY}.

  If \code{rank} is a numeric vector with more than one
  element, e.g. a range of ranks, then \code{\link{nmf}}
  performs the estimation procedure described in
  \code{\link{nmfEstimateRank}}.}

  \item{method}{specification of the NMF algorithm. The
  most common way of specifying the algorithm is to pass
  the access key (i.e. a character string) of an algorithm
  stored in the package's dedicated registry, but methods
  exists that handle other types of values, such as
  \code{function} or \code{list} object. See their
  descriptions in section \emph{Methods}.

  If \code{method} is missing the algorithm to use is
  obtained from the option
  \code{nmf.getOption('default.algorithm')}, unless it can
  be infer from the type of NMF model to fit, if this later
  is available from other arguments. Factory fresh default
  value is \sQuote{brunet}, which corresponds to the
  standard NMF algorithm from \cite{Brunet2004} (see
  section \emph{Algorithms}).

  Cases where the algorithm is inferred from the call are
  when an NMF model is passed in arguments \code{rank} or
  \code{seed} (see description for
  \code{nmf,matrix,numeric,NULL} in section
  \emph{Methods}).}

  \item{...}{extra arguments to allow extension of the
  generic. Arguments that are not used in the chain of
  internal calls to \code{nmf} methods are passed to the
  function that effectively implements the algorithm that
  fits an NMF model on \code{x}.}

  \item{name}{name associated with the NMF algorithm
  implemented by the function \code{method} [only used when
  \code{method} is a function].}

  \item{objective}{specification of the objective function
  associated with the algorithm implemented by the function
  \code{method} [only used when \code{method} is a
  function].

  It may be either \code{'euclidean'} or \code{'KL'} for
  specifying the euclidean distance (Frobenius norm) or the
  Kullback-Leibler divergence respectively, or a function
  with signature \code{(x="NMF", y="matrix", ...)} that
  computes the objective value for an NMF model \code{x} on
  a target matrix \code{y}, i.e. the residuals between the
  target matrix and its NMF estimate. Any extra argument
  may be specified, e.g. \code{function(x, y, alpha,
  beta=2, ...)}.}

  \item{mixed}{a logical that indicates if the algorithm
  implemented by the function \code{method} support
  mixed-sign target matrices, i.e. that may contain
  negative values [only used when \code{method} is a
  function].}

  \item{seed}{specification of the starting point or
  seeding method, which will compute a starting point,
  usually using data from the target matrix in order to
  provide a good guess.

  The seeding method may be specified in the following way:

  \describe{

  \item{a \code{character} string:}{ giving the name of a
  \emph{registered} seeding method. The corresponding
  method will be called to compute the starting point.

  Available methods can be listed via \code{nmfSeed()}. See
  its dedicated documentation for details on each available
  registered methods (\code{\link{nmfSeed}}). }

  \item{a \code{list}:}{ giving the name of a
  \emph{registered} seeding method and, optionally, extra
  parameters to pass to it.}

  \item{a single \code{numeric}:}{ that is used to seed the
  random number generator, before generating a random
  starting point.

  Note that when performing multiple runs, the L'Ecuyer's
  RNG is used in order to produce a sequence of random
  streams, that is used in way that ensures that parallel
  computation are fully reproducible. }

  \item{an object that inherits from
  \code{\linkS4class{NMF}}:}{ it should contain the data of
  an initialised NMF model, i.e. it must contain valid
  basis and mixture coefficient matrices, directly usable
  by the algorithm's workhorse function.}

  \item{a \code{function}:}{ that computes the starting
  point. It must have signature \code{(object="NMF",
  target="matrix", ...)} and return an object that inherits
  from class \code{NMF}. It is recommended to use argument
  \code{object} as a template for the returned object, by
  only updating the basis and coefficient matrices, using
  \code{\link{basis<-}} and \code{\link{coef<-}}
  respectively. }

  }}

  \item{rng}{rng specification for the run(s). This
  argument should be used to set the the RNG seed, while
  still specifying the seeding method argument \var{seed}.}

  \item{model}{specification of the type of NMF model to
  use.

  It is used to instantiate the object that inherits from
  class \code{\linkS4class{NMF}}, that will be passed to
  the seeding method. The following values are supported:
  \itemize{

  \item \code{NULL}, the default model associated to the
  NMF algorithm is instantiated and \code{...} is looked-up
  for arguments with names that correspond to slots in the
  model class, which are passed to the function
  \code{\link{nmfModel}} to instantiate the model.
  Arguments in \code{...} that do not correspond to slots
  are passed to the algorithm.

  \item a single \code{character} string, that is the name
  of the NMF model class to be instantiate.  In this case,
  arguments in \code{...} are handled in the same way as
  when \code{model} is \code{NULL}.

  \item a \code{list} that contains named values that are
  passed to the function \code{\link{nmfModel}} to
  instantiate the model.  In this case, \code{...} is not
  looked-up at all, and passed entirely to the algorithm.
  This means that all necessary model parameters must be
  specified in \code{model}.

  }

  \strong{Argument/slot conflicts:} In the case a parameter
  of the algorithm has the same name as a model slot, then
  \code{model} MUST be a list -- possibly empty --, if one
  wants this parameter to be effectively passed to the
  algorithm.

  If a variable appears in both arguments \code{model} and
  \code{\dots}, the former will be used to initialise the
  NMF model, the latter will be passed to the NMF
  algorithm.  See code examples for an illustration of this
  situation.}

  \item{nrun}{number of runs to perform. It specifies the
  number of runs to perform. By default only one run is
  performed, except if \code{rank} is a numeric vector with
  more than one element, in which case a default of 30 runs
  per value of the rank are performed, allowing the
  computation of a consensus matrix that is used in
  selecting the appropriate rank (see
  \code{\link{consensus}}).

  When using a random seeding method, multiple runs are
  generally required to achieve stability and avoid
  \emph{bad} local minima.}

  \item{.options}{this argument is used to set runtime
  options.

  It can be a \code{list} containing named options with
  their values, or, in the case only boolean/integer
  options need to be set, a character string that specifies
  which options are turned on/off or their value, in a
  unix-like command line argument way.

  The string must be composed of characters that correspond
  to a given option (see mapping below), and modifiers '+'
  and '-' that toggle options on and off respectively. E.g.
  \code{.options='tv'} will toggle on options \code{track}
  and \code{verbose}, while \code{.options='t-v'} will
  toggle on option \code{track} and toggle off option
  \code{verbose}.

  Modifiers '+' and '-' apply to all option character found
  after them: \code{t-vp+k} means \code{track=TRUE},
  \code{verbose=parallel=FALSE}, and \code{keep.all=TRUE}.
  The default behaviour is to assume that \code{.options}
  starts with a '+'.

  for options that accept integer values, the value may be
  appended to the option's character e.g. \code{'p4'} for
  asking for 4 processors or \code{'v3'} for showing
  verbosity message up to level 3.

  The following options are available (the characters after
  \dQuote{-} are those to use to encode \code{.options} as
  a string): \describe{

  \item{debug - d}{ Toggle debug mode (default:
  \code{FALSE}). Like option \code{verbose} but with more
  information displayed.}

  \item{keep.all - k}{ used when performing multiple runs
  (\code{nrun}>1): if \code{TRUE}, all factorizations are
  saved and returned (default: \code{FALSE}). Otherwise
  only the factorization achieving the minimum residuals is
  returned.}

  \item{parallel - p}{ this option is useful on multicore
  *nix or Mac machine only, when performing multiple runs
  (\code{nrun} > 1) (default: \code{TRUE}). If \code{TRUE},
  the runs are performed using the parallel foreach backend
  defined in argument \code{.pbackend}. If this is set to
  \code{'mc'} or \code{'par'} then \code{nmf} tries to
  perform the runs using multiple cores with package
  \code{link[doParallel]{doParallel}} -- which therefore
  needs to be installed.

  If equal to an integer, then \code{nmf} tries to perform
  the computation on the specified number of processors.
  When passing options as a string the number is appended
  to the option's character e.g. \code{'p4'} for asking for
  4 processors.

  If \code{FALSE}, then the computation is performed
  sequentially using the base function
  \code{\link{sapply}}.

  Unlike option 'P' (capital 'P'), if the computation
  cannot be performed in parallel, then it will still be
  carried on sequentially.

  \strong{IMPORTANT NOTE FOR MAC OS X USERS:} The parallel
  computation is based on the \code{doMC} and
  \code{multicore} packages, so the same care should be
  taken as stated in the vignette of \code{doMC}:
  \emph{\dQuote{it is not safe to use doMC from R.app on
  Mac OS X. Instead, you should use doMC from a terminal
  session, starting R from the command line.}} }

  \item{parallel.required - P}{ Same as \code{p}, but an
  error is thrown if the computation cannot be performed in
  parallel or with the specified number of processors.}

  \item{restore.seed - r}{ deprecated option since version
  0.5.99. Will throw a warning if used.}

  \item{simplifyCB - S}{ toggle simplification of the
  callback results. Default is \code{TRUE}}

  \item{track - t}{ enables error tracking (default:
  FALSE). If \code{TRUE}, the returned object's slot
  \code{residuals} contains the trajectory of the objective
  values, which can be retrieved via \code{residuals(res,
  track=TRUE)} This tracking functionality is available for
  all built-in algorithms. }

  \item{verbose - v}{ Toggle verbosity (default:
  \code{FALSE}). If \code{TRUE}, messages about the
  configuration and the state of the current run(s) are
  displayed. The level of verbosity may be specified with
  an integer value, the greater the level the more messages
  are displayed. Value \code{FALSE} means no messages are
  displayed, while value \code{TRUE} is equivalent to
  verbosity level 1. }

  }}

  \item{.pbackend}{specification of the
  \code{\link{foreach}} parallel backend to register and/or
  use when running in parallel mode. See options \code{p}
  and \code{P} in argument \code{.options} for how to
  enable this mode. Note that any backend that is
  internally registered is cleaned-up on exit, so that the
  calling foreach environment should not be affected by a
  call to \code{nmf} -- except when \code{.pbackend=NULL}.

  Currently it accepts the following values: \describe{

  \item{\sQuote{par}}{ use the backend(s) defined by the
  package \code{\link{doParallel}};} \item{a numeric
  value}{ use the specified number of cores with
  \code{doParallel} backend;} \item{\sQuote{seq}}{ use the
  foreach sequential backend \code{doSEQ};}
  \item{\code{NULL}}{ use currently registered backend;}
  \item{\code{NA}}{ do not compute using a foreach loop --
  and therefore not in parallel -- but rather use a call to
  standard \code{\link{sapply}}. This is useful for when
  developing/debugging NMF algorithms, as foreach loop
  handling may sometime get in the way.

  Note that this is equivalent to using
  \code{.options='-p'} or \code{.options='p0'}, but takes
  precedence over any option specified in \code{.options}:
  e.g. \code{nmf(..., .options='P10', .pbackend=NA)}
  performs all runs sequentially using \code{sapply}. Use
  \code{nmf.options(pbackend=NA)} to completely disable
  foreach/parallel computations for all subsequent
  \code{nmf} calls.}

  \item{\sQuote{mc}}{ identical to \sQuote{par} and defined
  to ensure backward compatibility.} }}

  \item{.callback}{Used when option \code{keep.all=FALSE}
  (default).  It allows to pass a callback function that is
  called after each run when performing multiple runs (i.e.
  with \code{nrun>1}). This is useful for example if one is
  also interested in saving summary measures or process the
  result of each NMF fit before it gets discarded. After
  each run, the callback function is called with two
  arguments, the \code{\linkS4class{NMFfit}} object that as
  just been fitted and the run number: \code{.callback(res,
  i)}. For convenience, a function that takes only one
  argument or has signature \code{(x, ...)} can still be
  passed in \code{.callback}. It is wrapped internally into
  a dummy function with two arguments, only the first of
  which is passed to the actual callback function (see
  example with \code{summary}).

  The call is wrapped into a tryCatch so that callback
  errors do not stop the whole computation (see below).

  The results of the different calls to the callback
  function are stored in a miscellaneous slot accessible
  using the method \code{$} for \code{NMFfit} objects:
  \code{res$.callback}. By default \code{nmf} tries to
  simplify the list of callback result using \code{sapply},
  unless option \code{'simplifyCB'} is \code{FASE}.

  If no error occurs \code{res$.callback} contains the list
  of values that resulted from the calling the callback
  function --, ordered as the fits. If any error occurs in
  one of the callback calls, then the whole computation is
  \strong{not} stopped, but the error message is stored in
  \code{res$.callback}, in place of the result.

  See the examples for sample code.}
}
\value{
  The returned value depends on the run mode:

  \item{Single run:}{An object of class
  \code{\linkS4class{NMFfit}}.}

  \item{Multiple runs, single method:}{When \code{nrun > 1}
  and \code{method} is not \code{list}, this method returns
  an object of class \code{\linkS4class{NMFfitX}}.}

  \item{Multiple runs, multiple methods:}{When \code{nrun >
  1} and \code{method} is a \code{list}, this method
  returns an object of class \code{\linkS4class{NMFList}}.}
}
\description{
  The function \code{nmf} is a S4 generic defines the main
  interface to run NMF algorithms within the framework
  defined in package \code{NMF}. It has many methods that
  facilitates applying, developing and testing NMF
  algorithms.

  The package vignette \code{vignette('NMF')} contains an
  introduction to the interface, through a sample data
  analysis.
}
\details{
  The \code{nmf} function has multiple methods that compose
  a very flexible interface allowing to: \itemize{ \item
  combine NMF algorithms with seeding methods and/or
  stopping/convergence criterion at runtime;

  \item perform multiple NMF runs, which are computed in
  parallel whenever the host machine allows it;

  \item run multiple algorithms with a common set of
  parameters, ensuring a consistent environment (notably
  the RNG settings). }

  The workhorse method is
  \code{nmf,matrix,numeric,NMFStrategy}, which is
  eventually called by all other methods. The other methods
  provides convenient ways of specifying the NMF
  algorithm(s), the factorization rank, or the seed to be
  used. Some allow to directly run NMF algorithms on
  different types of objects, such as \code{data.frame} or
  \code{\link[Biobase]{ExpressionSet}} objects.
}
\section{Methods}{
  \describe{

  \item{nmf}{\code{signature(x = "data.frame", rank =
  "ANY", method = "ANY")}: Fits an NMF model on a
  \code{data.frame}.

  The target \code{data.frame} is coerced into a matrix
  with \code{\link{as.matrix}}. }

  \item{nmf}{\code{signature(x = "matrix", rank =
  "numeric", method = "NULL")}: Fits an NMF model using an
  appropriate algorithm when \code{method} is not supplied.

  This method tries to select an appropriate algorithm
  amongst the NMF algorithms stored in the internal
  algorithm registry, which contains the type of NMF models
  each algorithm can fit. This is possible when the type of
  NMF model to fit is available from argument \code{seed},
  i.e. if it is an NMF model itself. Otherwise the
  algorithm to use is obtained from
  \code{nmf.getOption('default.algorithm')}.

  This method is provided for internal usage, when called
  from other \code{nmf} methods with argument \code{method}
  missing in the top call (e.g.
  \code{nmf,matrix,numeric,missing}). }

  \item{nmf}{\code{signature(x = "matrix", rank =
  "numeric", method = "list")}: Fits multiple NMF models on
  a common matrix using a list of algorithms.

  The models are fitted sequentially with \code{nmf} using
  the same options and parameters for all algorithms. In
  particular, irrespective of the way the computation is
  seeded, this method ensures that all fits are performed
  using the same initial RNG settings.

  This method returns an object of class
  \code{\linkS4class{NMFList}}, that is essentially a list
  containing each fit. }

  \item{nmf}{\code{signature(x = "matrix", rank =
  "numeric", method = "character")}: Fits an NMF model on
  \code{x} using an algorithm registered with access key
  \code{method}.

  Argument \code{method} is partially match against the
  access keys of all registered algorithms (case
  insensitive). Available algorithms are listed in section
  \emph{Algorithms} below or the introduction vignette. A
  vector of their names may be retrieved via
  \code{nmfAlgorithm()}. }

  \item{nmf}{\code{signature(x = "matrix", rank =
  "numeric", method = "function")}: Fits an NMF model on
  \code{x} using a custom algorithm defined the function
  \code{method}.

  The supplied function must have signature
  \code{(x=matrix, start=NMF, ...)} and return an object
  that inherits from class \code{\linkS4class{NMF}}. It
  will be called internally by the workhorse \code{nmf}
  method, with an NMF model to be used as a starting point
  passed in its argument \code{start}.

  Extra arguments in \code{...} are passed to \code{method}
  from the top \code{nmf} call. Extra arguments that have
  no default value in the definition of the function
  \code{method} are required to run the algorithm (e.g. see
  argument \code{alpha} of \code{myfun} in the examples).

  If the algorithm requires a specific type of NMF model,
  this can be specified in argument \code{model} that is
  handled as in the workhorse \code{nmf} method (see
  description for this argument). }

  \item{nmf}{\code{signature(x = "matrix", rank = "NMF",
  method = "ANY")}: Fits an NMF model using the NMF model
  \code{rank} to seed the computation, i.e. as a starting
  point.

  This method is provided for convenience as a shortcut for
  \code{nmf(x, nbasis(object), method, seed=object, ...)}
  It discards any value passed in argument \code{seed} and
  uses the NMF model passed in \code{rank} instead. It
  throws a warning if argument \code{seed} not missing.

  If \code{method} is missing, this method will call the
  method \code{nmf,matrix,numeric,NULL}, which will infer
  an algorithm suitable for fitting an NMF model of the
  class of \code{rank}. }

  \item{nmf}{\code{signature(x = "matrix", rank = "NULL",
  method = "ANY")}: Fits an NMF model using the NMF model
  supplied in \code{seed}, to seed the computation, i.e. as
  a starting point.

  This method is provided for completeness and is
  equivalent to \code{nmf(x, seed, method, ...)}. }

  \item{nmf}{\code{signature(x = "matrix", rank =
  "missing", method = "ANY")}: Method defined to ensure the
  correct dispatch to workhorse methods in case of argument
  \code{rank} is missing. }

  \item{nmf}{\code{signature(x = "matrix", rank =
  "numeric", method = "missing")}: Method defined to ensure
  the correct dispatch to workhorse methods in case of
  argument \code{method} is missing. }

  \item{nmf}{\code{signature(x = "matrix", rank = "matrix",
  method = "ANY")}: Fits an NMF model partially seeding the
  computation with a given matrix passed in \code{rank}.

  The matrix \code{rank} is used either as initial value
  for the basis or mixture coefficient matrix, depending on
  its dimension.

  Currently, such partial NMF model is directly used as a
  seed, meaning that the remaining part is left
  uninitialised, which is not accepted by all NMF
  algorithm. This should change in the future, where the
  missing part of the model will be drawn from some random
  distribution.

  Amongst built-in algorithms, only \sQuote{snmf/l} and
  \sQuote{snmf/r} support partial seeds, with only the
  coefficient or basis matrix initialised respectively. }

  \item{nmf}{\code{signature(x = "matrix", rank =
  "data.frame", method = "ANY")}: Shortcut for \code{nmf(x,
  as.matrix(rank), method, ...)}. }

  \item{nmf}{\code{signature(x = "formula", rank = "ANY",
  method = "ANY")}: This method implements the interface
  for fitting formula-based NMF models. See
  \code{\link{nmfModel}}.

  Argument \code{rank} target matrix or formula
  environment. If not missing, \code{model} must be a
  \code{list}, a \code{data.frame} or an \code{environment}
  in which formula variables are searched for. }

  }
}
\section{Optimized C++ vs. plain R}{
  Lee and Seung's multiplicative updates are used by
  several NMF algorithms. To improve speed and memory
  usage, a C++ implementation of the specific matrix
  products is used whenever possible. It directly computes
  the updates for each entry in the updated matrix, instead
  of using multiple standard matrix multiplication.

  The algorithms that benefit from this optimization are:
  'brunet', 'lee', 'nsNMF' and 'offset'. % and 'lnmf'
  However there still exists plain R versions for these
  methods, which implement the updates as standard matrix
  products. These are accessible by adding the prefix '.R#'
  to their name: '.R#brunet', '.R#lee', '.R#nsNMF' and
  '.R#offset'.
}

\section{Algorithms}{
  All algorithms are accessible by their respective access
  key as listed below. The following algorithms are
  available: \describe{

  \item{\sQuote{brunet}}{ Standard NMF, based on the
  Kullback-Leibler divergence, from \cite{Brunet et al. (2004)}. It
  uses simple multiplicative updates from \cite{Lee et al. (2001)},
  enhanced to avoid numerical underflow.

  Default stopping criterion: invariance of the
  connectivity matrix (see
  \code{\link{nmf.stop.connectivity}}). }

  \item{\sQuote{lee}}{ Standard NMF based on the Euclidean
  distance from \cite{Lee et al. (2001)}. It uses simple
  multiplicative updates.

  Default stopping criterion: invariance of the
  connectivity matrix (see
  \code{\link{nmf.stop.connectivity}}). }

  \item{ls-nmf}{ Least-Square NMF from \cite{Wang et al. (2006)}. It
  uses modified versions of Lee and Seung's multiplicative
  updates for the Euclidean distance, which incorporates
  weights on each entry of the target matrix, e.g. to
  reflect measurement uncertainty.

  Default stopping criterion: stationarity of the objective
  function (see \code{\link{nmf.stop.stationary}}). }

  \item{\sQuote{nsNMF}}{ Nonsmooth NMF from
  \cite{Pascual-Montano et al. (2006)}. It uses a modified version of
  Lee and Seung's multiplicative updates for the
  Kullback-Leibler divergence \cite{Lee et al. (2001)}, to fit a
  extension of the standard NMF model, that includes an
  intermediate smoothing matrix, meant meant to produce
  sparser factors.

  Default stopping criterion: invariance of the
  connectivity matrix (see
  \code{\link{nmf.stop.connectivity}}). }

  \item{\sQuote{offset}}{ NMF with offset from
  \cite{Badea (2008)}. It uses a modified version of Lee and
  Seung's multiplicative updates for Euclidean distance
  \cite{Lee et al. (2001)}, to fit an NMF model that includes an
  intercept, meant to capture a common baseline and shared
  patterns, in order to produce cleaner basis components.

  Default stopping criterion: invariance of the
  connectivity matrix (see
  \code{\link{nmf.stop.connectivity}}). }

  \item{\sQuote{pe-nmf}}{ Pattern-Expression NMF from
  \emph{Zhang2008}. It uses multiplicative updates to
  minimize an objective function based on the Euclidean
  distance, that is regularized for effective expression of
  patterns with basis vectors.

  Default stopping criterion: stationarity of the objective
  function (see \code{\link{nmf.stop.stationary}}). }

  \item{\sQuote{snmf/r}, \sQuote{snmf/l}}{ Alternating
  Least Square (ALS) approach from \cite{Kim et al. (2007)}. It
  applies the nonnegative least-squares algorithm from
  \cite{Van Benthem et al. (2004)} (i.e. fast combinatorial
  nonnegative least-squares for multiple right-hand), to
  estimate the basis and coefficient matrices alternatively
  (see \code{\link{fcnnls}}). It minimises an
  Euclidean-based objective function, that is regularized
  to favour sparse basis matrices (for \sQuote{snmf/l}) or
  sparse coefficient matrices (for \sQuote{snmf/r}).

  Stopping criterion: built-in within the internal
  workhorse function \code{nmf_snmf}, based on the KKT
  optimality conditions. }

  }
}

\section{Seeding methods}{
  The purpose of seeding methods is to compute initial
  values for the factor matrices in a given NMF model. This
  initial guess will be used as a starting point by the
  chosen NMF algorithm.

  The seeding method to use in combination with the
  algorithm can be passed to interface \code{nmf} through
  argument \code{seed}. The seeding seeding methods
  available in registry are listed by the function
  \code{\link{nmfSeed}} (see list therein).

  Detailed examples of how to specify the seeding method
  and its parameters can be found in the \emph{Examples}
  section of this man page and in the package's vignette.
}
\examples{
# Only basic calls are presented in this manpage.
# Many more examples are provided in the demo file nmf.R
\dontrun{
demo('nmf')
}

# random data
x <- rmatrix(20,10)

# run default algorithm with rank 2
res <- nmf(x, 2)

# specify the algorithm
res <- nmf(x, 2, 'lee')

# get verbose message on what is going on
res <- nmf(x, 2, .options='v')
\dontrun{
# more messages
res <- nmf(x, 2, .options='v2')
# even more
res <- nmf(x, 2, .options='v3')
# and so on ...
}
}
\references{
  Badea L (2008). "Extracting gene expression profiles
  common to colon and pancreatic adenocarcinoma using
  simultaneous nonnegative matrix factorization." _Pacific
  Symposium on Biocomputing. Pacific Symposium on
  Biocomputing_, *290*, pp. 267-78. ISSN 1793-5091, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/18229692>.

  Brunet J, Tamayo P, Golub TR and Mesirov JP (2004).
  "Metagenes and molecular pattern discovery using matrix
  factorization." _Proceedings of the National Academy of
  Sciences of the United States of America_, *101*(12), pp.
  4164-9. ISSN 0027-8424, <URL:
  http://dx.doi.org/10.1073/pnas.0308531101>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/15016911>.

  Kim H and Park H (2007). "Sparse non-negative matrix
  factorizations via alternating non-negativity-constrained
  least squares for microarray data analysis."
  _Bioinformatics (Oxford, England)_, *23*(12), pp.
  1495-502. ISSN 1460-2059, <URL:
  http://dx.doi.org/10.1093/bioinformatics/btm134>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/17483501>.

  Lee DD and Seung H (2001). "Algorithms for non-negative
  matrix factorization." _Advances in neural information
  processing systems_. <URL:
  http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization\#0>.

  Pascual-Montano A, Carazo JM, Kochi K, Lehmann D and
  Pascual-marqui RD (2006). "Nonsmooth nonnegative matrix
  factorization (nsNMF)." _IEEE Trans. Pattern Anal. Mach.
  Intell_, *28*, pp. 403-415.

  Van Benthem M and Keenan MR (2004). "Fast algorithm for
  the solution of large-scale non-negativity-constrained
  least squares problems." _Journal of Chemometrics_,
  *18*(10), pp. 441-450. ISSN 0886-9383, <URL:
  http://dx.doi.org/10.1002/cem.889>, <URL:
  http://doi.wiley.com/10.1002/cem.889>.

  Wang G, Kossenkov AV and Ochs MF (2006). "LS-NMF: a
  modified non-negative matrix factorization algorithm
  utilizing uncertainty estimates." _BMC bioinformatics_,
  *7*, pp. 175. ISSN 1471-2105, <URL:
  http://dx.doi.org/10.1186/1471-2105-7-175>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/16569230>.
}
\seealso{
  \code{\link{nmfAlgorithm}}
}
\keyword{methods}

