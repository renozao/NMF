\docType{class}
\name{NMF-class}
\alias{$,NMF-method}
\alias{$<-,NMF-method}
\alias{.DollarNames,NMF-method}
\alias{NMF-class}
\alias{misc}
\title{Generic Interface for Nonnegative Matrix Factorisation Models}
\usage{
  misc(object, ...)

  \S4method{$}{NMF}(x, name)

  \S4method{$}{NMF}(x, name) <- value

  \S4method{.DollarNames}{NMF}(x, pattern = "")
}
\arguments{
  \item{object}{an object that inherit from class
  \code{NMF}}

  \item{...}{extra arguments (not used)}

  \item{x}{ object from which to extract element(s) or in
  which to replace element(s).  }

  \item{name}{ A literal character string or a \link{name}
  (possibly \link{backtick} quoted).  For extraction, this
  is normally (see under \sQuote{Environments}) partially
  matched to the \code{\link{names}} of the object.  }

  \item{value}{typically an array-like \R object of a
  similar class as \code{x}.}

  \item{pattern}{ A regular expression.  Only matching
  names are returned.  }
}
\description{
  The class \code{NMF} is a \emph{virtual class} that
  defines a common interface to handle Nonnegative Matrix
  Factorization models (NMF models) in a generic way.
  Provided a minimum set of generic methods is implemented
  by concrete model classes, these benefit from a whole set
  of functions and utilities to perform common computations
  and tasks in the context of Nonnegative Matrix
  Factorization.

  The function \code{misc} provides access to miscellaneous
  data members stored in slot \code{misc} (as a
  \code{list}), which allow extensions of NMF models to be
  implemented, without defining a new S4 class.
}
\details{
  Class \code{NMF} makes it easy to develop new models that
  integrate well into the general framework implemented by
  the \emph{NMF} package.

  Following a few simple guidelines, new types of NMF
  models benefit from all the functionalities available for
  the built-in NMF models -- that derive themselves from
  class \code{NMF}. See section \emph{Implementing NMF
  models} below.

  See \code{\linkS4class{NMFstd}}, and references and links
  therein for details on the built-in implementations of
  the standard NMF model and its extensions.
}
\section{Slots}{
  \describe{

  \item{misc}{A list that is used internally to temporarily
  store algorithm parameters during the computation.}

  }
}
\section{Methods}{
  \describe{

  \item{$}{\code{signature(x = "NMF")}: shortcut for
  \code{x@misc[[name, exact=TRUE]]} respectively. }

  \item{$}{\code{signature(x = "NMF")}: shortcut for
  \code{x@misc[[name, exact=TRUE]]} respectively. }

  \item{$<-}{\code{signature(x = "NMF")}: shortcut for
  \code{x@misc[[name]] <- value} }

  \item{$<-}{\code{signature(x = "NMF")}: shortcut for
  \code{x@misc[[name]] <- value} }

  \item{.DollarNames}{\code{signature(x = "NMF")}:
  Auto-completion for \code{\linkS4class{Octave}} objects }

  \item{.DollarNames}{\code{signature(x = "NMF")}:
  Auto-completion for \code{\linkS4class{Octave}} objects }

  \item{.basis}{\code{signature(object = "NMF")}: Pure
  virtual method for objects of class
  \code{\linkS4class{NMF}}, that should be overloaded by
  sub-classes, and throws an error if called. }

  \item{.basis<-}{\code{signature(object = "NMF", value =
  "matrix")}: Pure virtual method for objects of class
  \code{\linkS4class{NMF}}, that should be overloaded by
  sub-classes, and throws an error if called. }

  \item{.coef}{\code{signature(object = "NMF")}: Pure
  virtual method for objects of class
  \code{\linkS4class{NMF}}, that should be overloaded by
  sub-classes, and throws an error if called. }

  \item{.coef<-}{\code{signature(object = "NMF", value =
  "matrix")}: Pure virtual method for objects of class
  \code{\linkS4class{NMF}}, that should be overloaded by
  sub-classes, and throws an error if called. }

  \item{[}{\code{signature(x = "NMF", i = "ANY", j = "ANY",
  drop = "ANY")}: This method provides a convenient way of
  sub-setting objects of class \code{NMF}, using a
  matrix-like syntax.

  See \code{\link{[,NMF,ANY,ANY,ANY-method}} for more
  details. }

  \item{basiscor}{\code{signature(x = "NMF", y =
  "matrix")}: Computes the correlations between the basis
  vectors of \code{x} and the columns of \code{y}. }

  \item{basiscor}{\code{signature(x = "NMF", y = "NMF")}:
  Computes the correlations between the basis vectors of
  \code{x} and \code{y}. }

  \item{basiscor}{\code{signature(x = "NMF", y =
  "missing")}: Computes the correlations between the basis
  vectors of \code{x}. }

  \item{basismap}{\code{signature(object = "NMF")}: Plots a
  heatmap of the basis matrix of the NMF model
  \code{object}. This method also works for fitted NMF
  models (i.e. \code{NMFfit} objects). }

  \item{c}{\code{signature(x = "NMF")}: Binds compatible
  matrices and NMF models together. }

  \item{coefficients}{\code{signature(object = "NMF")}:
  Alias to \code{coef,NMF}, therefore also pure virtual. }

  \item{connectivity}{\code{signature(object = "NMF")}:
  Computes the connectivity matrix for an NMF model, for
  which cluster membership is given by the most
  contributing basis component in each sample. See
  \code{\link{predict,NMF-method}}. }

  \item{consensus}{\code{signature(object = "NMF")}: This
  method is provided for completeness and is identical to
  \code{\link{connectivity}}, and returns the connectivity
  matrix, which, in the case of a single NMF model, is also
  the consensus matrix. }

  \item{consensushc}{\code{signature(object = "NMF")}:
  Compute the hierarchical clustering on the connectivity
  matrix of \code{object}. }

  \item{consensusmap}{\code{signature(object = "NMF")}:
  Plots a heatmap of the connectivity matrix of an NMF
  model. }

  \item{deviance}{\code{signature(object = "NMF")}:
  Computes the distance between a matrix and the estimate
  of an \code{NMF} model. }

  \item{dim}{\code{signature(x = "NMF")}: method for NMF
  objects for the base generic \code{\link{dim}}. It
  returns all dimensions in a length-3 integer vector: the
  number of row and columns of the estimated target matrix,
  as well as the factorization rank (i.e. the number of
  basis components). }

  \item{dimnames}{\code{signature(x = "NMF")}: Returns the
  dimension names of the NMF model \code{x}.

  It returns either NULL if no dimnames are set on the
  object, or a 3-length list containing the row names of
  the basis matrix, the column names of the mixture
  coefficient matrix, and the column names of the basis
  matrix (i.e. the names of the basis components). }

  \item{dimnames<-}{\code{signature(x = "NMF", value =
  "ANY")}: sets the dimension names of the NMF model
  \code{x}.

  \code{value} can be \code{NULL} which resets all
  dimension names, or a 1, 2 or 3-length list providing
  names at least for the rows of the basis matrix.

  See \code{\link{dimnames<-,NMF,ANY-method}} for more
  details. }

  \item{extractFeatures}{\code{signature(object = "NMF")}:
  Select basis-specific features from an NMF model, by
  applying the method \code{extractFeatures,matrix} to its
  basis matrix. }

  \item{featureScore}{\code{signature(object = "NMF")}:
  Computes feature scores on the basis matrix of an NMF
  model. }

  \item{fitted}{\code{signature(object = "NMF")}: Pure
  virtual method for objects of class
  \code{\linkS4class{NMF}}, that should be overloaded by
  sub-classes, and throws an error if called. }

  \item{ibterms}{\code{signature(object = "NMF")}: Default
  pure virtual method that ensure a method is defined for
  concrete NMF model classes. }

  \item{icterms}{\code{signature(object = "NMF")}: Default
  pure virtual method that ensure a method is defined for
  concrete NMF model classes. }

  \item{loadings}{\code{signature(x = "NMF")}: Method
  loadings for NMF Models

  The method \code{loadings} is identical to \code{basis},
  but do not accept any extra argument.

  See \code{\link{loadings,NMF-method}} for more details. }

  \item{metaHeatmap}{\code{signature(object = "NMF")}:
  Deprecated method that is substituted by
  \code{\link{coefmap}} and \code{\link{basismap}}. }

  \item{nmf.equal}{\code{signature(x = "NMF", y = "NMF")}:
  Compares two NMF models.

  Arguments in \code{...} are used only when
  \code{identical=FALSE} and are passed to
  \code{all.equal}. }

  \item{nmf.equal}{\code{signature(x = "NMF", y =
  "NMFfit")}: Compares two NMF models when at least one
  comes from a NMFfit object, i.e. an object returned by a
  single run of \code{\link{nmf}}. }

  \item{nmf.equal}{\code{signature(x = "NMF", y =
  "NMFfitX")}: Compares two NMF models when at least one
  comes from multiple NMF runs. }

  \item{nneg}{\code{signature(object = "NMF")}: Apply
  \code{nneg} to the basis matrix of an \code{\link{NMF}}
  object (i.e. \code{basis(object)}). All extra arguments
  in \code{...} are passed to the method
  \code{nneg,matrix}. }

  \item{predict}{\code{signature(object = "NMF")}: Default
  method for NMF models }

  \item{profcor}{\code{signature(x = "NMF", y = "matrix")}:
  Computes the correlations between the basis profiles of
  \code{x} and the rows of \code{y}. }

  \item{profcor}{\code{signature(x = "NMF", y = "NMF")}:
  Computes the correlations between the basis profiles of
  \code{x} and \code{y}. }

  \item{profcor}{\code{signature(x = "NMF", y =
  "missing")}: Computes the correlations between the basis
  profiles of \code{x}. }

  \item{rmatrix}{\code{signature(x = "NMF")}: Returns the
  target matrix estimate of the NMF model \code{x},
  perturbated by adding a random matrix generated using the
  default method of \code{rmatrix}: it is a equivalent to
  \code{fitted(x) + rmatrix(fitted(x), ...)}.

  This method can be used to generate random target
  matrices that depart from a known NMF model to a
  controlled extend. This is useful to test the robustness
  of NMF algorithms to the presence of certain types of
  noise in the data. }

  \item{rnmf}{\code{signature(x = "NMF", target =
  "numeric")}: Generates a random NMF model of the same
  class and rank as another NMF model.

  This is the workhorse method that is eventually called by
  all other methods. It generates an NMF model of the same
  class and rank as \code{x}, compatible with the
  dimensions specified in \code{target}, that can be a
  single or 2-length numeric vector, to specify a square or
  rectangular target matrix respectively.

  See \code{\link{rnmf,NMF,numeric-method}} for more
  details. }

  \item{rnmf}{\code{signature(x = "NMF", target =
  "missing")}: Generates a random NMF model of the same
  dimension as another NMF model.

  It is a shortcut for \code{rnmf(x, nrow(x), ncol(x),
  ...)}, which returns a random NMF model of the same class
  and dimensions as \code{x}. }

  \item{rposneg}{\code{signature(object = "NMF")}: Apply
  \code{rposneg} to the basis matrix of an
  \code{\link{NMF}} object. }

  \item{show}{\code{signature(object = "NMF")}: Show method
  for objects of class \code{NMF} }

  \item{sparseness}{\code{signature(x = "NMF")}: Compute
  the sparseness of an object of class \code{NMF}, as the
  sparseness of the basis and coefficient matrices computed
  separately.

  It returns the two values in a numeric vector with names
  \sQuote{basis} and \sQuote{coef}. }

  \item{summary}{\code{signature(object = "NMF")}: Computes
  summary measures for a single NMF model.

  The following measures are computed:

  See \code{\link{summary,NMF-method}} for more details. }

  }
}
\section{Implementing NMF models}{
  The class \code{NMF} only defines a basic data/low-level
  interface for NMF models, as a collection of generic
  methods, responsible with data handling, upon which
  relies a comprehensive set of functions, composing a rich
  higher-level interface.

  Actual NMF models are defined as sub-classes that
  inherits from class \code{NMF}, and implement the
  management of data storage, providing definitions for the
  interface's pure virtual methods.

  The minimum requirement to define a new NMF model that
  integrates into the framework of the \emph{NMF} package
  are the followings:

  \itemize{

  \item Define a class that inherits from class \code{NMF}
  and implements the new model, say class \code{myNMF}.

  \item Implement the following S4 methods for the new
  class \code{myNMF}: \describe{
  \item{fitted}{\code{signature(object = "myNMF", value =
  "matrix")}: Must return the estimated target matrix as
  fitted by the NMF model \code{object}.  }
  \item{basis}{\code{signature(object = "myNMF")}: Must
  return the basis matrix(e.g. the first matrix factor in
  the standard NMF model).  }
  \item{basis<-}{\code{signature(object = "myNMF", value =
  "matrix")}: Must return \code{object} with the basis
  matrix set to \code{value}.  }
  \item{coef}{\code{signature(object = "myNMF")}: Must
  return the matrix of mixture coefficients (e.g. the
  second matrix factor in the standard NMF model).  }
  \item{coef<-}{\code{signature(object = "myNMF", value =
  "matrix")}: Must return \code{object} with the matrix of
  mixture coefficients set to \code{value}.  } }

  The \emph{NMF} package provides "pure virtual"
  definitions of these methods for class \code{NMF} (i.e.
  with signatures \code{(object='NMF', ...)} and
  \code{(object='NMF', value='matrix')}) that throw an
  error if called, so as to force their definition for
  model classes.

  \item Optionally, implement method
  \code{rnmf}(signature(x="myNMF", target="ANY")). This
  method should call \code{callNextMethod(x=x,
  target=target, ...)} and fill the returned NMF model with
  its specific data suitable random values. }

  For concrete examples of NMF models implementations, see
  class \code{\linkS4class{NMFstd}} and its extensions
  (e.g. classes \code{\linkS4class{NMFOffset}} or
  \code{\linkS4class{NMFns}}).
}

\section{Creating NMF objects}{
  Strictly speaking, because class \code{NMF} is virtual,
  no object of class \code{NMF} can be instantiated, only
  objects from its sub-classes. However, those objects are
  sometimes shortly referred in the documentation and
  vignettes as "\code{NMF} objects" instead of "objects
  that inherits from class \code{NMF}".

  For built-in models or for models that inherit from the
  standard model class \code{\linkS4class{NMFstd}}, the
  factory method \code{nmfModel} enables to easily create
  valid \code{NMF} objects in a variety of common
  situations. See documentation for the the factory method
  \code{\link{nmfModel}} for more details.
}
\examples{
# show all the NMF models available (i.e. the classes that inherit from class NMF)
nmfModels()
# show all the built-in NMF models available
nmfModels(builtin.only=TRUE)

# class NMF is a virtual class so cannot be instantiated:
try( new('NMF') )

# To instantiate an NMF model, use the factory method nmfModel. see ?nmfModel
nmfModel()
nmfModel(3)
nmfModel(3, model='NMFns')
}
\references{
  Definition of Nonnegative Matrix Factorization in its
  modern formulation: \cite{Lee et al. (1999)}

  Historical first definition and algorithms: \cite{Paatero
  et al. (1994)}

  Lee DD and Seung HS (1999). "Learning the parts of
  objects by non-negative matrix factorization." _Nature_,
  *401*(6755), pp. 788-91. ISSN 0028-0836, <URL:
  http://dx.doi.org/10.1038/44565>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/10548103>.

  Paatero P and Tapper U (1994). "Positive matrix
  factorization: A non-negative factor model with optimal
  utilization of error estimates of data values."
  _Environmetrics_, *5*(2), pp. 111-126. <URL:
  http://dx.doi.org/10.1002/env.3170050203>, <URL:
  http://www3.interscience.wiley.com/cgi-bin/abstract/113468839/ABSTRACT>.
}
\seealso{
  Main interface to perform NMF in
  \code{\link{nmf-methods}}.

  Built-in NMF models and factory method in
  \code{\link{nmfModel}}.

  Method \code{\link{seed}} to set NMF objects with values
  suitable to start algorithms with.

  Other NMF-interface: \code{\link{.basis}},
  \code{\link{.basis<-}}, \code{\link{.coef}},
  \code{\link{.coef<-}}, \code{\link{basis}},
  \code{\link{basis<-}}, \code{\link{coef}},
  \code{\link{coef<-}}, \code{\link{coefficients}},
  \code{\link{loadings,NMF-method}},
  \code{\link{nmfModel}}, \code{\link{nmfModels}},
  \code{\link{rnmf}}, \code{\link{scoef}}
}
\keyword{methods}

