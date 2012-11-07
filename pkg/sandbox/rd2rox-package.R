

#' Advanced usage of package NMF
#' 
#' The functions documented here provide advanced functionalities useful when
#' developing within the framework implemented in the NMF package.
#' 

#' 
#' \describe{
#' 
#' \item{isNMFfit}{ tells if an object results from an NMF fit. That is it
#' checks if \code{object} inherits from class \code{\linkS4class{NMFfit}} or
#' form class \code{\linkS4class{NMFfitX}}, which are returned by the function
#' \code{\link{nmf}}.  If \code{object} is a \code{list} and
#' \code{recursive=TRUE}, then the check is performed on each element of the
#' list, and the return value is a vector (or a list if \code{object} is a list
#' of list) of the same length as \code{object}.  }
#' 
#' \item{is.nmf}{ tests if an object is an NMF model or a class that extends
#' the class NMF.  } }
#' 
#' @aliases NMF-advanced isNMFfit is.nmf
#' @param object any R object.
#' @param recursive if \code{TRUE} and \code{object} is a list then the check
#' is performed on each element of the list. Note that the recursivity only
#' applies in the case of lists that are not themselves NMFfit objects, unlike
#' \code{NMFfitXn} objects for which the result of \code{isNMFfit} will always
#' be \code{TRUE} (a single logical value).
#' 

#' @return For \code{isNMFfit}, a \code{logical} vector (or a list if
#' \code{object} is a list of list) of the same length as \code{object}.
#' @author Renaud Gaujoux
#' @seealso \code{\linkS4class{NMFfit}}, \code{\linkS4class{NMFfitX}}
#' @examples
#' 
#' 
#' 	# generate a random 50 x 10 matrix
#' 	V <- rmatrix(50, 10)
#' 	
#' 	# single run
#' 	res <- nmf(V, 3)
#' 	isNMFfit(res)
#' 	
#' 	# multiple runs - keeping single fit
#' 	resm <- nmf(V, 3, nrun=3)
#' 	isNMFfit(resm)
#' 	
#' 	# multiple runs - keeping all fits
#' 	resM <- nmf(V, 3, nrun=3, .opt='k') 
#' 	isNMFfit(resM)
#' 	
#' 	# with a list of results
#' 	isNMFfit(list(res, resm, resM, 'not a result'))
#' 	isNMFfit(list(res, list(resm, resM), 'not a result')) # list of list
#' 	isNMFfit(list(res, resm, resM, 'not a result'), recursive=FALSE)
#' 	
#' 
NULL




#' Layer to use the NMF package within Bioconductor
#' 
#' The package NMF provides an optional layer for working with common objects
#' and functions defined in the Bioconductor platform.
#' 
#' It provides:
#' 
#' \itemize{
#' 
#' \item computation functions that support \code{ExpressionSet} objects as
#' inputs.
#' 
#' \item aliases and methods for generic functions defined and widely used by
#' Bioconductor base packages.
#' 
#' \item specialized vizualization methods that adapt the titles and legend
#' using bioinformatics terminology.
#' 
#' \item functions to link the results with annotations, etc...
#' 
#' }
#' 
#' 
#' @name Bioconductor specific layer
#' @aliases bioc-NMF distance,ExpressionSet,NMF-method featureNames
#' featureNames,NMF-method featureNames<- featureNames<-,NMF-method sampleNames
#' sampleNames,NMF-method sampleNames<- sampleNames<-,NMF-method
#' sampleNames<-,NMF,ANY-method metagenes metagenes-methods
#' metagenes,NMF-method metagenes<- metagenes<-,NMF,matrix-method metaprofiles
#' metaprofiles-methods metaprofiles,NMF-method metaprofiles<-
#' metaprofiles<-,NMF,matrix-method nmeta nmeta-methods nmeta,NMF-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{distance}{\code{signature(target = "ExpressionSet", x = "NMF", method,
#' ...)}: returns the distance between the expression matrix and a NMF model,
#' according to a given measure.  If both argument \code{target} and \code{x}
#' are missing, this function returns the \code{function} defined by argument
#' \code{method}. The later can either be a \code{function} or a
#' \code{character} string that correspond to a registered distance metric. For
#' the moment only the metrics 'KL' and 'euclidean' are defined.
#' 
#' %See function \code{\link{distance}} for more details.  }
#' 
#' \item{featureNames}{\code{signature(object = "NMF")}: returns the row names
#' of the basis matrix, as \code{\link[=rownames-NMF]{rownames}}.  If is
#' defined for the generic function \code{\link[Biobase]{featureNames}} from
#' the \code{Biobase} package.  }
#' 
#' \item{featureNames<-}{\code{signature(object = "NMF", value = "ANY")}: sets
#' the row names of the basis matrix, as
#' \code{\link[=rownames-NMF]{rownames<-}} Argument \code{value} must be in a
#' format accepted by the \code{\link{rownames}} method defined for matrices.
#' If is defined for the generic function \code{\link[Biobase]{featureNames<-}}
#' from the \code{Biobase} package.  }
#' 
#' \item{metagenes}{\code{signature(object = "NMF")}: returns the metagenes
#' matrix according to the model defined in \code{object}.  It is an alias to
#' method \code{\link{basis}}.  }
#' 
#' \item{metagenes<-}{\code{signature(object = "NMF", value = "matrix")}: sets
#' the metagenes matrix in \code{object}, and returns the updated object.  It
#' is an alias to method \code{\link{basis<-}}.  }
#' 
#' \item{metaprofiles}{\code{signature(object = "NMF")}: returns the
#' metaprofiles matrix according to the model defined in \code{object}.  It is
#' an alias to method \code{\link{coef}}.  }
#' 
#' \item{metaprofiles<-}{\code{signature(object = "NMF", value = "matrix")}:
#' sets the metaprofiles matrix in \code{object}, and returns the updated
#' object.  It is an alias to method \code{\link{coef<-}}.  }
#' 
#' \item{nmeta}{\code{signature(object = "NMF")}: returns the number of
#' metagenes use in NMF model \code{object}.  It is an alias to
#' \code{\link{nbasis}}.  }
#' 
#' \item{sampleNames}{\code{signature(object = "NMF")}: returns the column
#' names of the mixture coefficient matrix, as
#' \code{\link[=colnames-NMF]{colnames}}.  If is defined for the generic
#' function \code{\link[Biobase]{sampleNames}} from the \code{Biobase} package.
#' }
#' 
#' \item{sampleNames<-}{\code{signature(object = "NMF", value = "ANY")}: sets
#' the columns names of the basis matrix, as
#' \code{\link[=colnames-NMF]{colnames<-}} Argument \code{value} must be in a
#' format accepted by the \code{\link{colnames}} method defined for matrices.
#' If is defined for the generic function \code{\link[Biobase]{sampleNames<-}}
#' from the \code{Biobase} package.  }
#' 
#' }
#' @seealso NMF, NMF-utils
#' @keywords methods
NULL





#' Compute the Correlation Between NMF Basis and Mixture Coefficient Matrices
#' 
#' Functions \code{basiscor} (resp. \code{profcor}) computes the correlation
#' matrix between the basis vectors (resp. the basis profiles) of two NMF
#' models, or of an NMF model and a given compatible matrix.
#' 
#' The computation uses the base function \code{\link{cor}}.
#' 
#' \describe{
#' 
#' \item{basiscor}{ computes the correlation matrix between the basis vectors
#' (i.e. the \emph{columns} of the basis vector matrix) or the \code{columns}
#' (for matrix arguments) of \var{x} and \var{y}.
#' 
#' The arguments' dimensions must be compatible, i.e same number of rows and
#' basis vectors (or columns for matrix arguments).
#' 
#' If \var{y} is missing, then the correlations are computed between \var{x}
#' and \code{y=x}.  }
#' 
#' \item{profcor}{: computes the correlation matrix between the basis profiles
#' (i.e. the \emph{rows} of the mixture coefficient matrix) or the rows (for
#' matrix arguments) of \var{x} and \var{y}.
#' 
#' The arguments' dimensions must be compatible, i.e same number of columns and
#' basis vectors (or rows for matrix arguments).
#' 
#' If \var{y} is missing, then the correlations are computed between \var{x}
#' and \code{y=x}.  }
#' 
#' }
#' 
#' @name Correlations between basis/profiles
#' @aliases basiscor basiscor-methods basiscor,NMF,NMF-method
#' basiscor,NMF,matrix-method basiscor,matrix,NMF-method
#' basiscor,NMF,missing-method profcor profcor-methods profcor,NMF,NMF-method
#' profcor,NMF,matrix-method profcor,matrix,NMF-method
#' profcor,NMF,missing-method
#' @docType methods
#' @param x A \code{matrix} or an object that inherits from class
#' \code{\linkS4class{NMF}}.
#' @param y A \code{matrix} or an object that inherits from class
#' \code{\linkS4class{NMF}}.  If \var{y} is missing then the correlations are
#' computed using \var{y=x}.
#' @param ...  extra arguments passed to function \code{\link{cor}}.
#' @author Renaud Gaujoux
#' @seealso \code{\link{basis}}, \code{\link[=coef]{coef,NMF-method}}
#' @keywords methods plot
#' @examples
#' 
#' 
#' # generate two random NMF model
#' a <- rnmf(100, 3, 20)
#' b <- rnmf(100, 3, 20)
#' 
#' # Compute auto-correlations
#' basiscor(a)
#' # Compute correlations with b
#' basiscor(a, b)
#' 
#' # try to recover the underlying NMF model 'a'
#' res <- nmf(fitted(a), 3)
#' 
#' # Compute correlations with the true model
#' basiscor(a, res)
#' profcor(a, res)
#' 
#' # Compute correlations with a random compatible matrix
#' W <- rmatrix(basis(a))
#' basiscor(a, W)
#' \dontshow{ basiscor(W, a) }
#' 
#' H <- rmatrix(coef(a))
#' profcor(a, H)
#' \dontshow{ profcor(H, a) }
#' 
#' 
NULL





#' Dimension names for NMF objects
#' 
#' The methods \code{dimnames}, \code{rownames}, \code{colnames} and
#' \code{basisnames} and their respective replacement form allow to get and set
#' the dimension names of the matrix factors in a NMF model.
#' 
#' They behave as their equivalent on \code{matrix} objects, and ensure that
#' the dimension names are handled in a consistent way on both factors --
#' especially \code{basisnames<-} which affects both matrix factors
#' simultaneously.
#' 
#' The methods \code{dimnames} and \code{basisnames} are implemented as S4
#' methods, while the methods \code{\link{rownames}} and \code{\link{colnames}}
#' are the default ones that make use of the result from \code{dimnames}.
#' 
#' \describe{
#' 
#' \item{basisnames, basisnames<-}{: returns (resp. simultaneously sets) the
#' names of the columns of the matrix of basis vectors and the rows of the
#' mixture coefficient matrix.  }
#' 
#' \item{colnames, colnames<-}{: returns/sets the names of the columns of the
#' mixture coefficient matrix.  Note that the standard arguments
#' \code{do.NULL}, \code{prefix} as described in \code{\link{colnames}} should
#' not be used (it will always return a character vector of length 1, which is
#' likely to be incorrect).  }
#' 
#' \item{rownames, rownames<-}{: returns/sets the names of the columns of the
#' basis vector matrix.  Note that the standard arguments \code{do.NULL},
#' \code{prefix} as described in \code{\link{rownames}} must not be used (it
#' will always return a character vector of length 1, which is likely to be
#' incorrect).  }
#' 
#' \item{dimnames}{\code{signature(x = "NMF")}: returns the dimension names of
#' the NMF model \code{x}.  It returns either NULL if no dimnames are set on
#' the object, or a 3-length list containing the row names of the basis matrix,
#' the column names of the mixture coefficient matrix, and the column names of
#' the basis matrix (i.e. the basis vector names).  }
#' 
#' \item{dimnames<-}{\code{signature(x = "NMF", value)}: sets the dimension
#' names of the NMF model \code{x}.  \code{value} can be \code{NULL} which
#' resets all dimension names, or a 1, 2 or 3-length list providing names at
#' least for the rows of the basis vector matrix.  The optional second element
#' of \code{value} (NULL if absent) is used to set the column names of the
#' mixture coefficient matrix.  The optional third element of \code{value}
#' (NULL if absent) is used to set both the column names of the basis vector
#' matrix and the row names of the mixture coefficient matrix.  }
#' 
#' }
#' 
#' @name dimnames-methods
#' @aliases dimnames-NMF rownames-NMF colnames-NMF basisnames
#' basisnames,NMF-method basisnames<- basisnames<-,NMF-method basisnames<-.NMF
#' dimnames,NMF-method dimnames<-,NMF-method dimnames<-.NMF
#' @docType methods
#' @param x an object of class \code{\linkS4class{NMF}}.
#' @param value a character vector, or \code{NULL} or, in the case of
#' \code{dimnames<-}, a list 2 or 3-length list of character vectors.  See
#' section \emph{Details} for more details.
#' @param ...  extra argument to pass to internal methods. Not used.
#' @keywords methods
#' @examples
#' 
#' 
#' # create a random NMF object
#' a <- rnmf(2, 5, 3)
#' 
#' # set dimensions
#' dims <- list( features=paste('f', 1:nrow(a), sep=''), samples=paste('s', 1:ncol(a), sep=''), basis=paste('b', 1:nbasis(a), sep='') )
#' dimnames(a) <- dims
#' dimnames(a)
#' basis(a)
#' coef(a)
#' 
#' # access the dimensions separately
#' rownames(a)
#' colnames(a)
#' basisnames(a)
#' 
#' # set only the first dimension (rows of basis): the other two dimnames are set to NULL
#' dimnames(a) <- dims[1]
#' dimnames(a)
#' basis(a)
#' coef(a)
#' 
#' # set only the two first dimensions (rows and columns of basis and coef respectively):
#' # the basisnames are set to NULL 
#' dimnames(a) <- dims[1:2]
#' dimnames(a)
#' basis(a)
#' 
#' # reset the dimensions
#' dimnames(a) <- NULL
#' dimnames(a)
#' basis(a)
#' coef(a)
#' 
#' # set each dimensions separately
#' rownames(a) <- paste('X', 1:nrow(a), sep='') # only affect rows of basis
#' basis(a)
#' 
#' colnames(a) <- paste('Y', 1:ncol(a), sep='') # only affect columns of coef
#' coef(a)
#' 
#' basisnames(a) <- paste('Z', 1:nbasis(a), sep='') # affect both basis and coef matrices
#' basis(a)
#' coef(a)
#' 
#' 
#' 
NULL





#' Dimension for NMF objects
#' 
#' The methods \code{dim}, \code{nrow}, \code{ncol} and \code{nbasis} return
#' the different dimensions of an NMF model.
#' 
#' The package NMF defines a S4 method \code{dim} for objects of class
#' \code{\linkS4class{NMF}} (resp. \code{\linkS4class{NMFfitXn}}). This allows
#' the base methods \code{\link{nrow}} and \code{\link{ncol}} to be used, to
#' get the number of rows and columns of the target matrix estimated by an NMF
#' model (resp. fit).
#' 
#' The method \code{nbasis} is a new S4 generic that returns the third element
#' of the vector returned \code{dim}.
#' 

#' 
#' \describe{
#' 
#' \item{dim}{ returns a 3-length vector containing the dimensions of the
#' target matrix fitted by \code{x} together with the factorization rank of
#' \code{x}.  e.g. It will return \code{c(2000, 30, 3)} for an \code{NMF} model
#' that fits a 2000 x 30 target matrix using 3 basis vectors.
#' 
#' For \code{\linkS4class{NMFfitXn}} objects it returns the common dimension of
#' the NMF problem fitted by all the runs -- based on the dimension of the
#' first fit in the list.
#' 
#' }
#' 
#' \item{nrow}{: returns the number of rows of the basis matrix.  It
#' corresponds to the number of rows of the fitted target matrix.  }
#' 
#' \item{ncol}{: returns the number of columns of the mixture coefficient
#' matrix.  It corresponds to the number of columns of the fitted target
#' matrix.  }
#' 
#' \item{nbasis}{: returns the number of columns of the basis matrix.  It
#' corresponds to the rank of the factorization.
#' 
#' When \code{x} is a matrix, then it returns the attribute 'nbasis', which is
#' for example attached to the consensus matrix returned by the method
#' \code{\link{consensus}}. This is used to keep track of data about the parent
#' fit and annotate plots.  }
#' 
#' }
#' 
#' @name dim-methods: Dimensions for NMF objects
#' @aliases nbasis nbasis,matrix-method nbasis,NMF-method
#' nbasis,NMFfitXn-method dim,NMF-method dim,NMFfitXn-method
#' @docType methods
#' @param x an object of class \code{\linkS4class{NMF}} or
#' \code{\linkS4class{NMFfitXn}}.
#' @return For \code{dim}, a 3-length integer vector.
#' 
#' For \code{nrow}, \code{ncol} and \code{nbasis} a single integer.
#' @seealso \code{\link[=dimnames,NMF-method]{dimnames}}
#' @keywords methods
#' @examples
#' 
#' 
#' # dimensions of an empty NMF model
#' dim( nmfModel() )
#' 
#' # dimensions of a random NMF model
#' x <- nmfModel(3, 100, 20)
#' dim(x)
#' nrow(x)
#' ncol(x)
#' nbasis(x)
#' 
#' # dimensions of a single NMF fit
#' V <- rmatrix(100, 20)
#' x <- nmf(V, 3)
#' dim(x)
#' nrow(x)
#' ncol(x)
#' nbasis(x)
#' 
#' # dimensions of a multiple NMF fit
#' x <- nmf(V, 3, nrun=3)
#' dim(x)
#' nrow(x)
#' ncol(x)
#' nbasis(x)
#' 
#' # dimensions of a multiple NMF fit (when keeping of the fits)
#' x <- nmf(V, 3, nrun=3, .opt='k')
#' dim(x)
#' nrow(x)
#' ncol(x)
#' nbasis(x)
#' 
#' 
#' 
NULL





#' Golub ExpressionSet from Brunet et al. Paper
#' 
#' The original data is related to Golub et al., and this version is the one
#' used and referenced in Brunet et al.  The samples are from 27 patients with
#' acute lymphoblastic leukemia (ALL) and 11 patients with acute myeloid
#' leukemia (AML).
#' 
#' The samples were assayed using Affymetrix Hgu6800 chips and the original
#' data on the expression of 7129 genes (Affymetrix probes) are available on
#' the Broad Institute web site (see references below).
#' 
#' The data in \code{esGolub} were obtained from the web site related to Brunet
#' et al.'s publication on an application of Nonnegative Matrix Factorization
#' (see link in section \emph{Source}).
#' 
#' They contain the 5,000 most highly varying genes according to their
#' coefficient of variation, and were installed in an object of class
#' \code{\link[Biobase]{ExpressionSet-class}}.
#' 
#' 
#' @name esGolub
#' @docType data
#' @format There are 3 covariates listed.
#' 
#' \itemize{
#' 
#' \item Samples: The original sample labels.  \item ALL.AML: Whether the
#' patient had AML or ALL. It is a \code{\link{factor}} with levels
#' \code{c('ALL', 'AML')}.  \item Cell: ALL arises from two different types of
#' lymphocytes (T-cell and B-cell).  This specifies which for the ALL patients;
#' There is no such information for the AML samples. It is a
#' \code{\link{factor}} with levels \code{c('T-cell', 'B-cell', NA)}.
#' 
#' }
#' 

#' @references
#' 
#' Brunet, J.P., Tamayo, P., Golub, T.R., and Mesirov, J.P. (2004)
#' \emph{Metagenes and molecular pattern discovery using matrix factorization}.
#' Proc Natl Acad Sci USA 101(12), 4164--4169.
#' 
#' T. R. Golub et al. (1999) \emph{Molecular Classification of Cancer: Class
#' Discovery and Class Prediction by Gene Expression Monitoring}.  Science,
#' 531-537, 1999
#' 
#' Original data from Golub et al.:\cr
#' \url{http://www-genome.wi.mit.edu/mpr/data_set_ALL_AML.html}
#' 

#' @source \url{http://www.broadinstitute.org/publications/broad872}
#' @keywords datasets
#' @examples
#' 
#'  	data(esGolub)
#'  	esGolub
#'  	\dontrun{pData(esGolub)}
#'  
NULL





#' Fast Combinatorial Non-Negative Least-Square
#' 
#' This function solves the following non-negative least square linear problem
#' using normal equations and the fast combinatorial strategy from Benthem and
#' Keenan (2004):
#' 
#' \deqn{% }{min ||Y - X K||_F, s.t. K>=0}\deqn{ \begin{array}{l}% }{min ||Y -
#' X K||_F, s.t. K>=0}\deqn{ \min \|Y - X K\|_F\\% }{min ||Y - X K||_F, s.t.
#' K>=0}\deqn{ \mbox{s.t. } K>=0% }{min ||Y - X K||_F, s.t. K>=0}\deqn{
#' \end{array}% }{min ||Y - X K||_F, s.t. K>=0}\deqn{ }{min ||Y - X K||_F, s.t.
#' K>=0} where \eqn{\|.\|_F} is the Frobenius norm.
#' 
#' The resulting algorithm is very fast to converge compared to other
#' approaches.
#' 
#' Within the \code{NMF} package, this algorithm is used internally by the
#' SNMF/R(L) algorithm from Kim and Park (2007) to solve general Nonnegative
#' Matrix Factorization (NMF) problems, using alternating non-negative
#' constrained least-squares. That is by iteratively and alternatively estimate
#' each matrix factor (see section \emph{References}).
#' 
#' It is provided separately so that it can be used to solve other types of
#' non-negative least squares problem. For faster computation, please the
#' internal -- non-exported -- function \code{NMF:::.fcnnls} The code is a port
#' from the original MATLAB code used in Kim and Park (2007) (see references).
#' 
#' Given two real matrices \eqn{Y} and \eqn{X}, of dimension \eqn{n \times p}{n
#' x p} and \eqn{n \times r}{n x r} respectively, this algorithm solves for the
#' optimal nonnegative matrix \eqn{K} (\eqn{r \times p}{r x p}) such that:
#' \deqn{% }{min ||Y - X K||_F, s.t. K>=0}\deqn{ \begin{array}{l}% }{min ||Y -
#' X K||_F, s.t. K>=0}\deqn{ \min \|Y - X K\|_F\\% }{min ||Y - X K||_F, s.t.
#' K>=0}\deqn{ \mbox{s.t. } K>=0% }{min ||Y - X K||_F, s.t. K>=0}\deqn{
#' \end{array}% }{min ||Y - X K||_F, s.t. K>=0}\deqn{ }{min ||Y - X K||_F, s.t.
#' K>=0} where \eqn{\|.\|_F} is the Frobenius norm.
#' 
#' It is based on the active/passive set method. It uses the unconstrained
#' solution \eqn{K_u} obtained from the unconstrained least squares problem,
#' i.e. \eqn{\min \|Y - X K\|_F^2}{min ||Y - X K||_F^2} , so as to determine
#' the initial passive sets.
#' 
#' @name fcnnls
#' @aliases fcnnls fcnnls,matrix,matrix-method fcnnls,numeric,matrix-method
#' fcnnls,ANY,numeric-method
#' @docType methods
#' @param x the coefficient matrix
#' @param y the target matrix to be approximated by \eqn{X K}.
#' @param verbose toggle verbosity (default is \code{FALSE}).
#' @param pseudo By default (\code{pseudo=FALSE}) the algorithm uses Gaussian
#' elimination to solve the successive internal linear problems, using the
#' \code{\link{solve}} function.  If \code{pseudo=TRUE} the algorithm uses
#' Moore-Penrose generalized \code{\link[corpcor]{pseudoinverse}} from the
#' \code{corpcor} package instead of \link{solve}.
#' @param ...  extra arguments passed to the internal function \code{.fcnnls}.
#' Currently not used.
#' @return The returned value is a list containing the following components:
#' \item{x}{ the estimated optimal matrix \eqn{K}.} \item{fitted}{ the fitted
#' matrix \eqn{X K}.} \item{residuals}{ the residual matrix \eqn{Y - X K}.}
#' \item{deviance}{ the residual sum of squares between the fitted matrix
#' \eqn{X K} and the target matrix \eqn{Y}. That is the sum of the square
#' residuals.} \item{passive}{ a \eqn{r x p} logical matrix containing the
#' passive set, that is the set of entries in \eqn{K} that are not null (i.e.
#' strictly positive).} \item{pseudo}{ a logical that is \code{TRUE} if the
#' computation was performed using the pseudoinverse. See argument
#' \code{pseudo}.}
#' @author Renaud Gaujoux
#' @seealso \code{\link{nmf}}
#' @references M. H. van Benthem and M. R. Keenan (2004).  \emph{Fast algorithm
#' for the solution of large-scale non-negativity-constrained least squares
#' problems}.  J. Chemometrics 2004, \bold{18}:441-450.
#' 
#' Kim, H. and Park, H. (2007).  \emph{Sparse non-negative matrix
#' factorizations via alternating non-negativity-constrained least squares for
#' microarray data analysis}.  Bioinformatics 2007; \bold{23(12)}:1495-502.
#' 
#' Original MATLAB code from Van Benthem and Keenan, slightly modified by H.
#' Kim:\cr \url{http://www.cc.gatech.edu/~hpark/software/fcnnls.m}
#' @keywords optimize multivariate regression
#' @examples
#' 
#' ## Define a random non-negative matrix matrix
#' n <- 200; p <- 20; r <- 3
#' V <- rmatrix(n, p)
#' 
#' ## Compute the optimal matrix K for a given X matrix
#' X <- rmatrix(n, r)
#' res <- fcnnls(X, V)
#' 
#' ## Compute the same thing using the Moore-Penrose generalized pseudoinverse
#' res <- fcnnls(X, V, pseudo=TRUE)
#' 
#' ## It also works in the case of single vectors
#' y <- runif(n)
#' res <- fcnnls(X, y)
#' # or
#' res <- fcnnls(X[,1], y)
#' 
#' 
NULL





#' Extract the Best Model and Result Object From NMF Runs
#' 
#' Methods \code{fit} and \code{minfit} extract the best NMF model and fit
#' objects respectively, from the result of function \code{\link{nmf}}. That is
#' the run that achieves the best approximation error across all runs performed
#' in the call to \code{nmf}.
#' 

#' 

#' 
#' Method \code{fit} returns the NMF model itself as an object of class
#' \code{\linkS4class{NMF}}, while method \code{minfit} returns the result of
#' the best run of function \code{nmf}. It is an object of class
#' \code{\linkS4class{NMFfit}} that, besides the NMF model, contains data about
#' the run (runtime, RNG seed, etc...).
#' 
#' They define a single interface to access the NMF model or fit result for
#' objects of class \code{NMFfit}, \code{NMFfitX1} and \code{NMFfitXn},
#' obtained from single runs, non-conservative multiple runs (i.e. that only
#' return the best fit), and conservative multiple runs (i.e. that return the
#' list of all the fits).
#' 
#' IMPORTANT: note that the behaviour of \code{fit,NMFfitX1} and
#' \code{fit,NMFfitXn} changed in version 0.5.3, as these functions returned
#' the best fit as an \code{NMFfit} object, which is what function
#' \code{minfit} now does.
#' 
#' @name fit-methods: Extracting NMF models
#' @aliases fit fit-methods fit,NMFfit-method fit,NMFfitX-method minfit
#' minfit-methods minfit,NMFfit-method minfit,NMFfitX-method
#' @docType methods
#' @param object an object that inherits from class \code{\linkS4class{NMFfit}}
#' or \code{\linkS4class{NMFfitX}}
#' @param ...  extra arguments to allow future extension.
#' @return
#' 
#' \itemize{
#' 
#' \item \code{fit} An object of class \code{\linkS4class{NMF}} that holds the
#' NMF model that achieves the best approximation error.
#' 
#' \item \code{minfit}: An object of class \code{\linkS4class{NMFfit}} that
#' contains data about the best run in addition to the actual NMF model.
#' 
#' }
#' @seealso \code{\link{residuals}}, \code{\link{nmf}}
#' @keywords methods
#' @examples
#' 
#' 
#' # generate a random target matrix
#' V <- rmatrix(100, 20)
#' 
#' # fit a single NMF model
#' res <- nmf(V, 3)
#' 
#' # extract best NMF model
#' fit(res)
#' # extract best run 
#' minfit(res)
#' # in the case of single runs it is the result itself
#' identical(minfit(res), res)
#' 
#' # perform non-conservative multiple runs
#' res <- nmf(V, 3, nrun=3)
#' # extract best NMF model
#' fit(res)
#' # extract best run 
#' minfit(res)
#' 
#' # perform conservative multiple runs
#' res <- nmf(V, 3, nrun=3, .opt='k')
#' # extract best NMF model
#' fit(res)
#' # extract best run 
#' minfit(res)
#' 
#' 
#' 
NULL





#' Plots a Heatmap of the Basis Vector and Mixture Coefficient Matrices from an
#' NMF model
#' 
#' Produces a heatmap of the basis vectors or mixture coefficients using a
#' \code{\link{aheatmap}}, with parameters tuned for displaying NMF results.
#' 
#' Both \code{basismap} and \code{coefmap} are S4 generic functions.
#' 
#' @name Plots: heatmaps
#' @aliases heatmap-NMF basismap basismap-methods basismap,NMF-method
#' basismap,NMFfitX-method coefmap coefmap-methods coefmap,NMF-method
#' coefmap,NMFfitX-method consensusmap consensusmap-methods
#' consensusmap,NMF-method consensusmap,NMFfitX-method
#' consensusmap,matrix-method consensusmap,NMF.rank-method
#' consensusmap,list-method
#' @docType methods
#' @param object An object that inherits from class \code{\linkS4class{NMF}}
#' (which includes class \code{\linkS4class{NMFfit}}) or class
#' \code{\linkS4class{NMFfitX}}.
#' @param tracks annotation track to add to the plot: \itemize{ \item'basis':
#' annotates the dominant basis component for each sample.  \item'consensus':
#' annotates the consensus clusters obtained by hierarchical clustering using
#' the consensus matrix as a similarity matrix and average linkage.  }
#' @param info if \code{TRUE} then information about the NMF fit is added to
#' the plot.  If not \code{FALSE}, the this argument is passed to the heatmap
#' plotting function \code{\link{aheatmap}}.
#' @param subsetRow Argument that specifies how to filter the features that
#' will appear in the heatmap.  When \code{FALSE} (default), all the features
#' are used.  Besides the values supported by argument \code{subsetRow} of
#' \code{\link{aheatmap}}, other possible values are:
#' 
#' \itemize{ \item \code{TRUE}: only the features that are basis-specific are
#' used.  The default selection method is from Kim and Park (2007). Other
#' selection methods can be specified as a single \code{character} string or a
#' custom \code{function} (cf. argument \code{method} for function
#' \code{\link{extractFeatures}}).
#' 
#' \item a single \code{character} string that specifies the filtering method
#' to be used to select the basis-specific features that should appear in the
#' heatmap (cf. argument \code{method} for function
#' \code{\link{extractFeatures}}).
#' 
#' \item a single \code{numeric} value: in this case the features are sorted in
#' decreasing order according to their contribution to each basis, and only the
#' top \code{filter} of each basis are kept.  Multiple occurrences are removed.
#' }
#' 

#' @param ...  Used to pass extra parameters to the subsequent call to the
#' internal drawing function \code{heatmap.2.plus}.
#' @section Methods: \describe{
#' 
#' \item{basismap}{\code{signature(object = "NMF")}: plots a heatmap of the
#' basis vector matrix. The features used in the heatmap can be filtered in
#' various different way using the argument \code{filter}.}
#' 
#' \item{coefmap}{\code{signature(object = "NMF")}: plots a heatmap of the
#' mixture coefficient matrix.}
#' 
#' \item{basismap, coefmap}{\code{signature(object = "NMFfitX")}: These methods
#' are used when \code{object} is the result of a multirun-NMF fit, as returned
#' by the function \code{\link{nmf}} when called with argument \code{nrun>1}.
#' It simply calls the corresponding method on the best NMF fit of
#' \code{object}: \code{basismap(object, ...)} is roughly equivalent to
#' \code{basismap(fit(object), ...)} and \code{coefmap(object, ...)} to
#' \code{coefmap(fit(object), ...)}. } }
#' @author Renaud Gaujoux
#' @seealso \code{\link{basis}}, \code{\link[=coef]{coef,NMF-method}},
#' \code{\link{consensusmap}}
#' @keywords methods plot
#' @examples
#' 
#' 
#' # generate a random target matrix
#' V <- rmatrix(100, 20)
#' res <- nmf(V, 3)
#' # Generates heatmaps
#' \dontrun{
#' basismap(res)
#' coefmap(res)
#' }
#' 
#' # filtered in different ways
#' \dontrun{
#' basismap(res, filter=TRUE)
#' basismap(res, filter=5)
#' basismap(res, filter='max')
#' basismap(res, filter=c(1,2,3, 10, 15, 20))
#' }
#' # the same not re-ordering the features
#' \dontrun{ basismap(res, filter=c(1,2,3, 10, 15, 20), Rowv=FALSE)}
#' 
#' ## multirun NMF
#' res <- nmf(V, 3, nrun=3)
#' \dontrun{consensusmap(res)}
#' # compare the heatmap when there is a true signal
#' Vs <- syntheticNMF(100, 3, 30, noise=TRUE)
#' res <- nmf(Vs, 3, nrun=10)
#' \dontrun{consensusmap(res)} 
#' 
#' 
#' 
NULL





#' Interface Class for Nonnegative Matrix Factorisation Models
#' 
#' This is a \emph{virtual class} that defines a common interface to handle
#' Nonnegative Matrix Factorisation models (NMF models) in a generic way.
#' 
#' It provides the definition for a minimum set of generic methods that are
#' used in common computations and tasks in the context of Nonnegative Matrix
#' Factorisations.
#' 
#' Class \code{NMF} makes it easy to develop new models that integrates well
#' into the general framework implemented by the \emph{NMF} package.
#' 
#' Following a few simple guidelines, new models benefit from all the
#' functionalities available to built-in NMF models -- that derive themselves
#' from class \code{NMF}.  See section \emph{Defining new NMF models} below.
#' 
#' See section \code{\linkS4class{NMFstd}}, references and links therein for
#' details on the standard NMF model and its -- built-in -- extensions.
#' 
#' 
#' @name NMF-class: Common interface class for NMF models
#' @aliases NMF-class consensus,NMF-method distance distance,matrix,NMF-method
#' distance,missing,missing-method fitted,NMF-method is.empty.nmf
#' is.empty.nmf,NMF-method modelname modelname,NMF-method show,NMF-method
#' summary,NMF-method
#' @docType class
#' @section Slots: This class contains a single slot, that is used internally
#' during the computations.
#' 
#' \describe{ \item{list("misc")}{A list that is used internally to temporarily
#' store algorithm parameters during the computation.}\item{:}{A list that is
#' used internally to temporarily store algorithm parameters during the
#' computation.} }
#' 
#' The purpose of this class is to define a common interface for NMF models as
#' a collection of generic methods.  Classes that inherits from class
#' \code{NMF} are responsible for the management of data storage and the
#' implementation of the interface's pure virtual methods.
#' @author Renaud Gaujoux
#' @seealso Main interface to perform NMF in \code{\link{nmf-methods}}.
#' 
#' Built-in NMF models and factory method in \code{\link{nmfModel}}.
#' 
#' Method \code{\link{seed}} to set NMF objects with values suitable to start
#' algorithms with.
#' @references
#' 
#' Definition of Nonnegative Matrix Factorization in its modern formulation:
#' 
#' Lee D.D. and Seung H.S. (1999).  Learning the parts of objects by
#' non-negative matrix factorization.  \emph{Nature}, \bold{401}, 788--791.
#' 
#' Historical first definition and algorithms:
#' 
#' Paatero, P., Tapper, U. (1994).  Positive matrix factorization: A
#' non-negative factor model with optimal utilization of error estimates of
#' data values.  \emph{Environmetrics}, \bold{2}, 111--126 ,
#' doi:10.1002/env.3170050203.
#' 
#' Reference for some utility functions:
#' 
#' Kim, H. and Park, H. (2007).  Sparse non-negative matrix factorizations via
#' alternating non-negativity-constrained least squares for microarray data
#' analysis.  \emph{Bioinformatics}.
#' 
#' Hoyer (2004).  Non-negative matrix factorization with sparseness
#' constraints.  \emph{Journal of Machine Learning Research}, \bold{5},
#' 1457-1469.
#' @keywords classes
#' @examples
#' 
#' 
#' # show all the NMF models available (i.e. the classes that inherit from class NMF)
#' nmfModels()
#' # show all the built-in NMF models available
#' nmfModels(builtin.only=TRUE)
#' 
#' # class NMF is a virtual class so cannot be instantiated: 
#' # the following generates an error
#' \dontrun{new('NMF')}
#' 
#' # To instantiate an NMF model, use the factory method nmfModel. see ?nmfModel
#' nmfModel()
#' nmfModel(3)
#' nmfModel(3, model='NMFns')
#' 
#' 
#' 
NULL





#' Comparing Results from Different NMF Runs
#' 

#' 
#' This functions allow to compare the results of different NMF runs. The
#' results do not need to be from the same algorithm, nor even of the same
#' dimension.
#' 

#' 
#' \describe{
#' 
#' \item{as.NMFList}{: wrap the arguments into a \code{\linkS4class{NMFList}}
#' object.}
#' 
#' \item{compare}{: shortcut for \code{summary(as.NMFList(object)} (cf.
#' \code{summary} method below)}.
#' 
#' \item{plot}{: plot on a single graph the residuals tracks for each element
#' in \code{x}.  See function \code{\link{nmf}} for details on how to enable
#' the tracking of residuals.}
#' 
#' \item{runtime}{: returns the computational time used to compute all the
#' results in the list, as stored in slot \code{runtime} of \code{object}.
#' 
#' The time is computed using the function \code{\link{system.time}} which
#' returns object of class \code{\link[=proc.time]{proc_time}}.
#' 
#' Note that argument \code{...} is not used.  }
#' 
#' \item{summary}{: \code{summary} method for objects of class \code{NMFList}.
#' 
#' It compute summary measures for each NMF result in the list and return them
#' in rows in a \code{data.frame}.  By default all the measures are included in
#' the result, and \code{NA} values are used where no data is available or the
#' measure does not apply to the result object (e.g. the dispersion for single
#' NMF runs is not meaningful).  This method is very useful to compare and
#' evaluate the performance of different algorithms.  }
#' 
#' }
#' 
#' @name Comparing the results of different NMF runs
#' @aliases nmf-compare as.NMFList as.NMFList-methods as.NMFList,ANY-method
#' compare compare-methods compare,list-method plot,NMFList,missing-method
#' plot,NMFList-method summary,NMFList-method
#' @docType methods
#' @param object A \code{list} or an object of class
#' \code{\linkS4class{NMFList}}.
#' @param select the columns to be output in the result \code{data.frame}.  The
#' column are given by their names (partially matched).  The column names are
#' the names of the summary measures returned by the \code{summary} methods of
#' the corresponding NMF results.
#' @param sort.by the sorting criteria, i.e. a partial match of a column name,
#' by which the result \code{data.frame} is sorted.  The sorting direction
#' (increasing or decreasing) is computed internally depending on the chosen
#' criteria (e.g. decreasing for the cophenetic coefficient, increasing for the
#' residuals).
#' @param unlist boolean to specify if the arguments should be unlisted before
#' wrapping them into a \code{NMFList} object or comparing them.
#' @param x An object of class \code{\linkS4class{NMFList}}.
#' @param ...  Used to pass extra arguments to subsequent calls: \itemize{
#' \item in \code{as.NMFList} the list of NMF results to wrap into a
#' \code{NMFList} object.  \item in \code{plot}: graphical parameters passed to
#' the \code{plot} function.  \item in \code{compare} and \code{summary}: extra
#' arguments passed to the \code{summary} method of each result object (cf.
#' \code{\link{summary,NMF-method}}).  }
#' @author Renaud Gaujoux
#' @seealso \linkS4class{NMFfitX1}, \linkS4class{NMFfitXn},
#' \link[=summary,NMF-method]{summary}
#' @references
#' 
#' Brunet, J.P. et al. (2004) \emph{Metagenes and molecular pattern discovery
#' using matrix factorization} Proc Natl Acad Sci USA 101(12), 4164--4169.
#' 
#' Kim, H. and Park, H. (2007).  \emph{Sparse non-negative matrix
#' factorizations via alternating non-negativity-constrained least squares for
#' microarray data analysis}.  Bioinformatics 2007; \bold{23(12)}:1495-502.
#' \url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
#' @examples
#' 
#' 
#' 
#' # generate a synthetic dataset with known classes: 50 features, 18 samples (5+5+8)
#' n <- 50; counts <- c(5, 5, 8);
#' V <- syntheticNMF(n, counts, noise=TRUE)
#' \dontrun{aheatmap(V)}
#' 
#' # build the class factor
#' groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))
#' 
#' # perform multiple runs of NMF (keep best only)
#' res <- nmf(V, 3, nrun=5)
#' res
#' 
#' # compute summary measures
#' summary(res)
#' 
#' # compute more summary measures
#' summary(res, target=V, class=groups)
#' 
#' # plot a heatmap of the consensus matrix with extra annotations
#' \dontrun{consensusmap(res, annCol=groups)}
#' 
#' # retrieve the predicted clusters of samples
#' predict(res)
#' 
#' # perform multiple runs of NMF and keep all the runs
#' res <- nmf(V, 3, nrun=5, .options='k')
#' res
#' 
#' # extract best fit
#' minfit(res)
#' # extract best NMF model
#' fit(res)
#' 
#' # compute/show computational times
#' runtime.all(res)
#' seqtime(res)
#' 
#' 
#' 
NULL





#' Defunct Functions in NMF Package
#' 

#' 
#' The functions or variables listed here are no longer part of package NMF as
#' they are not needed (any more).
#' 

#' 
#' \code{extra} was used to access slot \code{extra} in class \code{NMFfit}.
#' It has been replaced in version 0.3 by methods \code{$} and \code{$<-}.  See
#' \code{\link[=$,NMFfit-method]{$}}.
#' 
#' \code{clusters} is a defunct synonym for \code{\link{predict}} (Defunct in
#' 0.6).
#' 
#' \code{metaHeatmap,matrix} is a defunct synonym for \code{\link{aheatmap}}
#' that provides an enhanced heatmap drawing function.
#' 
#' @aliases defunct-NMF extra extra-methods extra,NMFfit-method clusters
#' clusters-methods clusters,NMF-method clusters,matrix-method
#' metaHeatmap,matrix-method has.track
#' @keywords internal
NULL





#' Deprecated Functions in NMF Package
#' 
#' These functions are provided for compatibility with older versions of NMF
#' package only, and may be defunct as soon as the next release.
#' 

#' 
#' \code{errorPlot} is a deprecated synonym for
#' \code{\link{plot,NMFfit-method}}.  \code{newNMF} is a deprecated synonym for
#' \code{\link{nmfModel}}.  \code{metaHeatmap} is a deprecated synonym for
#' \code{\link{basismap}}, \code{\link{coefmap}} or \code{\link{consensusmap}}
#' , to plot a heatmap of the basis matrix, mixture coefficient matrix or
#' consensus matrix respectively.
#' 
#' @name Deprecated functions
#' @aliases deprecated-NMF errorPlot errorPlot,NMFfit-method
#' errorPlot,ANY-method newNMF metaHeatmap metaHeatmap-methods
#' metaHeatmap,NMF-method metaHeatmap,NMF-method metaHeatmap,NMFfitX-method
#' @docType methods
#' @param x An object that inherits from class \code{\linkS4class{NMF}} or
#' \code{\linkS4class{NMFfit}}.
#' @param ...  extra graphical arguments passed to the \code{plot} function.
#' @keywords internal
NULL





#' Methods to Compare Two NMF Models
#' 
#' The method \code{nmf.equal} is used to compare two NMF models. The
#' comparison is performed entry-wise on both the basis and the mixture
#' coefficient matrices.
#' 
#' The two compared objects can be of class \code{\linkS4class{NMF}},
#' \code{\linkS4class{NMFfit}}, \code{\linkS4class{NMFfitX}} or \code{list},
#' allowing the direct comparison of models, fits or list of fits.  In the case
#' of \code{\linkS4class{NMFfitXn}} or \code{list} objects, either only the
#' best fits (default) or each fit separately are compared.
#' 

#' 
#' The comparison is performed using either the function
#' \code{\link{identical}} or the function \code{\link{all.equal}} on the
#' \code{\linkS4class{NMF}} objects stored in \code{x} and \code{y}.
#' 
#' Currently, comparing lists of different length will return \code{FALSE}.
#' 
#' @name nmf.equal-methods
#' @aliases nmf.equal nmf.equal-methods nmf.equal,NMFfit,NMF-method
#' nmf.equal,NMFfit,NMFfit-method nmf.equal,NMFfitX,NMF-method
#' nmf.equal,NMFfitX1,NMFfitX1-method nmf.equal,list,list-method
#' nmf.equal,NMF,NMF-method nmf.equal,NMF,NMFfit-method
#' nmf.equal,NMF,NMFfitX-method
#' @docType methods
#' @param x an object that inherits from class \code{\linkS4class{NMF}},
#' \code{\linkS4class{NMFfit}}, \code{\linkS4class{NMFfitX}} or a list of
#' these.
#' @param y an object that inherits from class \code{\linkS4class{NMF}},
#' \code{\linkS4class{NMFfit}}, \code{\linkS4class{NMFfitX}} or a list of
#' these.
#' @param identical a single \code{logical} (default to \code{TRUE}) that
#' specifies if the comparison should be made using the function
#' \code{\link{identical}} (i.e. exact equality) or the function
#' \code{\link{all.equal}} that allows for some numerical tolerance (controlled
#' by argument \code{tolerance}).
#' @param all If \code{all=FALSE} (default) the comparison is done between the
#' best fits of each list, i.e. the NMF model with the lowest deviance.
#' Otherwise, the comparison is made pairwise.
#' @param vector a logical used when \code{all=TRUE} that specifies if all the
#' pairwise comparisons must be returned as a vector. If \code{FALSE}
#' (default), the function will return \code{TRUE} only if all the fits in the
#' lists are identical/equals.
#' @param ...  extra arguments that are passed to \code{\link{all.equal}} (e.g.
#' \code{tolerance}) when \code{identical=FALSE}, or to the internal call to
#' \code{nmf.equal} when doing pairwise comparison.
#' @return
#' 
#' When \code{identical=TRUE} this function returns \code{TRUE} if the two
#' \code{NMF} models are identical, and \code{FALSE} otherwise.
#' 
#' When \code{identical=FALSE} this function returns the result of
#' \code{all.equal} applied to the two \code{NMF} models.  That is it returns
#' either \code{TRUE} if the objects are considered equals (given the numerical
#' tolerance), or a \code{character} vector giving a summary of the observed
#' differences.  The numerical tolerance and comparison parameters can be
#' passed through \code{...}.
#' 
#' When doing pairwise comparison (\code{all=TRUE}) and \code{vector=FALSE},
#' the function returns \code{TRUE} only if all the fits in the lists are
#' identical/equals.  When \code{vector=TRUE} the returned value is a logical
#' vector of the same length as the input lists.
#' @seealso \code{\link{identical}}, \code{\link{all.equal}}
#' @keywords methods
#' @examples
#' 
#' 
#' # create a random NMF model
#' a <- rnmf(100, 3, 20)
#' 
#' # compare a with itself
#' nmf.equal(a, a)
#' nmf.equal(a, a, identical=FALSE)
#' \dontshow{
#' stopifnot( nmf.equal(a, a) )
#' stopifnot( nmf.equal(a, a, identical=FALSE) )
#' }
#' 
#' # compare with a slightly different model
#' b <- a
#' basis(b) <- basis(b) + .Machine$double.eps^0.6
#' nmf.equal(a, b)
#' \dontshow{ stopifnot( !nmf.equal(a, b) ) }
#' nmf.equal(a, b, identical=FALSE) 
#' \dontshow{ stopifnot( nmf.equal(a, b, identical=FALSE) ) }
#' 
#' # compare with a model that is more different
#' b <- a
#' basis(b) <- basis(b) + .Machine$double.eps^0.4
#' nmf.equal(a, b)
#' \dontshow{ stopifnot( !nmf.equal(a, b) ) }
#' nmf.equal(a, b, identical=FALSE)
#' \dontshow{ stopifnot( !isTRUE(nmf.equal(a, b, identical=FALSE)) ) }
#' 
#' # the same works on results from function 'nmf'
#' V <- rmatrix(100, 20)
#' a <- nmf(V, 3)
#' 
#' # compare the result NMFfit object with the fitted NMF object
#' nmf.equal(a, fit(a))
#' 
#' # comparing with another run is unlikely to give the same result
#' b <- nmf(V, 3)
#' nmf.equal(a, b)
#' \dontshow{ stopifnot( !nmf.equal(a, b) ) }
#' # even when allowing for non identical objects
#' nmf.equal(a, b, FALSE)
#' \dontshow{ stopifnot( !isTRUE(nmf.equal(a, b, FALSE)) ) }
#' # but using the same random seed returns identical results
#' b <- nmf(V, 3, seed=getRNG(a))
#' nmf.equal(a, b)
#' \dontshow{ stopifnot( nmf.equal(a, b) ) }
#' 
#' # compare for multiple runs
#' a <- nmf(V, 3, nrun=3)
#' 
#' # nmf.equal compares the best fit
#' nmf.equal(a, fit(a))
#' \dontshow{ stopifnot( nmf.equal(a, fit(a)) ) }
#' # compare with another (random) NMF model
#' nmf.equal(a, rnmf(a))
#' # running another set of fits is unlikely to give an identical fit
#' b <- nmf(V, 3, nrun=3)
#' nmf.equal(a, b)
#' 
#' # but again: using the same random seed returns an identical fit
#' b <- nmf(V, 3, nrun=3, seed=getRNG1(a))
#' nmf.equal(a, b)
#' \dontshow{ stopifnot( nmf.equal(a, b) ) }
#' 
#' 
#' # compare all the runs from multiple runs
#' a <- nmf(V, 3, nrun=3, .opt='k')
#' 
#' # the comparison is done between the best fit
#' nmf.equal(a, fit(a))
#' \dontshow{ stopifnot( nmf.equal(a, fit(a)) ) }
#' # running another set of fits is unlikely to give identical results (for any of the runs)
#' b <- nmf(V, 3, nrun=3, .opt='k')
#' nmf.equal(a, b)
#' nmf.equal(a, b, all=TRUE)
#' nmf.equal(a, b, all=TRUE, vector=TRUE)
#' 
#' # but again: using the same random seed returns identical results (for each run)
#' b <- nmf(V, 3, nrun=3, seed=getRNG1(a), .opt='k')
#' nmf.equal(a, b)
#' nmf.equal(a, b, all=TRUE)
#' nmf.equal(a, b, all=TRUE, vector=TRUE)
#' \dontshow{ 
#' stopifnot( nmf.equal(a, b) )
#' stopifnot( nmf.equal(a, b, all=TRUE) )
#' stopifnot( all(nmf.equal(a, b, all=TRUE, vector=TRUE)) )
#' }
#' 
#' 
NULL





#' Estimate optimal rank for Nonnegative Matrix Factorization (NMF) models
#' 
#' A critical parameter in NMF algorithms is the factorization rank \eqn{r}.
#' It defines the number of basis effects used to approximate the target
#' matrix. Function \code{nmfEstimateRank} helps in choosing an optimal rank by
#' implementing simple approaches proposed in the literature.
#' 
#' Given a NMF algorithm and the target matrix, a common way of estimating
#' \eqn{r} is to try different values, compute some quality measures of the
#' results, and choose the best value according to this quality criteria. See
#' \emph{Brunet et al. (2004)} and \emph{Hutchins et al. (2008)}.
#' 
#' The function \code{nmfEstimateRank} allows to launch this estimation
#' procedure. It performs multiple NMF runs for a range of rank of
#' factorization and, for each, returns a set of quality measures together with
#' the associated consensus matrix.
#' 
#' @aliases nmfEstimateRank plot.NMF.rank
#' @param method A single NMF algorithm, in one of the format accepted by
#' interface \code{\link{nmf}}.
#' @param na.rm single logical that specifies if the rank for which the
#' measures are NA values should be removed from the graph or not (default to
#' \code{FALSE}).  This is useful when plotting results which include NAs due
#' to error during the estimation process. See argument \code{stop} for
#' \code{nmfEstimateRank} below.
#' @param nrun a \code{numeric} giving the number of run to perform for each
#' value in \code{range}.
#' @param range a \code{numeric} vector containing the ranks of factorization
#' to try.
#' @param ref reference object of class \code{NMF.rank}, as returned by
#' function \code{nmfEstimateRank}.  The measures contained in \code{ref} are
#' used and plotted as a reference.  The associated curves are drawn in
#' \emph{red}, while those from \code{x} are drawn in \emph{blue}.
#' @param verbose toggle verbosity.  This parameter only affects the verbosity
#' of the outer loop over the values in \code{rank}. To print verbose (resp.
#' debug) messages from each NMF run, one can use \code{.options='v'} (resp.
#' \code{.options='d'}) that will be passed to the \code{\link{nmf}} method.
#' @param stop logical flag for running the estimation process with fault
#' tolerance.  When \code{TRUE}, the whole execution will stop if any error is
#' raised.  When \code{FALSE} (default), the runs that raise an error will be
#' skipped, and the execution will carry on. The summary measures for the runs
#' with errors are set to NA values, and a warning is thrown.
#' @param what a \code{character} string that partially matches one of the
#' following item: \code{'all'}, \code{'cophenetic'}, \code{'rss'},
#' \code{'residuals'} , \code{'dispersion'}.  It specifies which measure must
#' be plotted (\code{what='all'} plots all the measures).
#' @param x For \code{nmfEstimateRank} a target object to be estimated, in one
#' of the format accepted by interface \code{\link{nmf}}.
#' 
#' For \code{plot.NMF.rank} an object of class \code{NMF.rank} as returned by
#' function \code{nmfEstimateRank}.
#' @param \dots For \code{nmfEstimateRank}, these are extra parameters passed
#' to interface \code{nmf}. Note that the same parameters are used for each
#' value of the rank.  See \code{\link{nmf}}.
#' 
#' For \code{plot.NMF.rank}, these are extra graphical parameter passed to the
#' standard function \code{plot}. See \code{\link{plot}}.
#' @return A S3 object (i.e. a list) of class \code{NMF.rank} with the
#' following slots: \item{measures }{a \code{data.frame} containing the quality
#' measures for each rank of factorizations in \code{range}. Each row
#' corresponds to a measure, each column to a rank. } \item{consensus }{ a
#' \code{list} of consensus matrices, indexed by the rank of factorization (as
#' a character string).}
#' 
#' \item{fit }{ a \code{list} of the fits, indexed by the rank of factorization
#' (as a character string).}
#' @author Renaud Gaujoux
#' @seealso nmf
#' @references
#' 
#' \emph{Metagenes and molecular pattern discovery using matrix factorization}
#' Brunet, J.~P., Tamayo, P., Golub, T.~R., and Mesirov, J.~P. (2004) Proc Natl
#' Acad Sci U S A 101(12), 4164--4169.
#' @examples
#' 
#' 
#' set.seed(123456)
#' n <- 50; r <- 3; m <- 20
#' V <- syntheticNMF(n, r, m, noise=TRUE)
#' 
#' # Use a seed that will be set before each first run
#' \dontrun{res.estimate <- nmfEstimateRank(V, seq(2,5), method='brunet', nrun=10, seed=123456)}
#' 
#' # plot all the measures
#' \dontrun{plot(res.estimate)}
#' # or only one: e.g. the cophenetic correlation coefficient
#' \dontrun{plot(res.estimate, 'cophenetic')}
#' 
#' 
NULL





#' Base Class for to store Nonnegative Matrix Factorisation results
#' 
#' Base class to handle the results of general \strong{Non-negative Matrix
#' Factorisation} algorithms (NMF).
#' 
#' It provides a general structure and generic functions to manage the results
#' of NMF algorithms.  It contains a slot with the fitted NMF model (see slot
#' \code{fit}) as well as data about the methods and parameters used to compute
#' the factorization.
#' 
#' The purpose of this class is to handle in a generic way the results of NMF
#' algorithms. Its slot \code{fit} contains the fitted NMF model as an object
#' of class \code{\linkS4class{NMF}}.
#' 
#' Other slots contains data about how the factorization has been computed,
#' such as the algorithm and seeding method, the computation time, the final
#' residuals, etc\dots{}
#' 
#' Class \code{NMFfit} acts as a wrapper class for its slot \code{fit}.  It
#' inherits from interface class \code{\linkS4class{NMF}} defined for generic
#' NMF models.  Therefore, all the methods defined by this interface can be
#' called directly on objects of class \code{NMFfit}. The calls are simply
#' dispatched on slot \code{fit}, i.e.  the results are the same as if calling
#' the methods directly on slot \code{fit}.
#' 
#' @name NMFfit-class
#' @aliases NMFfit-class algorithm algorithm<- algorithm,NMFfit-method
#' algorithm<-,NMFfit,ANY-method basis<-,NMFfit,matrix-method
#' basis,NMFfit-method coef,NMFfit-method coef<-,NMFfit,matrix-method
#' distance,matrix,NMFfit-method $ $<- $,NMFfit-method $<-,NMFfit-method fit<-
#' fit<-,NMFfit,NMF-method fitted,NMFfit-method initialize,NMFfit-method
#' modelname,NMFfit-method niter-methods niter niter<- niter,NMFfit-method
#' niter<-,NMFfit,numeric-method nrun,NMFfit-method objective objective<-
#' objective,NMFfit-method objective<-,NMFfit,character-method
#' objective<-,NMFfit,function-method residuals<- residuals<--methods
#' residuals<-,NMFfit-method runtime runtime,NMFfit-method
#' runtime.all,NMFfit-method seeding seeding-methods seeding,NMFfit-method
#' seeding<- seeding<-,NMFfit-method show,NMFfit-method summary,NMFfit-method
#' @docType class
#' @section Slots: \describe{ \item{list("fit")}{An object that inherits from
#' class \code{"NMF"}.  It contains the fitted NMF model.
#' 
#' Note that class \code{"NMF"} is a virtual class.  The default class for this
#' slot is \code{NMFstd}, that implements the standard NMF model.  }\item{:}{An
#' object that inherits from class \code{"NMF"}.  It contains the fitted NMF
#' model.
#' 
#' Note that class \code{"NMF"} is a virtual class.  The default class for this
#' slot is \code{NMFstd}, that implements the standard NMF model.  }
#' 
#' \item{list("residuals")}{A \code{"numeric"} vector that contains the final
#' residuals or the residuals track between the target matrix and its NMF
#' estimate(s).  Default value is \code{numeric()}.
#' 
#' See method \code{\link{residuals}} for details on accessor methods and main
#' interface \code{\link{nmf}} for details on how to compute NMF with residuals
#' tracking.}\item{:}{A \code{"numeric"} vector that contains the final
#' residuals or the residuals track between the target matrix and its NMF
#' estimate(s).  Default value is \code{numeric()}.
#' 
#' See method \code{\link{residuals}} for details on accessor methods and main
#' interface \code{\link{nmf}} for details on how to compute NMF with residuals
#' tracking.}
#' 
#' \item{list("method")}{A single \code{"character"} string that contains the
#' name of the algorithm used to fit the model. Default value is
#' \code{''}.}\item{:}{A single \code{"character"} string that contains the
#' name of the algorithm used to fit the model. Default value is \code{''}.}
#' 
#' \item{list("seed")}{A single \code{"character"} string that contains the
#' name of the seeding method used to seed the algorithm that computed the NMF.
#' Default value is \code{''}.  See \code{\link{nmf-methods}} for more details.
#' }\item{:}{A single \code{"character"} string that contains the name of the
#' seeding method used to seed the algorithm that computed the NMF. Default
#' value is \code{''}.  See \code{\link{nmf-methods}} for more details. }
#' 
#' \item{list("rng")}{An object that contains the RNG settings used for the
#' fit.  Currently the settings are stored as an integer vector, the value of
#' \code{\link{.Random.seed}} at the time the object is created.  It is
#' initialized by the \code{initialized} method.  See \code{\link{getRNG}} for
#' more details. }\item{:}{An object that contains the RNG settings used for
#' the fit.  Currently the settings are stored as an integer vector, the value
#' of \code{\link{.Random.seed}} at the time the object is created.  It is
#' initialized by the \code{initialized} method.  See \code{\link{getRNG}} for
#' more details. }
#' 
#' \item{list("distance")}{Either a single \code{"character"} string that
#' contains the name of the built-in objective function, or a \code{function}
#' that measures the residuals between the target matrix and its NMF estimate.
#' }\item{:}{Either a single \code{"character"} string that contains the name
#' of the built-in objective function, or a \code{function} that measures the
#' residuals between the target matrix and its NMF estimate. }
#' 
#' \item{list("parameters")}{A \code{"list"} that contains the extra parameters
#' specific to the algorithm used to fit the model. }\item{:}{A \code{"list"}
#' that contains the extra parameters specific to the algorithm used to fit the
#' model. }
#' 
#' \item{list("runtime")}{Object of class \code{"proc_time"} that contains
#' various measures of the time spent to fit the model.  See
#' \code{\link[base]{system.time}}}\item{:}{Object of class \code{"proc_time"}
#' that contains various measures of the time spent to fit the model.  See
#' \code{\link[base]{system.time}}}
#' 
#' \item{list("options")}{A \code{"list"} that contains the options used to
#' compute the object.}\item{:}{A \code{"list"} that contains the options used
#' to compute the object.}
#' 
#' \item{list("extra")}{A \code{"list"} that contains extra miscellaneous data
#' for internal usage only.  For example it can be used to store extra
#' parameters or temporary data, without the need to explicitly extend the
#' \code{NMFfit} class.  Currently built-in algorithms only use this slot to
#' store the number of iterations performed to fit the object.  Data that need
#' to be easily accessible by the end-user should be stored using the \code{$}
#' and \code{$<-} methods that access and set the \code{list} slot \code{misc}
#' inherited from class \code{NMF}.  }\item{:}{A \code{"list"} that contains
#' extra miscellaneous data for internal usage only.  For example it can be
#' used to store extra parameters or temporary data, without the need to
#' explicitly extend the \code{NMFfit} class.  Currently built-in algorithms
#' only use this slot to store the number of iterations performed to fit the
#' object.  Data that need to be easily accessible by the end-user should be
#' stored using the \code{$} and \code{$<-} methods that access and set the
#' \code{list} slot \code{misc} inherited from class \code{NMF}.  }
#' 
#' }
#' @author Renaud Gaujoux
#' @seealso Main interface to perform NMF in \code{\link{nmf-methods}}.
#' 
#' Method \code{\link{seed}} to set NMF objects with values suitable to start
#' algorithms with.
#' @keywords classes
#' @examples
#' 
#' 
#' # run default NMF algorithm on a random matrix
#' n <- 50; r <- 3; p <- 20
#' V <- rmatrix(n, p)  
#' res <- nmf(V, r)							
#' 
#' # result class is NMFfit
#' class(res)
#' 
#' # show result
#' res
#' 
#' # compute summary measures
#' summary(res)
#' 
#' 
NULL





#' Class to Store the Result from Multiple Runs of a NMF Algorithm when Only
#' the Best Fit is Kept
#' 
#' This class is used to return the result from a multiple run of a single NMF
#' algorithm performed with function \code{nmf} with the -- default -- option
#' \code{keep.all=FALSE} (cf. \code{\link{nmf}}).
#' 
#' It extends both classes \code{\linkS4class{NMFfitX}} and
#' \code{\linkS4class{NMFfit}}, and stores a the result of the best fit in its
#' \code{NMFfit} structure.
#' 
#' Beside the best fit, this class allows to hold data about the computation of
#' the multiple runs, such as the number of runs, the CPU time used to perform
#' all the runs, as well as the consensus matrix.
#' 
#' Due to the inheritance from class \code{NMFfit}, objects of class
#' \code{NMFfitX1} can be handled exactly as the results of single NMF run --
#' as if only the best run had been performed.
#' 
#' 
#' @name NMFfitX1-class
#' @aliases NMFfitX1-class consensus,NMFfitX1-method fit,NMFfitX1-method
#' minfit,NMFfitX1-method nrun,NMFfitX1-method runtime.all,NMFfitX1-method
#' show,NMFfitX1-method
#' @docType class
#' @section Slots: \describe{
#' 
#' \item{list("consensus")}{ object of class \code{"matrix"} used to store the
#' consensus matrix based on all the runs.}\item{:}{ object of class
#' \code{"matrix"} used to store the consensus matrix based on all the runs.}
#' 
#' \item{list("nrun")}{an \code{integer} that contains the number of runs
#' performed to compute the object.  }\item{:}{an \code{integer} that contains
#' the number of runs performed to compute the object.  }
#' 
#' \item{list("rng1")}{an object that contains RNG settings used for the first
#' run. See \code{\link{getRNG1}}.  }\item{:}{an object that contains RNG
#' settings used for the first run. See \code{\link{getRNG1}}.  }
#' 
#' \item{list("runtime.all")}{object of class \code{"proc_time"} that contains
#' various measures of the time spent to perform all the runs (inherited from
#' \code{NMFfitX})}\item{:}{object of class \code{"proc_time"} that contains
#' various measures of the time spent to perform all the runs (inherited from
#' \code{NMFfitX})} }
#' @author Renaud Gaujoux
#' @seealso \code{\linkS4class{NMFfitX}} , \code{\link{nmf-methods}},
#' \code{\link{nmf-multiple}}
#' @keywords classes
#' @examples
#' 
#' 
#' # generate a synthetic dataset with known classes
#' n <- 50; counts <- c(5, 5, 8);
#' V <- syntheticNMF(n, counts, noise=TRUE)
#' 
#' # build the class factor
#' groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))
#' 
#' # perform multiple runs of one algorithm, keeping only the best fit (default)
#' res <- nmf(V, 3, nrun=5) 
#' res
#' #NOTE: the implicit nmf options are .options=list(keep.all=FALSE) or .options='-k'
#' 
#' # compute summary measures
#' summary(res)
#' # get more info
#' summary(res, target=V, class=groups)
#' 
#' # show computational time
#' runtime.all(res)
#' 
#' 
NULL





#' Class to Store the Result from Multiple Runs of a NMF Algorithm when Keeping
#' All the Fits
#' 
#' This class is used to return the result from a multiple run of a single NMF
#' algorithm performed with function \code{nmf} with option
#' \code{keep.all=TRUE} (cf. \code{\link{nmf}}).
#' 
#' It extends both classes \code{\linkS4class{NMFfitX}} and \code{list}, and
#' stores the result of each run (i.e. a \code{NMFfit} object) in its
#' \code{list} structure.
#' 
#' IMPORTANT NOTE: This class is designed to be \strong{read-only}, even though
#' all the \code{list}-methods can be used on its instances. Adding or removing
#' elements would most probably lead to incorrect results in subsequent calls.
#' Capability for concatenating and merging NMF results is for the moment only
#' used internally, and should be included and supported in the next release of
#' the package.
#' 
#' 
#' @name NMFfitXn-class
#' @aliases NMFfitXn-class algorithm,NMFfitXn-method basis,NMFfitXn-method
#' coef,NMFfitXn-method modelname,NMFfitXn-method seeding,NMFfitXn-method
#' consensus,NMFfitXn-method entropy,NMFfitXn,ANY-method fit,NMFfitXn-method
#' minfit,NMFfitXn-method nrun,NMFfitXn-method purity,NMFfitXn,ANY-method
#' residuals,NMFfitXn-method seqtime,NMFfitXn-method show,NMFfitXn-method
#' @docType class
#' @section Slots: \describe{
#' 
#' \item{list(".Data")}{standard slot that contains the S3 \code{list} object
#' data.  See R documentation on S4 classes for more details.}\item{:}{standard
#' slot that contains the S3 \code{list} object data.  See R documentation on
#' S4 classes for more details.}
#' 
#' }
#' @author Renaud Gaujoux
#' @seealso \code{\linkS4class{NMFfitX}} , \code{\link{nmf-methods}},
#' \code{\link{nmf-multiple}}
#' @keywords classes
#' @examples
#' 
#' 
#' # generate a synthetic dataset with known classes
#' n <- 50; counts <- c(5, 5, 8);
#' V <- syntheticNMF(n, counts, noise=TRUE)
#' 
#' # build the class factor
#' groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))
#' 
#' # perform multiple runs of one algorithm, keeping all the fits
#' res <- nmf(V, 3, nrun=5, .options='k') # .options=list(keep.all=TRUE) also works
#' res
#' summary(res)
#' # get more info
#' summary(res, target=V, class=groups)
#' 
#' # compute/show computational times
#' runtime.all(res)
#' seqtime(res)
#' 
#' 
NULL





#' Virtual Class to Handle Results from Multiple Runs of a NMF Algorithms
#' 
#' This class defines a common interface to handle the results from multiple
#' runs of a single NMF algorithm, performed with the \code{\link{nmf}} method.
#' 
#' Currently, this interface is implemented by two classes,
#' \code{\linkS4class{NMFfitX1}} and \code{\linkS4class{NMFfitXn}}, which
#' respectively handle the case where only the best fit is kept, and the case
#' where the list of all the fits is returned.
#' 
#' See \code{\link{nmf-multiple}} for more details on the method arguments.
#' 
#' 
#' @name NMFfitX-class
#' @aliases NMFfitX-class consensus,NMFfitX-method featureNames,NMFfitX-method
#' sampleNames,NMFfitX-method nrun,NMFfitX-method show,NMFfitX-method
#' @docType class
#' @section Slots: \describe{
#' 
#' \item{list("runtime.all")}{Object of class \code{"proc_time"} that contains
#' various measures of the time spent to perform all the runs.}\item{:}{Object
#' of class \code{"proc_time"} that contains various measures of the time spent
#' to perform all the runs.}
#' 
#' }
#' @author Renaud Gaujoux
#' @seealso \code{\link{nmf-methods}}, \code{\link{nmf-multiple}},
#' \code{\linkS4class{NMFfitX1}}, \code{\linkS4class{NMFfitXn}}
#' @references
#' 
#' Brunet, J.P. et al. (2004) \emph{Metagenes and molecular pattern discovery
#' using matrix factorization} Proc Natl Acad Sci USA 101(12), 4164--4169.
#' 
#' Kim, H. and Park, H. (2007).  \emph{Sparse non-negative matrix
#' factorizations via alternating non-negativity-constrained least squares for
#' microarray data analysis}.  Bioinformatics 2007; \bold{23(12)}:1495-502.
#' \url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
#' @keywords classes
#' @examples
#' 
#' 
#' # generate a synthetic dataset with known classes
#' n <- 50; counts <- c(5, 5, 8);
#' V <- syntheticNMF(n, counts, noise=TRUE)
#' 
#' # perform multiple runs of one algorithm (default is to keep only best fit)
#' res <- nmf(V, 3, nrun=5)
#' str(res)
#' 
#' # plot a heatmap of the consensus matrix
#' \dontrun{ consensusmap(res) }
#' 
#' # perform multiple runs of one algorithm (keep all the fits)
#' res <- nmf(V, 3, nrun=5, .options='k')
#' str(res)
#' 
#' 
NULL





#' Undocumented or Internal methods and functions for package NMF
#' 
#' Dummy Rd file that gathers undocumented or internal methods and functions
#' for package NMF.
#' 
#' 
#' @aliases NMF-internals nmfAlgorithm nmfSeed
#' algorithm,NMFStrategyFunction-method
#' algorithm<-,NMFStrategyFunction,character-method
#' algorithm<-,NMFStrategyFunction,function-method show,NMFStrategy-method
#' show,NMFStrategyIterative-method modelname,NMFStrategy-method
#' name,NMFStrategy-method name<-,NMFStrategy,character-method
#' objective,NMFStrategy-method objective<-,NMFStrategy,character-method
#' objective<-,NMFStrategy,function-method name,NMFSeed-method
#' name<-,NMFSeed,character-method show,NMFSeed-method $,NMF-method
#' $<-,NMF-method nmfRegisterAlgorithm
#' nmfRegisterAlgorithm,NMFStrategy,character-method
#' nmfRegisterAlgorithm,NMFStrategy,missing-method
#' nmfRegisterAlgorithm,character,ANY-method
#' nmfRegisterAlgorithm,function,character-method std.divergence.update.h
#' std.divergence.update.w std.euclidean.update.h std.euclidean.update.w
#' nmf.stop.connectivity nmf.stop.iteration nmf.stop.stationary
#' nmf.stop.threshold
#' @author Renaud Gaujoux
#' @seealso \code{\linkS4class{NMF}}, \code{\link{nmf-methods}}, package's
#' vignette
#' @keywords internal
NULL





#' Class "NMFList" to Handle the Comparison of NMF Results
#' 
#' Class "NMFList" is used to wrap into a list the results of different NMF
#' runs with the objective to compare them.
#' 
#' While it handles indifferently any kind of NMF result, it is usually used to
#' compare NMF results from different algorithms.
#' 
#' 
#' @name NMFList-class
#' @aliases NMFList-class algorithm,NMFList-method runtime,NMFList-method
#' show,NMFList-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{as.NMFList(...)}.
#' @author Renaud Gaujoux
#' @seealso \code{\link{nmf-compare}}, \code{\link{nmf-multiple}},
#' \code{\link{nmf}}
#' @keywords classes
#' @examples
#' 
#' showClass("NMFList")
#' 
NULL





#' Main Interface to run NMF algorithms
#' 
#' This method implements the main interface to launch NMF algorithms within
#' the framework defined in package \code{NMF}. It allows to combine NMF
#' algorithms with seeding methods. The returned object can be directly passed
#' to visualisation or comparison methods.
#' 
#' For a tutorial on how to use the interface, please see the package's
#' vignette: \code{vignette('NMF')}
#' 
#' 
#' @name nmf-methods: Running NMF algorithms
#' @aliases nmf seed NMF Algorithms nmf-methods nmf,data.frame,ANY,ANY-method
#' nmf,ExpressionSet,ANY,ANY-method nmf,matrix,ANY,ANY-method
#' nmf,matrix,numeric,list-method nmf,matrix,numeric,character-method
#' nmf,matrix,numeric,function-method nmf,matrix,numeric,NMFStrategy-method
#' @docType methods
#' @param method The algorithm to use to perform NMF on \code{x}.  Different
#' formats are allowed: \code{character}, \code{function}.  If missing, the
#' method to use is retrieved from the NMF package's specific options by
#' \code{nmf.getOption("default.algorithm")} (the default built-in option is
#' \code{'brunet'}).  See section \emph{Methods} for more details on how each
#' format is used.
#' @param mixed Boolean that states if the algorithm requires a nonnegative
#' input matrix (\code{mixed=FALSE} which is the default value) or accepts
#' mixed sign input matrices (\code{mixed=TRUE}). An error is thrown if the
#' sign required is not fulfilled.  This parameter is useful to plug-in
#' algorithms such as semi-NMF, that typically does not impose nonnegativity
#' constraints on both the input and the basis component matrices.If
#' \code{NULL} then the NMF model is used
#' @param model
#' 
#' When \code{method} is a \code{function}, argument \code{model} must be
#' either a single \code{character} string (default to 'NMFstd') or a
#' \code{list} that specifies values for slots in the NMF model. The NMF model
#' to be instantiated can optionally be given by its class name in the first
#' element of the list [note: A NMF model is defined by a S4 class that extends
#' class \code{\linkS4class{NMF}}]. If no class name is specified then the
#' default model is used, see \code{\linkS4class{NMFstd}}.
#' 
#' When \code{method} is a single \code{character} string or a
#' \code{NMFStrategy} object, argument \code{model} must be \code{NULL}
#' (default) or a \code{list}.  Note that in this case the NMF model is defined
#' by the NMF strategy itself and cannot be changed.
#' 
#' \itemize{
#' 
#' \item If a single \code{character} string, argument \code{model} must be the
#' name of the class that defines the NMF model to be instantiated.  Arguments
#' in \code{...} are handled in the same way as when \code{model} is
#' \code{NULL}, see below.
#' 
#' \item If a \code{list} \code{all and only} its elements are used to
#' initialise the NMF model's slots.
#' 
#' \item if \code{NULL}, the arguments in \code{...} that have the same name as
#' slots in the NMF model associated with the NMF strategy of name
#' \code{method} are used to initialise these slots.
#' 
#' }
#' 
#' \strong{Important:} Values to initialise the NMF model's slots can be passed
#' in \code{\dots{}}.  However, if argument \code{model} is a \code{list} --
#' even empty -- then all and only its elements are used to initialise the
#' model, those in \code{...} are directly passed to the algorithm.
#' 
#' So to pass a parameter to the NMF algorithm, that has the same name as a
#' slot in the NMF model, argument \code{model} MUST be a list -- possibly
#' empty -- and contains all the values one wants to use for the NMF model
#' slots.
#' 
#' If a variable appears in both argument \code{model} and \code{\dots{}}, the
#' former will be used to initialise the NMF model, the latter will be passed
#' to the NMF algorithm.  See code examples for an illustration of this
#' situation.
#' @param name A \code{character} string to be used as a name for the custom
#' NMF algorithm.
#' @param nrun Used to perform multiple runs of the algorithm. It specifies the
#' number of runs to perform .  This argument is useful to achieve stability
#' when using a random seeding method.
#' @param objective Used when \code{method} is a \code{function}.  It must be A
#' \code{character} string giving the name of a built-in distance method or a
#' \code{function} to be used as the objective function.  It is used to compute
#' the residuals between the target matrix and its NMF estimate.
#' @param .callback Used when option \code{keep.all=FALSE} (default).  It
#' allows to pass a callback function that is called after each run when
#' performing multiple runs (i.e. with \code{nrun>1}).  This is useful for
#' example if one is also interested in saving summary measures or process the
#' result of each NMF fit before it gets discarded.  After each run, the
#' callback function is called with only one argument, the
#' \code{\linkS4class{NMFfit}} object that as just been fitted:
#' \code{.callback(res)} Therefore all other arguments should have default
#' values.
#' 
#' The results of the different calls to the callback function are stored in a
#' miscellaneous slot accessible by \code{result$.callback} (assuming one ran
#' \code{result <- nmf(...)}).  If no error occurs \code{result$.callback}
#' contains the list "simplified" by applying the \code{\link{sapply}}
#' function, which will try to convert a list with similar components into a
#' \code{vector}, a \code{matrix} or a \code{data.frame}.  If any error occurs
#' in one of the callback calls, then global computation is NOT stopped, but
#' the error is still stored in \code{result$.callback}, which is then a
#' \code{list}.
#' 
#' See the examples for a sample code.
#' @param .options this argument is used to set some runtime options. It can be
#' \code{list} containing the named options and their values, or, in the case
#' only boolean options need to be set, a character string that specifies which
#' options are turned on or off. The string must be composed of characters that
#' correspond to a given option. Characters '+' and '-' are used to explicitly
#' specify on and off respectively. E.g. \code{.options='tv'} will toggle on
#' options \code{track} and \code{verbose}, while \code{.options='t-v'} will
#' toggle on option \code{track} and off option \code{verbose}. Note that '+'
#' and '-' apply to all option character found after them. The default
#' behaviour is to assume that \code{.options} starts with a '+'.
#' 
#' The following options are available (note the characters that correspond to
#' each option, to be used when \code{.options} is passed as a string):
#' \describe{
#' 
#' \item{debug - d}{ Toggle debug mode. Like option \code{verbose} but with
#' more information displayed.}
#' 
#' \item{keep.all - k}{ used when performing multiple runs (\code{nrun}>1): if
#' toggled on, all factorizations are saved and returned, otherwise only the
#' factorization achieving the minimum residuals is returned.  }
#' 
#' \item{parallel - p}{ this option is useful on multicore *nix or Mac machine
#' only, when performing multiple runs (\code{nrun} > 1).  If toggled on, the
#' runs are performed using the parallel backend defined in argument
#' \code{.pbackend}. If this is set to \code{'mc'} then one tried to perform
#' the runs using multiple cores with package \code{link[package:doMC]{doMC}}
#' -- which therefore needs to be installed.
#' 
#' Unlike option 'P' (capital 'P'), if the computation cannot be performed in
#' parallel, then it will still be carried on sequentially.
#' 
#' \strong{IMPORTANT NOTE FOR MAC OS X USERS:} The parallel computation is
#' based on the \code{doMC} and \code{multicore} packages, so the same care
#' should be taken as stated in the vignette of \code{doMC}: \emph{\dQuote{it
#' is not safe to use doMC from R.app on Mac OS X. Instead, you should use doMC
#' from a terminal session, starting R from the command line.}} }
#' 
#' \item{parallel.required - P}{ Same as \code{p}, but an error is thrown if
#' the computation cannot be performed in parallel.  }
#' 
#' \item{restore.seed - r}{ used when seeding the NMF computation with a
#' numeric seed. When \code{TRUE} (default) the random seed
#' (\code{.Random.seed}) is restored to its value as before the call to the
#' \code{nmf} function.  }
#' 
#' \item{track - t}{ enables (resp. disables) error tracking. When \code{TRUE},
#' the returned object's slot \code{residuals} contains the trajectory of the
#' objective values.  This tracking functionality is available for all built-in
#' algorithms.  }
#' 
#' \item{verbose - v}{ Toggle verbosity. If on, messages about the
#' configuration and the state of the current run(s) are displayed.}
#' 
#' }% end describe .options
#' 

#' @param .pbackend define the parallel backend (from the
#' \code{\link[foreach]{foreach}} package) to use when running in parallel
#' mode.  See options \code{p} and \code{P} in argument \code{.options}.
#' Currently it accepts the following values: \code{'mc'} or a number that
#' specifies the number of cores to use, \code{'seq'} or \code{NULL} to use
#' sequential backend.
#' @param rank The factorization rank to achieve [i.e a single positive
#' \code{numeric}]
#' @param seed The seeding method to use to compute the starting point passed
#' to the algorithm.  See section \emph{Seeding methods} %and \code{function
#' \link{seed}} for more details on the possible classes and types for argument
#' \code{seed}.
#' @param x The target object to estimate. It can be a \code{matrix}, a
#' \code{data.frame}, an
#' \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}} object (this
#' requires the \code{Biobase} package to be installed).  See section
#' \emph{Methods} for more details.
#' @param \dots Extra parameters passed to the NMF algorithm's \code{run}
#' method or used to initialise the NMF model slots.  If argument \code{model}
#' is not supplied as a \code{list}, ANY of the arguments in \code{...} that
#' have the same name as slots in the NMF model to be instantiated will be used
#' to initialise these slots.  See also the \emph{Important} paragraph in
#' argument \code{model}.
#' @return The returned value depends on the run mode:
#' 
#' \item{Single run:}{An object that inherits from class
#' \code{\linkS4class{NMF}}.} \item{Multiple runs, single method:}{When
#' \code{nrun > 1} and \code{method} is NOT a \code{list}, this method returns
#' an object of class \code{\linkS4class{NMFfitX}}.} \item{Multiple runs,
#' multiple methods:}{When \code{nrun > 1} and \code{method} is a \code{list},
#' this method returns an object of class \code{\linkS4class{NMFList}}.}
#' @section Methods: \describe{
#' 
#' \item{x = "matrix", rank = "numeric", method = "list"}{ Performs NMF on
#' matrix \code{x} for each algorithm defined in the list \code{method}.}
#' 
#' \item{x = "data.frame", rank = "ANY", method = "ANY"}{ Performs NMF on a
#' \code{data.frame}: the target matrix is the converted data.frame
#' \code{as.matrix(x)} }
#' 
#' \item{x = "ExpressionSet", rank = "ANY", method = "ANY"}{ Performs NMF on an
#' \code{ExpressionSet}: the target matrix is the expression matrix
#' \code{exprs(x)}.
#' 
#' This method requires the \code{Biobase} package to be installed.  Special
#' methods for bioinformatics are provided in an optional layer, which is
#' automatically loaded when the \code{Biobase} is installed.  See
#' \code{\link{bioc-NMF}}. }
#' 
#' \item{x = "matrix", rank = "numeric", method = "character"}{ Performs NMF on
#' a \code{matrix} using an algorithm whose name is given by parameter
#' \code{method}.  The name provided must partially match the name of a
#' registered algorithm.  See section \emph{Algorithms} below or the package's
#' vignette for a list of the implemented algorithms and their respective
#' names.}
#' 
#' \item{x = "matrix", rank = "numeric", method = "function"}{ Performs NMF
#' using a custom algorithm defined by a \code{function}. It must have
#' signature \code{(x=matrix, start=NMF, ...)} and return an object that
#' inherits from class \code{NMF}. It should use its argument \code{start} as a
#' starting point.}
#' 
#' %\item{x = "matrix", rank = "numeric", method = "NMFStrategy"}{ Performs NMF
#' using an algorithm defined as a %\code{NMFStrategy} object. This version of
#' method \code{nmf} is the one that is eventually called by all the %above
#' versions.} }
#' @author Renaud Gaujoux
#' @seealso class \code{\linkS4class{NMF}}, \code{\link{utils-NMF}}, package's
#' vignette
#' @references
#' 
#' Lee, D.D. and Seung, H.S. (2000). \emph{Algorithms for non-negative matrix
#' factorization}.  In \emph{NIPS}, 556--562.
#' 
#' Brunet, J.P. et al. (2004). \emph{Metagenes and molecular pattern discovery
#' using matrix factorization}.  Proc Natl Acad Sci U S A, \bold{101}(12),
#' 4164--4169.\cr Original MATLAB code available from:\cr
#' \url{http://www.broadinstitute.org/cancer/pub/nmf}
#' 
#' Pascual-Montano, A. et al. (2006). \emph{Nonsmooth nonnegative matrix
#' factorization (nsnmf)}. IEEE transactions on pattern analysis and machine
#' intelligence, \bold{8}(3), 403--415.
#' 
#' Kim, H. and Park, H. (2007). \emph{Sparse non-negative matrix factorizations
#' via alternating non-negativity-constrained least squares for microarray data
#' analysis}. Bioinformatics. 2007; 23(12):1495-502.\cr Original MATLAB code
#' available from:\cr
#' \url{http://www.cc.gatech.edu/~hpark/software/nmfsh_comb.m}\cr
#' \url{http://www.cc.gatech.edu/~hpark/software/fcnnls.m}
#' 
#' Liviu Badea (2008). \emph{Extracting Gene Expression Profiles Common To
#' Colon And Pancreatic Adenocaricinoma Using Simultaneous Nonnegative Matrix
#' Factorization}. In Pacific Symposium on Biocomputing, \bold{13}, 279--290
#' 
#' %S. Li, X. Hou, and H. Zhang (2001). %Learning spatially localized,
#' parts-based representation.  %\emph{In Proc. CVPR}, 2001.
#' 
#' Zhang J, et al. (2008).  \emph{Pattern expression nonnegative matrix
#' factorization: algorithm and applications to blind source separation}.
#' Computational intelligence and neuroscience
#' 
#' C. Boutsidis and E. Gallopoulos (2007) \emph{SVD-based initialization: A
#' head start for nonnegative matrix factorization}.  Pattern Recognition.
#' doi:10.1016/j.patcog.2007.09.010\cr Original MATLAB code available from:\cr
#' \url{http://www.cs.rpi.edu/~boutsc/papers/paper1/nndsvd.m}
#' @keywords methods cluster math optimize
#' @examples
#' 
#' 
#' ## DATA
#' # generate a synthetic dataset with known classes: 100 features, 23 samples (10+5+8)
#' n <- 100; counts <- c(10, 5, 8); p <- sum(counts) 
#' V <- syntheticNMF(n, counts, noise=TRUE)
#' dim(V)
#' 
#' # build the class factor
#' groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))
#' 
#' ## RUN NMF ALGORITHMS
#' 
#' # run default algorithm
#' res <- nmf(V, 3)
#' res
#' summary(res, class=groups)
#' 
#' # run default algorithm multiple times (only keep the best fit)
#' res <- nmf(V, 3, nrun=10)
#' res
#' summary(res, class=groups)
#' 
#' # run default algorithm multiple times keeping all the fits
#' res <- nmf(V, 3, nrun=10, .options='k')
#' res
#' summary(res, class=groups)
#' 
#' \dontrun{
#' ## Note: one could have equivalently done
#' res <- nmf(V, 3, nrun=10, .options=list(keep.all=TRUE))
#' }
#' 
#' # run nonsmooth NMF algorithm
#' res <- nmf(V, 3, 'nsNMF')
#' res
#' summary(res, class=groups)
#' 
#' \dontrun{
#' ## Note: partial match also works
#' nmf(V, 3, 'ns')
#' }
#' 
#' \dontrun{
#' # Non default values for the algorithm's parameters can be specified in '...'
#' res <- nmf(V, 3, 'nsNMF', theta=0.8)
#' }
#' 
#' # compare some NMF algorithms (tracking the residual error)
#' res <- nmf(V, 3, list('brunet', 'lee', 'nsNMF'), seed=123456, .opt='t')
#' res
#' summary(res, class=groups)
#' # plot the track of the residual errors
#' \dontrun{plot(res)}
#' 
#' # run on an ExpressionSet (requires package Biobase)
#' \dontrun{
#' data(esGolub)
#' nmf(esGolub, 3)
#' }
#' 
#' ## USING SEEDING METHODS
#' 
#' # run default algorithm with the Non-negative Double SVD seeding method ('nndsvd')
#' nmf(V, 3, seed='nndsvd')
#' 
#' \dontrun{
#' ## Note: partial match also works
#' nmf(V, 3, seed='nn')
#' }
#' 
#' # run nsNMF algorithm, fixing the seed of the random number generator 
#' nmf(V, 3, 'nsNMF', seed=123456)
#' 
#' # run default algorithm specifying the starting point following the NMF standard model
#' start.std <- nmfModel(W=matrix(0.5, n, 3), H=matrix(0.2, 3, p))   
#' nmf(V, seed=start.std)
#' 
#' # to run nsNMF algorithm with an explicit starting point, this one
#' # needs to follow the 'NMFns' model:
#' start.ns <- nmfModel(model='NMFns', W=matrix(0.5, n, 3), H=matrix(0.2, 3, p))   
#' nmf(V, seed=start.ns)
#' # Note: the method name does not need to be specified as it is infered from the 
#' # when there is only one algorithm defined for the model.
#' 
#' # if the model is not appropriate (as defined by the algorihtm) an error is thrown 
#' # [cf. the standard model doesn't include a smoothing parameter used in nsNMF] 
#' \dontrun{nmf(V, method='ns', seed=start.std)}
#' 
#' ## Callback functions
#' # Pass a callback function to only save summary measure of each run
#' res <- nmf(V, 3, nrun=3, .callback=summary)
#' # the callback results are simplified into a matrix
#' res$.callback
#' 
#' # Pass a custom callback function
#' cb <- function(obj){ sparseness(obj) >= 0.5 }
#' res <- nmf(V, 3, nrun=3, .callback=cb)
#' res$.callback
#' 
#' # Passs a callback function which throws an error
#' cb <- function(){ i<-0; function(object){ i <<- i+1; if( i == 1 ) stop('SOME BIG ERROR'); summary(object) }}
#' res <- nmf(V, 3, nrun=3, .callback=cb())
#' 
#' 
NULL





#' Subset method for objects of class NMF
#' 
#' This method provide a convenient way of sub-setting objects of class
#' \code{NMF}, using a matrix-like syntax.
#' 
#' It allows to consistently subset one or both matrix factors in the NMF
#' model, as well as retrieving part of the basis components or part of the
#' mixture coefficients with a reduced amount of code.
#' 
#' Read section \emph{Value} for more details on the returned value.
#' 
#' 
#' @name Sub-setting NMF objects
#' @aliases subset-NMF [,NMF-method
#' @docType methods
#' @param x an object of class \code{\linkS4class{NMF}} to be subset.
#' @param i index used to subset on the \strong{rows} of the basis matrix (i.e.
#' the features).  It can be a \code{numeric}, \code{logical}, or
#' \code{character} vector (whose elements must match the row names of
#' \code{x}). In the case of a \code{logical} vector the entries are recycled
#' if necessary.
#' @param j index used to subset on the \strong{columns} of the mixture
#' coefficient matrix (i.e. the samples).  It can be a \code{numeric},
#' \code{logical}, or \code{character} vector (whose elements must match the
#' column names of \code{x}). In the case of a \code{logical} vector the
#' entries are recycled if necessary.
#' @param ...  used to specify a third index to subset on the basis components,
#' i.e. on both the columns and rows of the basis matrix and mixture
#' coefficient respectively.  It can be a \code{numeric}, \code{logical}, or
#' \code{character} vector (whose elements must match the basis names of
#' \code{x}). In the case of a \code{logical} vector the entries are recycled
#' if necessary.
#' 
#' Only the first extra subset index is used. A warning is thrown if more than
#' one extra argument is passed.
#' @param drop single \code{logical} value used to drop the \code{NMF-class}
#' wrapping and extract parts of the factor matrices.
#' 
#' When \code{drop=FALSE} it returns the \code{NMF} object \code{x} with the
#' basis matrix and/or mixture coefficient matrix subset accordingly to the
#' values in \code{i}, \code{j}, and \code{...}.
#' 
#' When \code{drop=TRUE} it returns the \code{matrix} that is subset "the more"
#' (see section \emph{Value}).
#' 
#' Note that in the case where both indexes \code{i} and \code{j} are provided,
#' argument \code{drop} is ignored: \code{x[i,j, drop=TRUE]} (resp.
#' \code{x[i,j,k, drop=TRUE]}) is identical to \code{x[i,j, drop=FALSE]} (resp.
#' \code{x[i,j,k, drop=FALSE]}).
#' @return
#' 
#' The returned value depends on the number of subset index passed and the
#' value of argument \code{drop}:
#' 
#' \itemize{ \item No index as in \code{x[]} or \code{x[,]}: the value is the
#' object \code{x} unchanged.
#' 
#' \item One single index as in \code{x[i]}: the value is the basis matrix
#' subset by \code{i}.  Precisely the call \code{x[i]} is equivalent to
#' \code{basis(x)[, i, drop=TRUE]}. If argument \code{drop} is present then it
#' is used: \code{x[i, drop=TRUE.or.FALSE]} <=> \code{basis(x)[, i,
#' drop=TRUE.or.FALSE]}.
#' 
#' \item More than one index with \code{drop=FALSE} (default) as in
#' \code{x[i,j]}, \code{x[i,]}, \code{x[,j]}, \code{x[i,j,k]}, \code{x[i,,k]},
#' etc...: the value is a \code{NMF} object whose basis and/or mixture
#' coefficient matrices have been subset accordingly. The third index \code{k}
#' affects simultaneously the columns of the basis matrix AND the rows of the
#' mixture coefficient matrix.
#' 
#' \item More than one index with \code{drop=TRUE} and \code{i} xor \code{j}
#' missing: the value returned is the matrix that is the more affected by the
#' subset index. That is that \code{x[i, , drop=TRUE]} and \code{x[i, , k,
#' drop=TRUE]} return the basis matrix subset by \code{[i,]} and \code{[i,k]}
#' respectively, while \code{x[, j, drop=TRUE]} and \code{x[, j, k, drop=TRUE]}
#' return the mixture coefficient matrix subset by \code{[,j]} and \code{[k,j]}
#' respectively.
#' 
#' }
#' @keywords methods
#' @examples
#' 
#' 
#' # create a dummy NMF object that highlight the different way of subsetting
#' a <- nmfModel(W=outer(seq(1,5),10^(0:2)), H=outer(10^(0:2),seq(-1,-10)))
#' basisnames(a) <- paste('b', 1:nbasis(a), sep='')
#' rownames(a) <- paste('f', 1:nrow(a), sep='')
#' colnames(a) <- paste('s', 1:ncol(a), sep='')
#' 
#' # or alternatively:
#' # dimnames(a) <- list( features=paste('f', 1:nrow(a), sep=''), samples=paste('s', 1:ncol(a), sep=''), basis=paste('b', 1:nbasis(a)) )
#' 
#' # look at the resulting NMF object 
#' a
#' basis(a)
#' coef(a)
#' 
#' # extract basis components
#' a[1]
#' a[1, drop=FALSE] # not dropping matrix dimension
#' a[2:3]
#' 
#' # subset on the features
#' a[1,]
#' a[2:4,]
#' # dropping the NMF-class wrapping => return subset basis matrix
#' a[2:4,, drop=TRUE]
#' 
#' # subset on the samples
#' a[,1]
#' a[,2:4]
#' # dropping the NMF-class wrapping => return subset coef matrix
#' a[,2:4, drop=TRUE]
#' 
#' # subset on the basis => subsets simultaneously basis and coef matrix
#' a[,,1]
#' a[,,2:3]
#' a[4:5,,2:3]
#' a[4:5,,2:3, drop=TRUE] # return subset basis matrix
#' a[,4:5,2:3, drop=TRUE] # return subset coef matrix
#' 
#' # 'drop' has no effect here
#' a[,,2:3, drop=TRUE]
#' 
NULL





#' Factory Method for NMF Models
#' 

#' 
#' \code{nmfModels} returns the types of NMF model currently defined that can
#' be used to create NMF objects.
#' 
#' \code{nmfModel} is a generic function which provides a convenient way to
#' build NMF models.
#' 
#' It provides a unique interface to create \code{NMF} objects that can follow
#' different NMF models, and is designed to resolve potential inconsistencies
#' in the matrices dimensions.
#' 

#' 
#' NMF model types are defined as S4 classes that inherit from class
#' \code{\linkS4class{NMF}}.
#' 
#' \code{nmfModel} methods act as factory methods to help in the creation of
#' NMF model objects in common situations: creating an empty model, a model
#' with given dimensions, a model with dimensions compatible with a given
#' target matrix, ...
#' 
#' All methods return an object that inherits from class \code{NMF}, except for
#' the call with no argument, which lists the NMF models defined in the session
#' (built-in and user-defined).
#' 
#' The returned \code{NMF} objects are suitable for seeding NMF algorithms via
#' argument \code{seed} of the \code{\link{nmf}} method.  In this case the
#' factorisation rank is implicitly set by the number of columns of the basis
#' vector matrix.
#' 
#' @name nmfModel - Factory method for NMF models
#' @aliases nmfModel nmfModel-methods nmfModel,missing,ANY-method
#' nmfModel,missing,missing-method nmfModel,NULL,ANY-method
#' nmfModel,numeric,matrix-method nmfModel,numeric,missing-method
#' nmfModel,numeric,numeric-method nmfModel,matrix,ANY-method
#' nmfModel,matrix,matrix-method nmfModel,ExpressionSet,ANY-method
#' nmfModel,ANY,ExpressionSet-method nmfModels
#' @docType methods
#' @param builtin.only If \code{TRUE}, only the models provided by the package
#' NMF itself are returned, discarding the user-defined models.
#' @param rank the target factorization rank
#' @param target the target matrix dimension
#' @param ncol the number of columns of the target matrix. This argument is
#' optional.  If not missing, it is used if \code{target} is not of length 2.
#' It takes precedence over the number of columns of \code{H} -- if this latter
#' is provided --, to define the target number of columns.
#' @param model the type of NMF model to instantiate, as the name of the
#' respective S4 class that derives from class \code{\linkS4class{NMF}}.
#' @param W the template basis matrix. Its dimensions must be compatible with
#' the target dimensions defined by \code{target} and/or \code{ncol} and/or
#' \code{H}. See argument \code{force.dim} for details on how its dimensions
#' are possibly reduced.
#' @param H the template mixture coefficient matrix.  Its dimensions must be
#' compatible with the target dimensions defined by \code{target} and/or
#' \code{ncol} and/or \code{W}. See argument \code{force.dim} for details on
#' how its dimensions are possibly reduced.
#' @param ...  extra \emph{named} arguments used to initialize other slots of
#' the NMF, specific to the class given in argument \code{model}.
#' @param force.dim if \code{TRUE} (default), the set of dimension of the
#' instantiated model \code{M} are chosen as follows:
#' 
#' \code{ nrow(M) = min(n[1], nrow(W)) ncol(M) = min(n[2], ncol, ncol(H))
#' nbasis(M) = min(rank, ncol(W), mrow(H)) }
#' 
#' The basis and mixture coefficient matrices are reduced accordingly, by
#' dropping extra rows and/or columns.
#' @param order.basis if \code{TRUE} (default) the rows of the mixture
#' coefficient matrix are ordered in the same order as the column of the basis
#' matrix. This is only relevant when \code{W} and \code{H} have the same
#' column and row names respectively, eventually in different orders.
#' @return For \code{nmfModels}: a character vector containing the names of the
#' S4 classes that define specific NMF models.
#' 
#' For \code{nmfModel}: an object that inherits from class
#' \code{\linkS4class{NMF}}.
#' @section Main Factory method:
#' 
#' The main factory engine of NMF models is implemented by the method with
#' signature \code{numeric, numeric}. Other factory methods provide convenient
#' ways of creating NMF models from e.g. a given target matrix or known
#' basis/coef matrices (see section \emph{Other Factory Methods}).
#' 
#' This factory method creates an object of class \code{model}, using the extra
#' parameters \code{...} to initialise slots that are specific to the given
#' model.
#' 
#' All NMF models implement get/set methods to access the matrix factors (see
#' \code{\link{basis}}), which can be initialised via arguments \code{W} and
#' \code{H}. For example, all the built-in models derive from class
#' \code{\linkS4class{NMFstd}}, which has two slots, \var{W} and \var{H}, to
#' hold the two factors.
#' 
#' If only argument \code{rank} is provided, the method creates a NMF model of
#' dimension 0x\code{rank}x0. That is that the basis and mixture coefficient
#' matrices, \var{W} and \var{H}, have dimension 0x\code{rank} and
#' \code{rank}x0 respectively.
#' 
#' If target dimensions are also provided in argument \code{target} as a
#' 2-length vector, then the method creates a \code{NMF} object compatible to
#' fit a target matrix of dimension \code{target[1]}x\code{target[2]}. That is
#' that the basis and mixture coefficient matrices, \var{W} and \var{H}, have
#' dimension \code{target[1]}x\code{rank} and \code{rank}x\code{target[2]}
#' respectively. The target dimensions can also be specified using both
#' arguments \code{target} and \code{ncol} to define the number of rows and the
#' number of columns of the target matrix respectively. If no other argument is
#' provided, these matrices are filled with NAs.
#' 
#' If arguments \code{W} and/or \code{H} are provided, the method creates a NMF
#' model where the basis and mixture coefficient matrices, \var{W} and \var{H},
#' are initialised using the values of \code{W} and/or \code{H}.
#' 
#' The dimensions given by \code{target}, \code{W} and \code{H}, must be
#' compatible. However if \code{force.dim=TRUE}, the method will reduce the
#' dimensions to the achieve dimension compatibility whenever possible.
#' 
#' When \code{W} and \code{H} are both provided, the \code{NMF} object created
#' is suitable to seed a NMF algorithm in a call to the \code{\link{nmf}}
#' method. Note that in this case the factorisation rank is implicitly set by
#' the number of basis vectors.
#' @keywords methods
#' @examples
#' 
#' 
#' # Build an empty NMF model 
#' nmfModel()
#' 
#' # create a NMF object based on one random matrix: the missing matrix is deduced
#' # Note this only works when using factory method NMF 
#' n <- 50; r <- 3; 
#' w <- rmatrix(n, r) 
#' nmfModel(W=w)
#' 
#' # create a NMF object based on random (compatible) matrices
#' p <- 20
#' h <- rmatrix(r, p)
#' nmfModel(W=w, H=h)
#' # or more conveniently without the argument names
#' nmfModel(w, h)
#' 
#' # Specify the type of NMF model (e.g. 'NMFns' for non-smooth NMF)
#' mod <- nmfModel(w, h, model='NMFns')
#' mod
#' 
#' # One can use an NMF model as a seed (initialization) when fitting a target matrix:
#' # NB: when called only with the seed the rank and NMF algorithm are selected based 
#' # on the input NMF model. e.g. here rank is r and the algorithm nsNMF is used
#' # See ?nmf. 
#' V <- rmatrix(n, p)
#' nmf(V, seed=mod)
#' 
#' # create an empty NMF model compatible with a given target matrix
#' nmfModel(V)
#' 
#' # create a r-ranked NMF model with a given target matrix
#' nmfModel(r, V)
#' 
#' # create a r-ranked NMF model with a given target dimensions n x p as a 2-length vector
#' nmfModel(r, c(n,p)) # directly
#' nmfModel(r, dim(V)) # or from an existing matrix <=> nmfModel(r, V)
#' # or alternatively passing each dimension separately
#' nmfModel(r, n, p)
#' 
#' # trying to create a NMF object based on incompatible matrices generates an error
#' h <- rmatrix(r+1, p)
#' \dontrun{
#' new('NMFstd', W=w, H=h)
#' nmfModel(w, h)
#' nmfModel(r+1, W=w, H=h)
#' }
#' 
#' # The factory method can be force the model to match some target dimensions
#' nmfModel(r, W=w, H=h)
#' nmfModel(r, n-1, W=w, H=h)
#' 
#' 
#' 
NULL





#' Handling Results from Multiple NMF Runs
#' 
#' The NMF package provides an easy way to perform multiple runs of a given NMF
#' algorithm on a target matrix.
#' 
#' The result from the \code{\link{nmf}} method is a
#' \code{\linkS4class{NMFfitX}} object that holds either all or only the best
#' run, depending on the running options:
#' 
#' \code{# keep only the best run} \code{object <- nmf(X, r, nrun=20)} \code{#
#' keep all the runs} \code{object <- nmf(X, r, nrun=20, .options='k')}
#' 
#' The methods documented here are used to handle such results. They are
#' usually independent of the type of result and can be used without change in
#' either situation (all runs kept or only the best one).
#' 
#' Note that when only the best result is kept, the result object conveniently
#' inherits from all the methods available for single runs. Therefore it can be
#' handled as if it had been computed by a single NMF run and all the methods
#' defined for such results can be used (cf. \code{\linkS4class{NMFfit}} and
#' \code{\link{utils-NMF}}).
#' 
#' See \code{\linkS4class{NMFfitXn}} and \code{\linkS4class{NMFfitX1}} for
#' details on the classes that implement respectively the case where all the
#' runs are kept and only the best run is kept.
#' 

#' 
#' \describe{
#' 
#' \item{consensus}{:
#' 
#' Computes the consensus matrix associated to the multiple NMF runs described
#' by \code{object}.  It is computed as the mean connectivity matrix of all the
#' runs.
#' 
#' It's been proposed by \emph{Brunet et al. (2004)} to help visualising and
#' measuring the stability of the clusters obtained by NMF approaches.
#' 
#' For objects of class \code{NMF} (e.g. results of a single NMF run, or NMF
#' models), the consensus matrix reduces to the connectivity matrix.
#' 
#' Note that argument \code{...} is not used.  }
#' 
#' \item{consensusmap}{\code{signature(object = "NMFfitX")}: plots a heatmap of
#' the consensus matrix of \code{object}.  See \code{\link{consensus}}.}
#' 
#' \item{cophcor}{ Computes the cophenetic correlation coefficient of consensus
#' matrix \code{object}, generally obtained from multiple NMF runs.
#' 
#' The cophenetic correlation coeffificient is based on the consensus matrix
#' (i.e. the average of connectivity matrices) and was proposed by \emph{Brunet
#' et al. (2004)} to measure the stability of the clusters obtained from NMF.
#' 
#' It is defined as the Pearson correlation between the samples' distances
#' induced by the consensus matrix (seen as a similarity matrix) and their
#' cophenetic distances from a hierachical clustering based on these very
#' distances (by default an average linkage is used).  See \emph{Brunet et al.
#' (2004)}.
#' 
#' Note that argument \code{...} is not used.  }
#' 
#' \item{dispersion}{ Computes the dispersion coefficient of consensus matrix
#' \code{object}, generally obtained from multiple NMF runs.
#' 
#' The dispersion coefficient is based on the consensus matrix (i.e. the
#' average of connectivity matrices) and was proposed by \emph{Kim and Park
#' (2007)} to measure the reproducibility of the clusters obtained from NMF .
#' It is defined as: \deqn{\rho = \sum_{i,j=1}^n 4 (C_{ij} - \frac{1}{2})^2 .}
#' , where \eqn{n} is the total number of samples.
#' 
#' We have \eqn{0 \leq \rho \leq 1} and \eqn{\rho = 1} only for a perfect
#' consensus matrix, where all entries 0 or 1.  A perfect consensus matrix is
#' obtained only when all the connectivity matrices are the same, meaning that
#' the algorithm gave the same clusters at each run.  See \emph{Kim and Park
#' (2007)}
#' 
#' Note that argument \code{...} is not used.  }
#' 
#' \item{nrun}{ returns the number of NMF runs performed to compute
#' \code{object}.
#' 
#' In the case of a \code{NMFfitXn} object it returns its length -- as it is
#' also a list.  In the case of a \code{NMFfitX1} object it returns the value
#' of its slot \code{nrun}.  In the case of a \code{NMFfit} object it always
#' returns 1 (this method exists to create a uniform access interface to NMF
#' results).  In the case of a \code{matrix} object, it looks for a value
#' attached the matrix as an attribute (e.g. by method
#' \code{\link{consensus}}). It returns NULL if no suitable value was found.
#' This is used to keep track of data about the parent fit and annotate plots.
#' 
#' Note that argument \code{...} is not used.  }
#' 
#' \item{plot.NMF.consensus}{: plots a heatmap of the consensus matrix
#' \code{x}. See \code{\link{consensus}}.}
#' 
#' \item{predict}{ returns a \code{factor} that gives the predicted cluster
#' index for each sample (resp. for each feature) based on the \emph{best NMF
#' factorization} stored in \code{object}.
#' 
#' The index correspond to the basis vector that most contributes to the sample
#' (resp. to which the feature contributes the most).  See
#' \code{\link{predict}} for more details.  }
#' 
#' \item{runtime.all}{: returns the computational time used to compute all the
#' runs and create \code{object}.  The time is computed using base function
#' \code{\link{system.time}} which returns object of class
#' \code{\link[=proc.time]{proc_time}}.
#' 
#' For \code{NMFfitXn} objects, there is also another time measure returned by
#' the \code{seqtime} method, which computes the sequential computational time,
#' that is the sum of the computational time used by each run.
#' 
#' Note that argument \code{...} is not used.  }
#' 
#' \item{seqtime}{: returns the sequential CPU time spent of all the runs in
#' the \code{object} -- which must be an instance of class
#' \code{\linkS4class{NMFfitXn}}.  It is the sum of the CPU time used to
#' compute each run. It returns \code{NULL} if the \code{object} is empty.
#' 
#' Note that argument \code{...} is not used.  }
#' 
#' \item{summary}{: \code{summary} method for objects of class \code{NMFfitX}.
#' 
#' It computes a set of measures to help evaluate the quality of the \emph{best
#' factorization} of the set. The result is similar to the result from the
#' \code{summary} method of \code{NMFfit} objects. See \code{\linkS4class{NMF}}
#' for details on the computed measures.  In addition, the cophenetic
#' correlation coefficient and the dispersion coefficient of the consensus
#' matrix are returned, as well as the total computational time.  See the
#' related methods above.  }
#' 
#' }
#' 
#' @name Handling the results of multiple NMF runs
#' @aliases nmf-multiple consensus consensus-methods cophcor cophcor-methods
#' cophcor,matrix-method cophcor,NMFfitX-method dispersion dispersion-methods
#' dispersion,matrix-method dispersion,NMFfitX-method plot.NMF.consensus nrun
#' nrun-methods nrun,matrix-method predict,NMFfitXn-method runtime.all
#' runtime.all-methods runtime.all,NMFfitX-method runtime.all,NMFfitXn-method
#' seqtime seqtime-methods summary,NMFfitX-method
#' @docType methods
#' @param null used in method \code{runtime.all} for \code{NMFfitXn} objects to
#' specify if the result should be \code{NULL} when the object has no time data
#' is stored the total computation time. In this case, if \code{null=FALSE}
#' (default), the method returns the sequential time (cf. \code{seqtime} below)
#' instead of \code{NULL}. It also emits a warning which can be toggle with
#' argument \code{warning}.
#' @param object A \code{matrix} or an object that inherits from class
#' \code{\linkS4class{NMFfitX}} or \code{\linkS4class{NMFfit}} -- depending on
#' the method.
#' @param warning used in method \code{runtime.all} for \code{NMFfitXn} objects
#' to specify if a warning should be emitted when the object has no time data
#' the total computation time and the sequential time is returned instead of
#' \code{NULL} (cf. argument \code{null}).
#' @param x An object that inherits from the S3 class \code{NMF.consensus} as
#' returned by the method \code{\link{consensus}}. Technically it is nothing
#' else than a \code{matrix}.
#' @param ...  Used to pass extra arguments to subsequent calls: \itemize{
#' \item in \code{consensusmap}: Used to pass extra parameters to the
#' subsequent call to the heatmap drawing function \code{\link{aheatmap}}.
#' \item in \code{predict}: extra arguments passed to function
#' \code{\link{predict,NMF-method}} \item in \code{summary}: extra arguments
#' like \code{target} or \code{class} passed to the method
#' \code{\link{summary,NMFfit-method}}.  }
#' @author Renaud Gaujoux
#' @seealso \linkS4class{NMFfitX1}, \linkS4class{NMFfitXn},
#' \link[=summary,NMF-method]{summary}
#' @references
#' 
#' Brunet, J.P. et al. (2004).  \emph{Metagenes and molecular pattern discovery
#' using matrix factorization}.  Proc Natl Acad Sci USA 101(12), 4164--4169.
#' 
#' Kim, H. and Park, H. (2007).  \emph{Sparse non-negative matrix
#' factorizations via alternating non-negativity-constrained least squares for
#' microarray data analysis}.  Bioinformatics 2007; \bold{23(12)}:1495-502.
#' \url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
#' @examples
#' 
#' 
#' 
#' # generate a synthetic dataset with known classes: 50 features, 18 samples (5+5+8)
#' n <- 50; counts <- c(5, 5, 8);
#' V <- syntheticNMF(n, counts, noise=TRUE)
#' \dontrun{aheatmap(V)}
#' 
#' # build the class factor
#' groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))
#' 
#' # perform multiple runs of NMF (keep best only)
#' res <- nmf(V, 3, nrun=5)
#' res
#' 
#' # compute summary measures
#' summary(res)
#' 
#' # compute more summary measures
#' summary(res, target=V, class=groups)
#' 
#' # plot a heatmap of the consensus matrix with extra annotations
#' \dontrun{consensusmap(res, annCol=groups)}
#' 
#' # retrieve the predicted clusters of samples
#' predict(res)
#' 
#' # perform multiple runs of NMF and keep all the runs
#' res <- nmf(V, 3, nrun=5, .options='k')
#' res
#' 
#' # extract best NMF model
#' fit(res)
#' # extract best result object
#' minfit(res)
#' 
#' # compute/show computational times
#' runtime.all(res)
#' seqtime(res)
#' 
#' 
#' 
NULL





#' Nonsmooth Nonnegative Matrix Factorization
#' 
#' Class that implements the \emph{Nonsmooth Nonnegative Matrix Factorization}
#' (nsNMF) model, required by the Nonsmooth NMF algorithm.
#' 
#' The Nonsmooth NMF algorithm is defined by Pascual-Montano et al. (2006) as a
#' modification of the standard divergence based NMF algorithm (see section
#' Details and references below).  It aims at obtaining sparser factor
#' matrices, by the introduction of a smoothing matrix.
#' 

#' 
#' The Nonsmooth NMF algorithm is a modification of the standard divergence
#' based NMF algorithm (see \code{\linkS4class{NMF}}).  Given a non-negative
#' \eqn{n \times p}{n x p} matrix \eqn{V} and a factorization rank \eqn{r}, it
#' fits the following model: \deqn{V \equiv W S(\theta) H,}{V &#126; W S(theta)
#' H,} where: \itemize{ \item \eqn{W} and \eqn{H} are such as in the standard
#' model, that is non-negative matrices of dimension \eqn{n \times r}{n x r}
#' and \eqn{r \times p}{r x p} respectively; \item \eqn{S} is a \eqn{r \times
#' r} square matrix whose entries depends on an extra parameter \eqn{0\leq
#' \theta \leq 1} in the following way: \deqn{S = (1-\theta)I +
#' \frac{\theta}{r} 11^T , } where \eqn{I} is the identity matrix and \eqn{1}
#' is a vector of ones. }
#' 
#' The interpretation of S as a smoothing matrix can be explained as follows:
#' Let \eqn{X} be a positive, nonzero, vector. Consider the transformed vector
#' \eqn{Y = S X}. If \eqn{\theta = 0}, then \eqn{Y = X} and no smoothing on
#' \eqn{X} has occurred.  However, as \eqn{\theta \to 1}{theta tends to 1}, the
#' vector \eqn{Y} tends to the constant vector with all elements almost equal
#' to the average of the elements of \eqn{X}. This is the smoothest possible
#' vector in the sense of non-sparseness because all entries are equal to the
#' same nonzero value, instead of having some values close to zero and others
#' clearly nonzero.
#' 
#' @name NMFns-class: Model class for Non-smooth NMF
#' @aliases NMFns-class fitted,NMFns-method show,NMFns-method smoothing
#' smoothing-methods smoothing,NMFns-method
#' @docType class
#' @param x an object of class \code{NMFns}
#' @param object an object of class \code{NMFns}
#' @param W the \code{matrix} of basis vectors, i.e. the first matrix factor in
#' the non-smooth NMF model.
#' @param H the \code{matrix} of mixture coefficients, i.e. the third matrix
#' factor the non-smooth NMF model.
#' @param S the smoothing \code{matrix}, i.e. the middle matrix factor in the
#' non-smooth NMF model.
#' @param theta a single \code{numeric}
#' @param ...  extra parameters passed to method \code{smoothing}. So typically
#' used to pass a value for \code{theta}.
#' @section Algorithm:
#' 
#' The Nonsmooth NMF algorithm uses a modified version of the multiplicative
#' update equations in Lee & Seung's method for Kullback-Leibler divergence
#' minimization.  The update equations are modified to take into account the --
#' constant -- smoothing matrix.  The modification reduces to using matrix
#' \eqn{W S} instead of matrix \eqn{W} in the update of matrix \eqn{H}, and
#' similarly using matrix \eqn{S H} instead of matrix \eqn{H} in the update of
#' matrix \eqn{W}.  %See \code{\link{NMF Algorithms}} for more details on the
#' built-in NMF algorithms.
#' 
#' After matrix \eqn{W} have been updated, each of its columns is scaled so
#' that it sums up to 1.
#' @author Renaud Gaujoux
#' @seealso \code{\linkS4class{NMF}} , \code{\link{nmf-methods}}
#' @references
#' 
#' Alberto Pascual-Montano et al. (2006). Nonsmooth Nonnegative Matrix
#' Factorization (nsNMF). \emph{IEEE Transactions On Pattern Analysis And
#' Machine Intelligence} , Vol. 28, No. 3, March 2006 403
#' @keywords classes
#' @examples
#' 
#' 
#' # create a completely empty NMF object
#' new('NMFns')
#' 
#' # create a NMF object based on random (compatible) matrices
#' n <- 50; r <- 3; p <- 20
#' w <- rmatrix(n, r) 
#' h <- rmatrix(r, p)
#' nmfModel(model='NMFns', W=w, H=h)
#' 
#' # apply Nonsmooth NMF algorithm to a random target matrix
#' V <- rmatrix(n, p)
#' \dontrun{nmf(V, r, 'ns')}
#' 
#' 
#' 
NULL





#' Nonnegative Matrix Factorization with Offset
#' 
#' Class that implements the \emph{Nonnegative Matrix Factorization with
#' Offset} model, required by the NMF with Offset algorithm.
#' 
#' The NMF with Offset algorithm is defined by Badea (2008) as a modification
#' of Lee & Seung's euclidean based NMF algorithm (see section Details and
#' references below).  It aims at obtaining 'cleaner' factor matrices, by the
#' introduction of an offset matrix, explicitly modelling a feature specific
#' baseline -- constant across samples.
#' 
#' 
#' @name NMFOffset-class: Model class for NMF with offset
#' @aliases NMFOffset-class fitted,NMFOffset-method initialize,NMFOffset-method
#' offset,NMFOffset-method rnmf,NMFOffset,numeric-method show,NMFOffset-method
#' @docType class
#' @section Objects from the Class:
#' 
#' Object of class \code{NMFOffset} can be created using the standard way with
#' operator \code{\link{new}}
#' 
#' However, as for all the classes that extend class
#' \code{\linkS4class{NMFstd}}, objects of class \code{NMFOffset} should be
#' created using factory method \code{\link{nmfModel}} :
#' 
#' \code{new('NMFOffset')}
#' 
#' \code{nmfModel(model='NMFOffset')}
#' 
#' \code{nmfModel(model='NMFOffset', W=w, offset=rep(1, nrow(w)))}
#' 
#' See \code{\link{nmfModel}} for more details on how to use the factory
#' method.
#' @author Renaud Gaujoux
#' @seealso \code{\linkS4class{NMF}}, \code{\link{nmf-methods}}
#' @references Badea (2008).  Extracting Gene Expression Profiles Common To
#' Colon And Pancreatic Adenocaricinoma Using Simultaneous Nonnegative Matrix
#' Factorization.  In \emph{Pacific Symposium on Biocomputing} , \bold{13},
#' 279-290
#' @keywords classes
#' @examples
#' 
#' 
#' # create a completely empty NMF object
#' new('NMFOffset')
#' 
#' # create a NMF object based on random (compatible) matrices
#' n <- 50; r <- 3; p <- 20
#' w <- rmatrix(n, r) 
#' h <- rmatrix(r, p)
#' nmfModel(model='NMFOffset', W=w, H=h, offset=rep(0.5, nrow(w)))
#' 
#' # apply Nonsmooth NMF algorithm to a random target matrix
#' V <- rmatrix(n, p)
#' \dontrun{nmf(V, r, 'offset')}
#' 
#' 
#' 
NULL





#' NMF Package Overview
#' 
#' The NMF package provides methods to perform Nonnegative Matrix Factorization
#' (NMF) , as well as a framework to develop and test new NMF algorithms.
#' 
#' A number of standard algorithms and seeding methods are implemented. Tuned
#' visualisation and post-analysis methods help in the evaluation of the
#' algorithms' performances or in the interpretation of the results.
#' 
#' 
#' @name NMF-package
#' @aliases NMF-package NMF
#' @docType package
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @seealso \code{\link{NMF-class}}, \code{\link{nmf}},
#' \code{\link[Biobase]{Biobase}}
#' @references
#' 
#' If you use the 'NMF' package in publications cite:
#' 
#' Renaud Gaujoux, Cathal Seoighe (2010).\cr \emph{A flexible R package for
#' nonnegative matrix factorization}.\cr BMC Bioinformatics 2010,
#' \bold{11}:367\cr \url{http://www.biomedcentral.com/1471-2105/11/367}
#' 
#' For pointers on further readings, please see references therein.
#' 
#' Definition of Nonnegative Matrix Factorization in its modern formulation:
#' 
#' Lee D.D. and Seung H.S. (1999).  Learning the parts of objects by
#' non-negative matrix factorization.  \emph{Nature}, \bold{401}, 788--791.
#' 
#' Historical first definition and algorithms:
#' 
#' Paatero, P., Tapper, U. (1994).  Positive matrix factorization: A
#' non-negative factor model with optimal utilization of error estimates of
#' data values.  \emph{Environmetrics}, \bold{2}, 111--126 ,
#' doi:10.1002/env.3170050203.
#' @keywords package package
#' @examples
#' 
#' 
#' # run default NMF algorithm on a random matrix
#' V <- rmatrix(500, 20)
#' res <- nmf(V, 3)
#' res
#' 
#' # compute some quality measures
#' summary(res)
#' 
#' # Visualize the results as heatmaps
#' \dontrun{coefmap(res)} # mixture coefficients
#' \dontrun{basismap(res)} # basis vectors
#' 
#' # run default NMF algorithm on a random matrix with actual patterns
#' set.seed(123456)
#' V <- syntheticNMF(500, 3, 20, noise=TRUE)
#' res <- nmf(V, 3)
#' res
#' 
#' # compute some quality measures
#' summary(res)
#' 
#' # Visualize the results as heatmaps
#' \dontrun{coefmap(res)} # mixture coefficients
#' \dontrun{basismap(res)} # basis vectors
#' 
#' 
#' 
NULL





#' Deprecated Class to store results from multiple runs of NMF algorithms
#' 
#' This class is deprecated and replaced by class \code{\linkS4class{NMFfitX}}
#' and its extensions. It remains only for backward compatibility and will be
#' defunct in the next release.
#' 
#' It extends the base class \code{list} to store the result from a multiple
#' run of NMF algorithms.
#' 
#' The elements are of class \code{NMF}.
#' 
#' 
#' @name NMFSet-class
#' @docType class
#' @section Slots: \describe{ \item{list("consensus")}{Object of class
#' \code{"matrix"} used to store the consensus matrix when multiple runs have
#' been performed with option \code{keep.all=FALSE}. In this case, only the
#' best factorization is returned, so the object is of length 1. However the
#' consensus matrix across all runs is still computed and stored in this
#' slot.}\item{:}{Object of class \code{"matrix"} used to store the consensus
#' matrix when multiple runs have been performed with option
#' \code{keep.all=FALSE}. In this case, only the best factorization is
#' returned, so the object is of length 1. However the consensus matrix across
#' all runs is still computed and stored in this slot.}
#' 
#' \item{list("nrun")}{an \code{integer} that contains the number of runs when
#' NMF is performed with option \code{keep.all=FALSE}.
#' 
#' See \code{\link{nmf}}.  }\item{:}{an \code{integer} that contains the number
#' of runs when NMF is performed with option \code{keep.all=FALSE}.
#' 
#' See \code{\link{nmf}}.  }
#' 
#' \item{list("runtime")}{Object of class \code{"proc_time"} that contains
#' various measures of the time spent to perform all the runs.}\item{:}{Object
#' of class \code{"proc_time"} that contains various measures of the time spent
#' to perform all the runs.}
#' 
#' \item{list(".Data")}{standard slot that contains the S3 \code{list} object
#' data.  See R documentation on S4 classes for more details.}\item{:}{standard
#' slot that contains the S3 \code{list} object data.  See R documentation on
#' S4 classes for more details.}
#' 
#' }
#' @author Renaud Gaujoux
#' @seealso \code{\linkS4class{NMFfitX}}, \code{\linkS4class{NMF}} ,
#' \code{\link{nmf-methods}}, \code{\link{nmf-multiple}}
#' @keywords classes
NULL





#' Implement of the standard NMF model
#' 
#' Class that implements the standard model of Nonnegative Matrix
#' Factorisation.
#' 
#' It provides a general structure and generic functions to manage
#' factorizations that follow NMF standard model.
#' 
#' Let \eqn{V} be a \eqn{n \times m} non-negative matrix and \eqn{r} a positive
#' integer.  In its standard form (see references below), a NMF of \eqn{V} is
#' commonly defined as a pair of matrices \eqn{(W, H)} such that:
#' 
#' \deqn{V \equiv W H,}
#' 
#' where: \itemize{ \item \eqn{W} and \eqn{H} are \eqn{n \times r} and \eqn{r
#' \times m} matrices respectively with non-negative entries; \item
#' \eqn{\equiv} is to be understood with respect to some loss function.  Common
#' choices of loss functions are based on Frobenius norm or Kullback-Leibler
#' divergence.  }
#' 
#' Integer \eqn{r} is called the \emph{factorization rank}.  Depending on the
#' context of application of NMF, the columns of \eqn{W} and \eqn{H} take
#' different names: \describe{ \item{columns of }{basis vector, metagenes,
#' factors, source, image basis}\item{list(list("W"))}{basis vector, metagenes,
#' factors, source, image basis} \item{columns of }{mixture coefficients,
#' metagenes expression profiles, weights}\item{list(list("H"))}{mixture
#' coefficients, metagenes expression profiles, weights} }
#' 
#' NMF approach has been successfully applied to several fields. Package NMF
#' was implemented trying to use names as generic as possible for objects and
#' methods.  The following terminology is used: \describe{ \item{samples}{the
#' columns of the target matrix \eqn{V}} \item{features}{the rows of the target
#' matrix \eqn{V}} \item{basis matrix}{the first matrix factor \eqn{W}}
#' \item{basis vectors}{the columns of first matrix factor \eqn{W}}
#' \item{mixture matrix}{the second matrix factor \eqn{H}} \item{mixtures
#' coefficients}{the columns of second matrix factor \eqn{H}} }
#' 
#' However, because package NMF was primilary implemented to work with gene
#' expression microarray data, it also provides a layer to easily and
#' intuitively work with objects from the Bioconductor base framework.  See
#' \link{bioc-NMF} for more details.
#' 

#' 
#' @name NMFstd-class: Model class for standard NMF
#' @aliases NMFstd-class fitted,NMFstd-method basis,NMFstd-method
#' basis<-,NMFstd,matrix-method coef,NMFstd-method coef<-,NMFstd,matrix-method
#' @docType class
#' @section Slots: \describe{ \item{list("W")}{A \code{"matrix"} that contains
#' the \emph{first} matrix factor of the factorisation }\item{:}{A
#' \code{"matrix"} that contains the \emph{first} matrix factor of the
#' factorisation } \item{list("H")}{A \code{"matrix"} that contains the
#' \emph{second} matrix factor of the factorisation }\item{:}{A \code{"matrix"}
#' that contains the \emph{second} matrix factor of the factorisation } }
#' @author Renaud Gaujoux
#' @seealso Main interface to perform NMF in \code{\link{nmf-methods}}.
#' 
#' Method \code{\link{seed}} to set NMF objects with values suitable to start
#' algorithms with.
#' @references
#' 
#' Definition of Nonnegative Matrix Factorization in its modern formulation:
#' 
#' Lee D.D. and Seung H.S. (1999).  Learning the parts of objects by
#' non-negative matrix factorization.  \emph{Nature}, \bold{401}, 788--791.
#' 
#' Historical first definition and algorithms:
#' 
#' Paatero, P., Tapper, U. (1994).  Positive matrix factorization: A
#' non-negative factor model with optimal utilization of error estimates of
#' data values.  \emph{Environmetrics}, \bold{2}, 111--126 ,
#' doi:10.1002/env.3170050203.
#' 
#' Reference for some utility functions:
#' 
#' Kim, H. and Park, H. (2007).  Sparse non-negative matrix factorizations via
#' alternating non-negativity-constrained least squares for microarray data
#' analysis.  \emph{Bioinformatics}.
#' 
#' Hoyer (2004).  Non-negative matrix factorization with sparseness
#' constraints.  \emph{Journal of Machine Learning Research}, \bold{5},
#' 1457-1469.
#' @keywords classes
#' @examples
#' 
#' 
#' # create a completely empty NMF object (i.e. 0 features, 0 basis components, 0 samples)
#' new('NMFstd')
#' 
#' # create a NMF object based on one random matrix: the missing matrix is deduced
#' # Note this only works when using factory method NMF 
#' n <- 50; r <- 3; 
#' w <- rmatrix(n, r) 
#' nmfModel(W=w)
#' 
#' # create a NMF object based on random (compatible) matrices
#' p <- 20
#' h <- rmatrix(r, p)
#' nmfModel(W=w, H=h)
#' 
#' # create a NMF object based on incompatible matrices: generate an error
#' h <- rmatrix(r+1, p)
#' \dontrun{
#' new('NMFstd', W=w, H=h)
#' nmfModel(w, h)
#' }
#' 
#' # Giving target dimensions to the factory method allow for coping with dimension
#' # incompatibilty (a warning is thrown in such case) 
#' nmfModel(r, W=w, H=h)
#' 
#' 
#' 
#' 
NULL





#' Class and Utility Methods for NMF objects
#' 
#' Define generic interface methods for class \code{\linkS4class{NMF}}, which
#' is the base -- virtual -- class of the results from any NMF algorithms
#' implemented within package NMF's framework.
#' 

#' 
#' \describe{
#' 
#' \item{connectivity}{ Computes the connectivity matrix for the samples based
#' on their mixture coefficients.
#' 
#' The connectivity matrix of a clustering is a matrix \eqn{C} containing only
#' 0 or 1 entries such that: \deqn{C_{ij} = \left\{\begin{array}{l}1\mbox{ if
#' sample }i\mbox{ belongs to the same cluster as sample }j\\0\mbox{
#' otherwise}\end{array}\right..}{% C_{ij} = 1 if sample i belongs to the same
#' cluster as sample j, 0 otherwise}
#' 
#' }
#' 
#' \item{entropy}{ The entropy is a measure of performance of a clustering
#' method, in recovering classes defined by factor a priori known (i.e. one
#' knows the true class labels).  Suppose we are given \eqn{l} categories,
#' while the clustering method generates \eqn{k} clusters. Entropy is given by:
#' \deqn{Entropy = - \frac{1}{n \log_2 l} \sum_{q=1}^k \sum_{j=1}^l n_q^j
#' \log_2 \frac{n_q^j}{n_q}} , where:
#' 
#' - \eqn{n} is the total number of samples;
#' 
#' - \eqn{n} is the total number of samples in cluster \eqn{q};
#' 
#' - \eqn{n_q^j} is the number of samples in cluster \eqn{q} that belongs to
#' original class \eqn{j} (\eqn{1 \leq j \leq l}).
#' 
#' The smaller the entropy, the better the clustering performance.
#' 
#' See \emph{Kim and Park (2007)}.  }
#' 
#' \item{evar}{ Computes the explained variance of the NMF model \code{object}.
#' 
#' For a target \eqn{V} It is defined as: \deqn{evar = 1 -
#' \frac{RSS}{\sum_{i,j} v_{ij}^2}}{ evar = 1 - RSS/sum v_{ij}^2},
#' 
#' where RSS is the residual sum of squares.
#' 
#' It is usefull to compare the performance of different models and their
#' ability to accurately reproduce the original target matrix.  Note that a
#' possible caveat is that some methods explicitly aim at minimizing the RSS
#' (i.e. maximizing the explained variance), while others do not.  }
#' 
#' \item{extractFeatures}{ Identify the most basis-specific features, using
#' different methods.  See details of argument \code{method}.  }
#' 
#' \item{featureScore}{ Computes the feature scores as suggested in \emph{Kim
#' and Park (2007)}.
#' 
#' The score for feature \eqn{i} is defined as: \deqn{S_i = 1 + \frac{1}{\log_2
#' k} \sum_{q=1}^k p(i,q) \log_2 p(i,q),} where \eqn{p(i,q)} is the probability
#' that the \eqn{i}-th feature contributes to basis \eqn{q}: \deqn{p(i,q) =
#' \frac{W(i,q)}{\sum_{r=1}^k W(i,r)} }
#' 
#' The feature scores are real values within the range [0,1].  The higher the
#' feature score the more basis-specific the corresponding feature.
#' 
#' }
#' 
#' \item{hasTrack}{Tells if an \code{NMFfit} object has a recorded error track.
#' It should return \code{TRUE} if \code{object} was computed with the option
#' \code{track=TRUE} or the flag \code{'t'}.}
#' 
#' \item{nmfApply}{ \code{apply}-like method for objects of class \code{NMF}.
#' 
#' When argument \code{MARGIN=1}, it calls the base method \code{apply} to
#' apply function \code{FUN} to the \emph{rows} of the basis component matrix.
#' 
#' When \code{MARGIN=2}, it calls the base method \code{apply} to apply
#' function \code{FUN} on the \emph{columns} of the mixture coefficient matrix.
#' 
#' When \code{MARGIN=3}, it calls the base method \code{apply} to apply
#' function \code{FUN} on the \emph{columns} of the basis matrix.
#' 
#' When \code{MARGIN=4}, it calls the base method \code{apply} to apply
#' function \code{FUN} on the \emph{rows} of the mixture coefficient matrix.
#' 
#' See \code{\link[base]{apply}} for more details on the output format.  }
#' 
#' \item{plot}{ plots the residuals track of the run that computed object
#' \code{x}.  See function \code{\link{nmf}} for details on how to enable the
#' tracking of residuals.}
#' 
#' \item{purity}{ Computes the purity of a clustering given a known factor.
#' 
#' The purity is a measure of performance of a clustering method, in recovering
#' the classes defined by a factor a priori known (i.e. one knows the true
#' class labels).  Suppose we are given \eqn{l} categories, while the
#' clustering method generates \eqn{k} clusters. Purity is given by:
#' \deqn{Purity = \frac{1}{n} \sum_{q=1}^k \max_{1 \leq j \leq l} n_q^j} ,
#' where:
#' 
#' - \eqn{n} is the total number of samples;
#' 
#' - \eqn{n_q^j} is the number of samples in cluster \eqn{q} that belongs to
#' original class \eqn{j} (\eqn{1 \leq j \leq l}).
#' 
#' The purity is therefore a real number in \eqn{[0,1]}.  The larger the
#' purity, the better the clustering performance.
#' 
#' See \emph{Kim and Park (2007)}.  }
#' 
#' \item{randomize}{ permute the row entries within each column of \code{x},
#' using a different permutation for each column.
#' 
#' The extra arguments in \code{...} are passed to the \code{sample} function,
#' and will be used for each column.
#' 
#' The result is a \code{matrix} of the same dimension as \code{x} (or
#' \code{exprs(x)} in the case \code{x} is an \code{ExpressionSet} object.).  }
#' 
#' \item{rss}{ returns the Residual Sum of Squares (RSS) between the target
#' object \code{target} and its estimation by the \code{object}. \emph{Hutchins
#' et al. (2008)} used the variation of the RSS in combination with \emph{Lee
#' and Seung}'s algorithm to estimate the correct number of basis vectors. The
#' optimal rank is chosen where the graph of the RSS first shows an inflexion
#' point. See references.
#' 
#' Note that this way of estimation may not be suitable for all models. Indeed,
#' if the NMF optimization problem is not based on the Frobenius norm, the RSS
#' is not directly linked to the quality of approximation of the NMF model.
#' 
#' }
#' 
#' \item{sparseness}{ Generic mathod that computes the sparseness of an object
#' as defined in \emph{Hoyer (2004)}.
#' 
#' This sparseness measure quantifies how much energy of a vector is packed
#' into only few components.  It is defined by: \deqn{Sparseness(x) =
#' \frac{\sqrt{n} - \frac{\sum |x_i|}{\sqrt{\sum x_i^2}}}{\sqrt{n}-1}} , where
#' \eqn{n} is the length of \code{x}.
#' 
#' The sparseness is a real number in \eqn{[0,1]}. It is equal to 1 if and only
#' if \code{x} contains a single nonzero component, and is equal to 0 if and
#' only if all components of \code{x} are equal.  It interpolates smoothly
#' between these two extreme values.  The closer to 1 is the sparseness the
#' sparser is the vector.
#' 
#' The basic definition is for a \code{numeric} vector. The sparseness of a
#' \code{matrix} is the mean sparseness of its column vectors. The sparseness
#' of an object of class \code{NMF}, is the a 2-length vector that contains the
#' sparseness of the basis and mixture coefficient matrices.  }
#' 
#' \item{syntheticNMF}{ Generate a synthetic matrix according to an underlying
#' NMF model.  It can be used to quickly test NMF algorithms.  }
#' 
#' }
#' 
#' @name Utility functions and methods
#' @aliases utils-NMF connectivity connectivity-methods connectivity,NMF-method
#' entropy entropy-methods entropy,factor,factor-method
#' entropy,NMF,factor-method entropy,NMF,factor-method
#' entropy,table,missing-method entropy,NMF,ANY-method evar evar-methods
#' evar,NMF-method extractFeatures extractFeatures,NMF-method featureScore
#' featureScore,NMF-method featureScore,matrix-method hasTrack nmfApply
#' nmfApply,NMF-method plot,NMFfit,missing-method plot,NMFfit-method purity
#' purity-methods purity,factor,factor-method purity,NMF,factor-method
#' purity,table,missing-method purity,NMF,ANY-method randomize rss
#' rss,NMF-method sparseness sparseness-methods sparseness,matrix-method
#' sparseness,NMF-method sparseness,numeric-method syntheticNMF
#' @docType methods
#' @param class A \code{factor} giving a known class membership for each
#' sample.
#' 
#' In methods \code{entropy} and \code{purity}, argument \code{class} is coerce
#' to a factor if necessary.
#' @param format the output format of the extracted features.  Possible values
#' are:
#' 
#' \itemize{
#' 
#' \item \code{list} (default) a list with one element per basis vector, each
#' containing the indices of the basis-specific features.  \item \code{combine}
#' a single integer vector containing the indices of the basis-specific
#' features for ALL the basis.  \item \code{subset} the object \code{object}
#' subset to contain only the basis-specific features.  }
#' 

#' @param FUN the function to be applied: see 'Details'. In the case of
#' functions like \code{+}, \code{%*%}, etc., the function name must be
#' backquoted or quoted.  See \code{link[base]{apply}} for more details.
#' @param MARGIN a vector giving the subscripts which the function will be
#' applied over. \code{1} indicates rows, \code{2}' indicates columns,
#' \code{c(1,2)} indicates rows and columns.  See \code{link[base]{apply}} for
#' more details.
#' @param method Method used to compute the feature scores and selecting the
#' features.
#' 
#' Possible values are: \itemize{
#' 
#' \item \code{kim} (default) to use Kim and Park (2007) scoring schema and
#' feature selection method.  The features are first scored using the function
#' \code{featureScore}.  Then only the features that fulfil both following
#' criteria are retained:
#' 
#' - score greater than \eqn{\hat{\mu} + 3 \hat{\sigma}}, where \eqn{\hat{\mu}}
#' and \eqn{\hat{\sigma}} are the median and the median absolute deviation
#' (MAD) of the scores respectively;
#' 
#' - the maximum contribution to a basis component is greater than the median
#' of all contributions (i.e. of all elements of W)
#' 
#' See \emph{Kim and Park (2007)}.
#' 
#' \item \code{max} where the score is the maximum contribution of each feature
#' to the basis vectors and the selection method is the one described in
#' \emph{Carmona-Saez (2006)}.  Briefly, for each basis vector, the features
#' are first sorted in descending order by their contribution to the basis
#' vector. Then, one selects only the first consecutive features from the
#' sorted list whose highest contribution in the basis matrix is found in the
#' considered basis (see section \emph{References}).  }
#' @param n Number of rows of the synthetic target matrix.
#' @param noise if \code{TRUE}, a random noise is added the target matrix.
#' @param object A \code{matrix} or an object that inherits from class
#' \code{\linkS4class{NMF}} or \code{\linkS4class{NMFfit}} -- depending on the
#' method.
#' @param offset a vector giving the offset to add to the synthetic target
#' matrix.  Its length should be equal to the number of rows \code{n}.
#' @param p Number of columns of the synthetic target matrix. Not used if
#' parameter \code{r} is a vector (see description of argument \code{r}).
#' @param r Underlying factorization rank. If a single \code{numeric} is given,
#' the classes are randomly generated from a multinomial distribution.  If a
#' numerical vector is given, then it should contain the counts in the
#' different classes (i.e integers). In such a case argument \code{p} is not
#' used and the number of columns is forced to be the sum of the counts.
#' @param return.factors If \code{TRUE}, the underlying matrices \code{W} and
#' \code{H} are also returned.
#' @param target the target object estimated by model \code{object}. It can be
#' a \code{matrix} or an \code{ExpressionSet}.
#' @param x
#' 
#' for \code{randomize}: the \code{matrix} or \code{ExpressionSet} object whose
#' entries will be randomised.
#' 
#' for \code{plot}: An object that inherits from class
#' \code{\linkS4class{NMFfit}}.
#' 
#' otherwise: An object that inherits from class \code{\linkS4class{NMF}}.
#' 

#' @param ...  Used to pass extra parameters to subsequent calls: \itemize{
#' \item in \code{nmfApply}: optional arguments to function \code{FUN}.  \item
#' in \code{randomize}: passed to the \code{sample} function.  }
#' @author Renaud Gaujoux
#' @seealso \code{\linkS4class{NMF}, \link[=summary,NMF-method]{summary}}
#' @references
#' 
#' Brunet, J.P. et al. (2004).  \emph{Metagenes and molecular pattern discovery
#' using matrix factorization}.  Proc Natl Acad Sci USA 101(12), 4164--4169.
#' 
#' Kim, H. and Park, H. (2007).  \emph{Sparse non-negative matrix
#' factorizations via alternating non-negativity-constrained least squares for
#' microarray data analysis}.  Bioinformatics 2007; \bold{23(12)}:1495-502.
#' \url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
#' 
#' Hoyer, P. O. (2004).  \emph{Non-negative Matrix Factorization with
#' Sparseness Constraints}.  \emph{Journal of Machine Learning Research} 5
#' (2004) 1457--1469.
#' 
#' Carmona-Saez, Pedro et al. (2006).  \emph{Biclustering of gene expression
#' data by non-smooth non-negative matrix factorization}.  BMC Bioinformatics
#' 7(1), 78.
#' @examples
#' 
#' 
#' 
#' # generate a synthetic dataset with known classes: 50 features, 18 samples (5+5+8)
#' n <- 50; counts <- c(5, 5, 8);
#' V <- syntheticNMF(n, counts, noise=TRUE)
#' \dontrun{aheatmap(V)}
#' 
#' # build the class factor
#' groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))
#' 
#' # perform default NMF
#' res <- nmf(V, 2)
#' res
#' 
#' \dontrun{coefmap(res, class=groups)}
#' \dontrun{basismap(res)}
#' # compute entropy and purity
#' entropy(res, class=groups)
#' purity(res, class=groups)
#' 
#' # perform NMF with the right number of basis components
#' res <- nmf(V, 3)
#' 
#' \dontrun{coefmap(res)}
#' \dontrun{basismap(res, 'features')}
#' entropy(res, class=groups)
#' purity(res, class=groups)
#' 
#' 
NULL





#' Package Specific Option Management
#' 
#' Allow the user to get/set/define package NMF specific options in the same
#' way as with base functions \code{\link[base]{options}} and
#' \code{\link[base]{getOption}}.
#' 

#' 
#' \describe{
#' 
#' \item{nmf.getOption}{ get the value of a single option.  }
#' 
#' \item{nmf.options}{ gets/sets/defines options in the way of base function
#' \code{\link[base]{options}}.  Invoking 'nmf.options()' with no arguments
#' returns a list with the current values of the options.  To access the value
#' of a single option, one should use \code{nmf.getOption("error.track")},
#' e.g., rather than \code{nmf.options("error.track")} which is a \code{list}
#' of length one.  }
#' 
#' \item{nmf.options.reset}{ Reset all \emph{built-in options} to their default
#' values.  Note that only built-in are reset. The options defined by the user
#' during the current session will keep their values.  }
#' 
#' }
#' 
#' @aliases options-NMF nmf.options nmf.getOption nmf.options.reset
#' nmf.options.runtime
#' @param \dots any options can be defined, using 'name = value' or by passing
#' a list of such tagged values.  However, only the ones below are used in
#' package NMF. Further, 'nmf.options('name') == nmf.options()['name']', see
#' the example.
#' @param name a character string holding an option name.
#' @param runtime a boolean used to specify if main interface function
#' \code{\link{nmf}} should store the option into the initial
#' \code{\linkS4class{NMF}} object before performing the computation.
#' @return For \code{nmf.getOption}, the current value set for option
#' \code{name}, or \code{NULL} if the option is unset.
#' 
#' For \code{nmf.options()}, a list of all set options sorted by name.  For
#' \code{options(name)}, a list of length one containing the set value, or
#' \code{NULL} if it is unset.  For uses setting one or more options, a list
#' with the previous values of the options changed (returned invisibly).
#' @section Options set in package NMF:
#' 
#' \describe{
#' 
#' \item{debug}{ logical. Similar to option \code{'verbose'} (see below), but
#' reports more information. } \item{default.algorithm}{ \code{character}. The
#' default NMF algorithm used by the \code{\link{nmf}} method when called
#' without argument \code{method}. } \item{default.seed}{ \code{character}. The
#' default seeding method used by the \code{\link{nmf}} method when called
#' without argument \code{seed}. } \item{error.track}{ logical. Should the
#' estimation error be tracked during the computations? If set to \code{TRUE}
#' then the error track can be plotted using method \code{\link{plot}}.  The
#' step size of the error track is set via option \code{track.interval} (see
#' below). } \item{parallel.backend}{ \code{character} or \code{numeric}.  The
#' default parallel backend used by the \code{\link{nmf}} method when called
#' with argument \code{nrun} greater than 1.
#' 
#' Currently it accepts the following values: \code{'mc'} or a number that
#' specifies the number of cores to use, \code{'seq'} or \code{NULL} to use
#' sequential backend for \code{foreach}, and the empty string \code{''} to
#' completely disable the parallel computation.
#' 
#' } \item{track.interval}{ numeric. The number of iterations performed between
#' two consecutive error points.  For performance reason, this value should be
#' too small, as the computation of the estimation error can be time consuming
#' (Default value is 30). } \item{verbose}{ logical. Should R report extra
#' information about the computations? }
#' 
#' }
#' @author Renaud Gaujoux
#' @seealso \code{\link{options}}
#' @examples
#' 
#' 
#' 	# save all options value 
#' 	op <- nmf.options(); 
#' 	utils::str(op) # op may contain functions.
#' 
#'     nmf.getOption("track.interval") == nmf.options()$track.interval # the latter needs more memory
#'     
#'     x <- rmatrix(50, 10) # create a random target matrix
#'     # or define a synthetic data with a hidden pattern using function syntheticNMF (see ?syntheticNMF) 
#'     \dontrun{x <- syntheticNMF(50, 5, 10, noise=TRUE)}
#'     
#'     # perform default NMF computation
#'     res <- nmf(x, 3)
#'     
#'     # Toogle on verbose mode
#'     nmf.options(verbose = TRUE)    
#'     res <- nmf(x, 3)    
#' 
#' 	# Toogle on debug mode
#'     nmf.options(debug = TRUE)    
#'     res <- nmf(x, 3)
#' 
#'     # set the error track step size, and save previous value
#'     old.o <- nmf.options(track.interval = 5)
#'     old.o
#'     
#'     # check options
#'     utils::str(nmf.options())
#'     # reset to default values
#'     nmf.options.reset()
#'     utils::str(nmf.options())
#' 	
#' 
NULL





#' Plotting Basis Profiles
#' 
#' The function \code{profplot} draw plots of the basis profiles, i.e. the rows
#' of the mixture coefficient matrix of NMF models. A given profile is composed
#' of the contribution of the corresponding basis to each sample.
#' 
#' When using NMF for clustering in particular, one looks for strong
#' associations between the basis and a priori known groups of samples.
#' Plotting the profiles may highlight such patterns.
#' 
#' The function can also be used to compare the profiles from two NMF models or
#' mixture coefficient matrices. In this case, it draws a scatter plot of the
#' paired profiles.
#' 
#' 
#' @name Plotting Basis Profiles
#' @aliases profplot
#' @docType methods
#' @param x a matrix or an NMF object from which is extracted the mixture
#' coefficient matrix. It is extracted from the best fit if \code{x} is the
#' results from multiple NMF runs.
#' @param y a matrix or an NMF object from which is extracted the mixture
#' coefficient matrix.It is extracted from the best fit if \code{y} is the
#' results from multiple NMF runs.
#' @param scale a logical that specifies whether the columns of the matrices
#' should be scaled into proportions (i.e. to sum up to one) before plotting.
#' Default is \code{FALSE}.
#' @param legend a logical that specifies whether drawing the legend or not, or
#' coordinates specifications passed to argument \code{x} of
#' \code{\link{legend}}, that specifies the position of the legend.
#' @param Colv specifies the way the columns of \code{x} are ordered before
#' plotting. It is used only when \code{y} is missing.  It can be: \itemize{
#' \item a single numeric value, specifying the index of a row of \code{x},
#' that is used to order the columns by \code{x[, order(x[abs(Colv),])]}.
#' Decreasing order is specified with a negative index.  \item an integer
#' vector directly specifying the order itself, in which case the columns are
#' ordered by \code{x[, Colv]} \item a factor used to order the columns by
#' \code{x[, order(Colv)]} and as argument \code{annotation} if this latter is
#' missing or not \code{NA}.  \item any other object with a suitable
#' \code{order} method. The columns are by \code{x[, order(Colv)]} }
#' @param labels a character vector containing labels for each sample (i.e.
#' each column of \code{x}). These are used for labelling the x-axis.
#' @param annotation a factor annotating each sample (i.e. each column of
#' \code{x}). If not missing, a coloured raw is plotted under the x-axis and
#' annotates each sample accordingly. If argument \code{Colv} is a factor, then
#' it is used to annotate the plot, unless \code{annotation=NA}.
#' @param ...  graphical parameters passed to matplot.
#' @seealso \code{\link{profcor}}
#' @keywords aplot
#' @examples
#' 
#' 
#' if( interactive() ){
#' 
#' # create a random target matrix
#' v <- rmatrix(50, 10)
#' 
#' # fit a single NMF model
#' res <- nmf(v, 3)
#' profplot(res)
#' 
#' # ordering according to first profile
#' profplot(res, Colv=1) # increasing
#' profplot(res, Colv=-1) # decreasing
#' 
#' # fit a multi-run NMF model
#' res2 <- nmf(v, 3, nrun=5)
#' profplot(res2)
#' 
#' # draw a profile correlation plot: this show how the basis components are 
#' # returned in an unpredictable order 
#' profplot(res, res2)
#' 
#' # looking at all the correlations allow to order the components in a "common" order
#' profcor(res, res2)
#' 
#' }
#' 
#' 
NULL





#' Creates a Random Matrix Using Any Given Distribution Function
#' 
#' This function provides a short-cut to create a random matrix whose entries
#' are drawn any given random distribution, as soon as this one is implemented
#' as a random variate generation function similar to: \code{runif},
#' \code{rnorm}, etc ...
#' 
#' It essentially wraps the following common call:
#' 
#' \code{matrix(dist(nrow*ncol, ...), nrow, ncol)}
#' 
#' into the following shorter call -- that should also be less prone to errors:
#' 
#' \code{rmatrix(nrow, ncol, dist, ...)}
#' 
#' 
#' @aliases rmatrix rmatrix-methods rmatrix,numeric-method
#' rmatrix,matrix-method rmatrix,NMF-method
#' @param x a \code{numeric} value giving the number of rows of the result
#' matrix, a matrix, or an \code{\linkS4class{NMF}} model.
#' @param y a \code{numeric} value giving the number of columns of the result
#' matrix.  If \var{y} is missing, then the created matrix is square (i.e.
#' \code{y=x}).
#' @param dist a random variate generation function (e.g. \code{\link{runif}},
#' \code{\link{rnorm}}, etc...)  from which to draw the matrix entries.  It
#' must be a function whose first parameter is the number of values \code{n} to
#' be drawn, and return a \code{numeric} vector of length \code{n}.
#' @param byrow logical. If \code{FALSE} (the default) the matrix is filled by
#' columns, otherwise the matrix is filled by rows. See \code{\link{matrix}}.
#' @param dimnames A \code{dimnames} attribute for the matrix: \code{NULL} or a
#' \code{list} of length 2 giving the row and column names respectively.  An
#' empty list is treated as \code{NULL}, and a list of length one as row names.
#' The list can be named, and the list names will be used as names for the
#' dimensions. See \code{\link{matrix}}.
#' @param \dots \itemize{ \item If \code{x} is a \code{numeric}: extra
#' parameters passed to function \code{dist}.  \item If \code{x} is a
#' \code{matrix}: extra parameters passed to function the internal call
#' \code{rmatrix(nrow(x), ncol(x), ...)}.  }
#' @return
#' 
#' a matrix whose entries are drawn from distribution \code{dist} and of
#' dimension \code{x} X \code{y} if \var{x} is a \code{numeric}, or of
#' dimension \code{nrow(x)} X \code{ncol(y)} if \var{x} is a \code{matrix} or
#' an \code{NMF} model.
#' @author Renaud Gaujoux
#' @seealso \code{\link{runif}}, \code{\link{rnorm}} or any other
#' \code{\link[=distributions]{rxxx}} similar random variate generation
#' function.
#' @keywords distribution
#' @examples
#' 
#' 
#' ## Generate a random matrix of a given size
#' rmatrix(5, 3)
#' \dontshow{ stopifnot( identical(dim(rmatrix(5, 3)), c(5L,3L)) ) }
#' 
#' ## Generate a random matrix of the same dimension of a template matrix
#' a <- matrix(1, 3, 4)
#' rmatrix(a)
#' \dontshow{ stopifnot( identical(dim(rmatrix(a)), c(3L,4L)) ) }
#' 
#' ## Generate a random matrix of the dimension of the target matrix of an NMF model 
#' a <- nmfModel(2, 10, 5)
#' rmatrix(a)
#' \dontshow{ stopifnot( identical(dim(rmatrix(a)), c(10L,5L)) ) }
#' 
#' ## Specificy the distribution to use
#' 
#' # the default is uniform
#' a <- rmatrix(1000, 50)
#' \dontrun{ hist(a) }
#' 
#' # use normal ditribution
#' a <- rmatrix(1000, 50, rnorm)
#' \dontrun{ hist(a) }
#' 
#' # extra arguments can be passed to the random variate generation function 
#' a <- rmatrix(1000, 50, rnorm, mean=2, sd=0.5)
#' \dontrun{ hist(a) }
#' 
#' 
NULL





#' Get the Random Number Generator State
#' 
#' The function \code{getRNG} is a S4 generic function that returns either the
#' current RNG settings if called with no argument, or the RNG settings used to
#' generate an object.  In the case of results from multiple NMF runs, it
#' returns the RNG settings used to compute the best fit.
#' 
#' The function \code{getRNG1} is a S4 generic function that returns the RNG
#' settings used for the first of the multiple NMF runs performed to generate
#' an \code{\linkS4class{NMFfitX}} object.  This is useful to reproduce or
#' compare results from different multiple runs.
#' 
#' The function \code{RNGdigest} returns a hash that uniquely identifies given
#' RNG settings.  %The hash can be -- and is -- used to compare two
#' \code{\link{rstream}} %objects and test if they represent the same random
#' streams, in the same state.  The hash is computed based on the result of
#' \code{getRNG}.
#' 
#' The function \code{rng.equal} tests whether two -- embedded -- RNG objects
#' represent the same random stream, in the same state. It does so by comparing
#' the hashes returned by \code{RNGdigest}.
#' 
#' The function \code{rng1.equal} tests whether two results from multiple NMF
#' runs used identical initial RNG settings.
#' 
#' The current methods for \code{getRNG} extracts the RNG settings from slot
#' \code{rng} for S4 objects, or from element \code{rng} or \code{noise$rng}
#' for a \code{list} (wich includes S3 objects), or from attribute \code{rng}
#' otherwise.
#' 
#' @name Utilities for Random Number Generators (RNG)
#' @aliases RNG-NMF getRNG getRNG-methods getRNG,missing-method
#' getRNG,list-method getRNG,ANY-method getRNG,NMFfitXn-method
#' getRNG,integer-method getRNG1 getRNG1-methods getRNG1,NMFfit-method
#' getRNG1,NMFfitX-method getRNG1,NMFfitX1-method getRNG1,NMFfitXn-method
#' RNGdigest rng.equal rng1.equal
#' @docType methods
#' @param object either \code{missing}, a \code{list}, or any object, usually
#' with relevant stored RNG information, such as \code{\linkS4class{NMFfit}}
#' objects. If missing, the returned value is based on the current RNG
#' settings.
#' 
#' For \code{getRNG1}, an \code{NMFfitX} object.
#' @param x any object handled by \code{getRNG}. For \code{rng1.equal},
#' \code{x} must be an \code{\linkS4class{NMFfitX}} object, such as returned by
#' multiple NMF runs.
#' @param y any object handled by \code{getRNG}. If \code{y} is missing, then
#' the current RNG settings are used. For \code{rng1.equal}, \code{y} must be
#' an \code{\linkS4class{NMFfitX}} object, such as returned by multiple NMF
#' runs.
#' @return
#' 
#' \code{getRNG} returns an integer vector (see \code{\link{.Random.seed}})
#' %For \code{getRNG}, the class of the returned value depends on the type of
#' the %extracted object (see details). %For \code{NMFfit} objects, the result
#' is an \code{\linkS4class{rstream}} object %that contains the RNG settings as
#' before the computation starts, i.e. before the %seeding step.
#' 
#' \code{RNGdigest} returns the RNG hash as a single character string.
#' 
#' \code{rng.equal} and \code{rng1.equal} returns \code{TRUE} if the two --
#' embedded -- RNG objects are identical, \code{FALSE} otherwise. See details
#' for the differences between these two functions.
#' @seealso \code{\link{.Random.seed}}, \code{\link{set.seed}}.
#' @keywords methods
#' @examples
#' 
#' 
#' # current random seed (by default it is a 626-length numeric vector)
#' getRNG()
#' 
#' # fit an NMF model on a random target matrix
#' V <- rmatrix(100,20)
#' s <- getRNG()
#' res <- nmf(V, 3)
#' 
#' # the random seed changed since by default nmf use a randomly generating starting point 
#' rng.equal(s)
#' # but the starting RNG settings are stored in the object and can be accessed by getRNG
#' getRNG(res)
#' # the setting used are the one in use before running NMF  
#' rng.equal(s, res)
#' \dontshow{ stopifnot(rng.equal(s, res)) }
#' 
#' # show the digest version
#' RNGdigest(res)
#' 
#' # For multiple runs, the RNG settings used for the first run is also stored
#' res <- nmf(V, 3, nrun=4)
#' # RNG used for the best fit
#' getRNG(res)
#' # RNG used for the first fit
#' getRNG1(res)
#' # they may differ if the best fit is not the first one
#' rng.equal(res, getRNG1(res))
#' 
#' 
#' 
NULL





#' Undocumented Internal Methods and Functions for RNG Handling
#' 
#' The functions documented here used internally by the package NMF to manage
#' RNG settings.
#' 
#' 
#' @aliases RNG-internals initialize,rstream2-method
#' initialize,rstream.user-method RNGinfo RNGlib RNGlibs RNGprovider
#' RNGrecovery RNGseed RNGseq setRNG setRNG-methods setRNG,ANY-method
#' setRNG,character-method setRNG,numeric-method
#' @author Renaud Gaujoux
#' @seealso \code{\link[NMF]{getRNG-methods}} and the package's vignette
#' @keywords internal
NULL





#' Generates Random NMF Models
#' 
#' Generates NMF models with random values drawn from a uniform distribution.
#' It returns an NMF model with basis and mixture coefficient matrices filled
#' with random values. The main purpose of the function \code{rnmf} is to
#' provide a common interface to generate random seeds used by the
#' \code{\link{nmf}} function.
#' 
#' If necessary, extensions of the standard NMF model or custom models must
#' define a method "rnmf,<NMF.MODEL.CLASS>,numeric" for initializing their
#' specific slots other than the basis and mixture coefficient matrices. In
#' order to benefit from the complete built-in interface, the overloading
#' methods should call the generic version using function
#' \code{\link{callNextMethod}}. See for example the method
#' \code{\link[=rnmf,NMFOffset,numeric-method]{rnmf,NMFOffset,numeric}}:
#' \code{showMethods(rnmf, class='NMFOffset', include=TRUE))}.
#' 
#' 
#' @name rnmf-methods: Generating random NMF models
#' @aliases rnmf rnmf-methods rnmf,ANY,matrix-method
#' rnmf,ANY,ExpressionSet-method rnmf,NMF,missing-method
#' rnmf,NMF,numeric-method rnmf,numeric,numeric-method
#' @docType methods
#' @return An NMF model, i.e. an object that inherits from class
#' \code{\linkS4class{NMF}}.
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"ANY\", target = \"matrix\")")}{ Generates an NMF
#' model that is compatible with the dimensions of the matrix \code{target}.
#' The entries are uniformly drawn between \code{0} and \code{max(target)}. If
#' \code{x} is an NMF model, the returned value is an object of the same class
#' as \code{x}. If \code{x} is a single numeric value then the returned value
#' is a standard NMF model with rank \code{x} (see
#' \code{\linkS4class{NMFstd}}). }
#' 
#' \item{list("signature(x = \"ANY\", target = \"ExpressionSet\")")}{ Generates
#' an NMF model that is compatible with the dimensions of the target expression
#' matrix \code{exprs{target}} (see
#' \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}).
#' 
#' This is a shortcut for \code{rnmf(x, exprs(target), ...)}.
#' 
#' }
#' 
#' \item{list("signature(x = \"NMF\", target = \"missing\", ...)")}{ Generates
#' an NMF model of the same class and dimension as \code{x}. The entries are
#' uniformly drawn between \code{0} and \code{max} (optionally specified in
#' \code{...}) that defaults to 1. }
#' 
#' \item{list("signature(x = \"NMF\", target = \"numeric\", ncol=NULL,
#' keep.names=TRUE, ...)")}{ Generates an NMF model of the same class as
#' \code{x}, compatible with the dimensions \code{target}, that can be a single
#' or 2-length numeric vector, to specify a square or rectangular target matrix
#' respectively. The second dimension can also be passed via argument
#' \code{ncol}. By default the dimnames of \code{x} are set on the returned NMF
#' model.  This behaviour is disabled with argument \code{keep.names=FALSE}.
#' See \code{\link{nmfModel}}. The entries are uniformly drawn between \code{0}
#' and \code{max} (optionally specified in \code{...}) that defaults to 1. }
#' 
#' \item{list("signature(x = \"numeric\", target = \"numeric\", ncol=NULL,
#' ...)")}{ This is a shortcut for \code{rnmf(nmfModel(x, target, ncol)),
#' ...)}.
#' 
#' It generates a standard NMF model compatible with the dimensions
#' \code{target}, that can be a single or 2-length numeric vector, to specify a
#' square or rectangular target matrix respectively. See
#' \code{\link{nmfModel}}. The entries are uniformly drawn between \code{0} and
#' \code{max} (optionally specified in \code{...}) that defaults to 1. }
#' 
#' }%end_describe
#' @author Renaud Gaujoux
#' @seealso \code{\linkS4class{NMF}}, \code{\linkS4class{NMFOffset}},
#' \code{\link{rmatrix}}
#' @keywords methods
#' @examples
#' 
#' 
#' # generate a random NMF model with rank 3 that fits a 100x20 matrix  
#' rnmf(3, 100, 20)
#' \dontshow{ stopifnot( identical(dim(rnmf(3, 100, 20)), c(100L,20L,3L)) ) }
#' # generate a random NMF model with rank 3 that fits a 100x100 matrix
#' rnmf(3, 100)
#' \dontshow{ stopifnot( identical(dim(rnmf(3, 100)), c(100L,100L,3L)) ) }
#' 
#' # generate a random NMF model based on an existing NMF model
#' a <- nmfModel(3, 100, 20)
#' b <- rnmf(a)
#' \dontshow{ stopifnot( !nmf.equal(a,b) ) }
#' 
#' # generate a random NMF model of the same class as an existing NMF model but
#' # of different dimension
#' a <- nmfModel(3, 100, 20, model='NMFns')
#' c <- rnmf(a, 50, 10)
#' \dontshow{ stopifnot( identical(dim(c), c(50L,10L,3L)) ) }
#' \dontshow{ stopifnot( is(c, 'NMFns') ) }
#' 
#' 
#' 
NULL





#' Methods for the Interface Defined in Package stats
#' 
#' The package NMF defines methods \code{deviance}, \code{predict} and
#' \code{residuals} for common NMF objects returned by the function
#' \code{\link{nmf}}.
#' 

#' 
#' \describe{
#' 
#' \item{deviance}{ returns the final approximation error between an NMF model
#' and its target matrix, as the value of the minimized objective function.  In
#' the case of multiple NMF runs, it returns the approximation error of the
#' best fit across the runs.
#' 
#' If not computed by the NMF algorithm itself, the value is automatically
#' computed at the end of the fitting process by the function
#' \code{\link{nmf}}, using the objective function associated with the NMF
#' algorithm.}
#' 
#' \item{predict}{ Computes the dominant basis component for each sample (resp.
#' feature) based on its associated entries in the mixture coefficient matrix
#' (i.e in \eqn{H}) (resp. basis component matrix (i.e in \eqn{W})).
#' 
#' When \code{what='samples'} the computation is performed on the mixture
#' coefficient matrix, or on the transposed basis matrix when
#' \code{what='features'}.
#' 
#' For each column, the dominant basis component is computed as the row index
#' for which the entry is the maximum within the column.
#' 
#' If argument \code{prob=FALSE} (default), the result is a \code{factor}.
#' Otherwise it returns a list with two elements: element \code{predict}
#' contains the computed indexes ( as a \code{factor}) and element \code{prob}
#' contains the vector of the associated probabilities, that is the relative
#' contribution of the maximum entry within each column.
#' 
#' }
#' 
#' \item{residuals}{ returns the -- final -- residuals between the target
#' matrix and the NMF result \code{object}. They are computed using the
#' objective function associated to the NMF algorithm that returned
#' \code{object}.  When called with \code{track=TRUE}, the whole residuals
#' track is returned, if available. Note that method \code{\link{nmf}} does not
#' compute the residuals track, unless explicitly required.
#' 
#' It is a S4 methods defined for the associated generic functions from package
#' \code{stats} (See \link[stats]{residuals}).
#' 
#' Note that stricly speaking, this method does not fulfill its contract as
#' defined by the package \code{stats}, but rather acts as function
#' \code{deviance}.  The will be changed in a later release to make it
#' consistent with its original purpose.  } }
#' 
#' @name Methods from stats: predict, deviance, residuals
#' @aliases stats-NMF deviance,NMFfit-method deviance,NMFfitXn-method predict
#' predict-methods predict,matrix-method predict,NMF-method
#' predict,NMFfitX-method residuals residuals,NMFfit-method
#' @docType methods
#' @param object A \code{matrix} or an object that inherits from class
#' \code{\linkS4class{NMF}} or \code{\linkS4class{NMFfit}} -- depending on the
#' method.
#' @param prob Should the probability associated with each cluster prediction
#' be computed and returned.
#' @param track if \code{TRUE}, the whole residuals track is returned.
#' Otherwise only the last residuals computed is returned.
#' @param what Specifies on which matrix the computation should be based:
#' \itemize{ \item \code{'columns', 'samples'}: mixture coefficients \item
#' \code{'rows', 'features'}: basis components \item \code{'consensus',
#' 'cmap'}: consensus matrix } See section \emph{Details} for more information.
#' @param ... argument passed to \code{predict,NMF} if \code{what} is not
#' \code{'consensus'}.
#' @seealso \code{\linkS4class{NMF}, \link[=summary,NMF-method]{summary}}
#' @examples
#' 
#' 
#' 
#' # generate random data matrix
#' X <- rmatrix(50, 20)
#' 
#' # run NMF
#' res <- nmf(X, 2)
#' res
#' 
#' # get the predicted column clusters
#' predict(res)
#' # get the predicted row clusters
#' predict(res, 'rows')
#' 
#' # get the final approximation error
#' deviance(res)
#' # this is the same as the residuals and no track is computed by default
#' residuals(res)
#' residuals(res, track=TRUE)
#' # compute the track
#' res <- nmf(X, 2, .options='t')
#' residuals(res, track=TRUE)
#' 
#' # perform multiple NMF runs
#' res <- nmf(X, 2, nrun=3)
#' 
#' # get the predicted consensus column clusters
#' predict(res, 'consensus')
#' 
#' 
NULL



