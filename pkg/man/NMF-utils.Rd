\name{NMF-utils}
\Rdversion{1.1}
\docType{methods}
\alias{NMF-utils}

\alias{connectivity}
\alias{connectivity-methods}
\alias{connectivity,NMF-method}

\alias{entropy}
\alias{entropy-methods}
\alias{entropy,factor,factor-method}
\alias{entropy,NMF,factor-method}
\alias{entropy,NMF,factor-method}
\alias{entropy,table,missing-method}
\alias{entropy,NMF,ANY-method}

\alias{evar}
\alias{evar-methods}
\alias{evar,NMF-method}

\alias{extractFeatures}
\alias{extractFeatures,NMF-method}
\alias{featureScore}
\alias{featureScore,NMF-method}
\alias{featureScore,matrix-method}

\alias{metaHeatmap}
\alias{metaHeatmap-methods}
\alias{metaHeatmap,matrix-method}
\alias{metaHeatmap,NMF-method}

\alias{nmfApply}
\alias{nmfApply,NMF-method}

\alias{plot,NMFfit,missing-method}
\alias{plot,NMFfit-method}

%\alias{predict.nmf}
\alias{predict}
\alias{predict-methods}
\alias{predict,NMF-method}
\alias{predict,matrix-method}

\alias{purity}
\alias{purity-methods}
\alias{purity,factor,factor-method}
\alias{purity,NMF,factor-method}
\alias{purity,table,missing-method}
\alias{purity,NMF,ANY-method}

\alias{randomize}

\alias{residuals}
\alias{residuals,NMFfit-method}
\alias{rss}
\alias{rss,NMF-method}

\alias{sparseness}
\alias{sparseness-methods}
\alias{sparseness,matrix-method}
\alias{sparseness,NMF-method}
\alias{sparseness,numeric-method}

\alias{syntheticNMF}

\title{ Class and Utility Methods for NMF objects}

\description{
Define generic interface methods for class \code{\linkS4class{NMF}}, which is 
the base -- virtual -- class of the results from any NMF algorithms implemented 
within package NMF's framework.

}

\usage{

\S4method{connectivity}{NMF}(x, ...)

\S4method{entropy}{NMF,factor}(x, class, ...)

\S4method{evar}{NMF}(object, target)

\S4method{residuals}{NMFfit}(object, track=FALSE, ...)

rss(object, ...)
\S4method{rss}{NMF}(object, target)

\S4method{featureScore}{NMF}(object, method=c('kim', 'max'))

\S4method{extractFeatures}{NMF}(object, method=c('kim', 'max')
								, format=c('list', 'combine', 'subset'))

\S4method{metaHeatmap}{NMF}(object, what=c('samples', 'features'), filter=FALSE, ...)

\S4method{nmfApply}{NMF}(object, MARGIN, FUN, ...)

\S4method{plot}{NMFfit}(x, ...)

\S4method{predict}{NMF}(object, what = c('samples', 'features'), prob=FALSE)

\S4method{purity}{NMF,factor}(x, class, ...)

randomize(x, ...)

\S4method{sparseness}{NMF}(x)

syntheticNMF(n, r, p, offset=NULL, noise=FALSE, return.factors=FALSE)

}

\arguments{

	\item{class}{ A \code{factor} giving a known class membership for each sample.
	
	In methods \code{entropy} and \code{purity}, argument \code{class} is coerce to a 
	factor if necessary.
	}
	
	\item{filter}{ if \code{TRUE}, only the features that are basis-specific are used. 
	Those features are those returned by function \code{extractFeatures}. }
	
	\item{format}{ the output format of the extracted features.
	Possible values are:
	
	\itemize{
	
	\item \code{list} (default) a list with one element per basis vector, each containing
	the indices of the basis-specific features. 
	\item \code{combine} a single integer vector containing the indices of the basis-specific 
	features for ALL the basis.
	\item \code{subset} the object \code{object} subset to contain only the basis-specific 
	features.
	}
	
	}
	
	\item{FUN}{ the function to be applied: see 'Details'. In the case of
	functions like \code{+}, \code{\%*\%}, etc., the function name must be
	backquoted or quoted.
	See \code{link[base]{apply}} for more details.}

	\item{MARGIN}{ a vector giving the subscripts which the function will be
		applied over. \code{1} indicates rows, \code{2}' indicates columns,
		\code{c(1,2)} indicates rows and columns.
		See \code{link[base]{apply}} for more details.
		}

	
	\item{method}{ Method used to compute the feature scores and selecting the features.
	
	Possible values are:
		\itemize{
		
		\item \code{kim} (default) to use Kim and Park (2007) scoring schema and feature selection 
		method.
		The features are first scored using the function \code{featureScore}. 
		Then only the features that fulfil both following criteria are retained:
	 	
	 	- score greater than \eqn{\hat{\mu} + 3 \hat{\sigma}}, where \eqn{\hat{\mu}} 
	 	and \eqn{\hat{\sigma}} are the median and the median absolute deviation (MAD) 
	 	of the scores respectively;
	 	
	 	- the maximum contribution to a basis component is greater than the median 
	 	of all contributions (i.e. of all elements of W)
	
		See \emph{Kim and Park (2007)}.
					
		\item \code{max} where the score is the maximum contribution of each feature 
		to the basis vectors and the selection method is the one described in \emph{Carmona-Saez (2006)}.
		Briefly, for each basis vector, the features are first sorted in descending order by their 
		contribution to the basis vector. Then, one selects only the first consecutive features from 
		the sorted list whose highest contribution in the basis matrix is found in 
		the considered basis (see section \emph{References}).
		}
	}
	
	\item{n}{ Number of rows of the synthetic target matrix. }	
	
	\item{noise}{ if \code{TRUE}, a random noise is added the target matrix. }
	
	\item{object}{ A \code{matrix} or an object that inherits from class 
	\code{\linkS4class{NMF}} or \code{\linkS4class{NMFfit}} -- depending on the method. }
	
	\item{offset}{ a vector giving the offset to add to the synthetic target matrix. 
	Its length should be equal to the number of rows \code{n}.}
	
	\item{prob}{ Should the probability associated with each cluster prediction be 
	computed and returned.
	}
	
	\item{p}{ Number of columns of the synthetic target matrix. Not used if parameter 
	\code{r} is a vector (see description of argument \code{r}).}
	
%	\item{palette}{ A character ctring that defines a color palette used 
%	internally with function \code{\link[RColorBrewer]{brewer.pal}}. }
	
	\item{r}{ Underlying factorization rank. If a single \code{numeric} is given, 
	the classes are randomly generated from a multinomial distribution. 
	If a numerical vector is given, then it should contain the counts in the different 
	classes (i.e integers). In such a case argument \code{p} is not used and the number of columns 
	is forced to be the sum of the counts.}
	
	\item{return.factors}{ If \code{TRUE}, the underlying matrices \code{W} and 
	\code{H} are also returned.}
	
%	\item{scale.to.one}{ if \code{TRUE}, the columns are scaled to sum up to 1. }
	
	\item{target}{ the target object estimated by model \code{object}. It can be 
	a \code{matrix} or an \code{ExpressionSet}.
	}
	
	\item{track}{ if \code{TRUE}, the whole residuals track is returned. 
	Otherwise only the last residuals computed is returned.}
	
	\item{what}{ Specifies on which matrix (basis components or mixture coefficients) 
	the computation should be performed.}
	
	\item{x}{ 
	
		for \code{randomize}: the \code{matrix} or \code{ExpressionSet} object whose 
		entries will be randomised.
		
		for \code{plot}: An object that inherits from class \code{\linkS4class{NMFfit}}.
		
		otherwise: An object that inherits from class \code{\linkS4class{NMF}}. 
				
		}	
	
	\item{...}{ Used to pass extra parameters to subsequent calls:
		\itemize{
			\item in \code{metaHeatmap}: Graphical parameters passed to function \code{\link[gplots]{heatmap.2}}
			\item in \code{nmfApply}: optional arguments to function \code{FUN}.
			\item in \code{randomize}: passed to the \code{sample} function.
			\item in \code{residuals}: not used.
		}
	}
}


\details{

\describe{
	
	\item{connectivity}{
	Computes the connectivity matrix for the samples based on their mixture coefficients.
	
 	The connectivity matrix of a clustering is a matrix \eqn{C} containing 
 	only 0 or 1 entries such that:
 	\deqn{C_{ij} = \left\{\begin{array}{l}1\mbox{ if sample }i\mbox{ belongs to the same cluster as sample }j\\0\mbox{ otherwise}\end{array}\right..}{%
 	C_{ij} = 1 if sample i belongs to the same cluster as sample j, 0 otherwise}
	
	}
		
	\item{entropy}{
	The entropy is a measure of performance of a clustering method, in 
	recovering classes defined by factor a priori known	(i.e. one knows the 
	true class labels).
	Suppose we are given \eqn{l} categories, while the clustering method 
	generates \eqn{k} clusters. Entropy is given by:
 	\deqn{Entropy = - \frac{1}{n \log_2 l} \sum_{q=1}^k \sum_{j=1}^l n_q^j \log_2 \frac{n_q^j}{n_q}}
	, where:
	
 - \eqn{n} is the total number of samples;
 
 - \eqn{n} is the total number of samples in cluster \eqn{q};
 
 - \eqn{n_q^j} is the number of samples in cluster \eqn{q} that belongs to 
 original class \eqn{j} (\eqn{1 \leq j \leq l}).

	 The smaller the entropy, the better the clustering performance.
	 
	 See \emph{Kim and Park (2007)}.
	}

	\item{evar}{
	Computes the explained variance of the NMF model \code{object}.
	
 	For a target \eqn{V} It is defined as:
 	\deqn{evar = 1 - \frac{RSS}{\sum_{i,j} v_{ij}^2}}{ evar = 1 - RSS/sum v_{ij}^2},
 	
 	where RSS is the residual sum of squares.
	
	It is usefull to compare the performance of different models and their ability 
	to accurately reproduce the original target matrix.
	Note that a possible caveat is that some methods explicitly aim at minimizing 
	the RSS (i.e. maximizing the explained variance), while others do not.
	}
	
	\item{extractFeatures}{
	Identify the most basis-specific features, using different methods.
	See details of argument \code{method}.
	}
	
	\item{featureScore}{
	Computes the feature scores as suggested in \emph{Kim and Park (2007)}.

	The score for feature \eqn{i} is defined as:
 	\deqn{S_i = 1 + \frac{1}{\log_2 k} \sum_{q=1}^k p(i,q) \log_2 p(i,q),}
 	where \eqn{p(i,q)} is the probability that the \eqn{i}-th feature 
 	contributes to basis \eqn{q}:
 	\deqn{p(i,q) = \frac{W(i,q)}{\sum_{r=1}^k W(i,r)} }

 	The feature scores are real values within the range [0,1]. 
 	The higher the feature score the more basis-specific the corresponding feature.

	}

	\item{metaHeatmap}{ Produces a heatmap of the basis components or mixture 
	coefficients using a \code{heatmap}-like custom function, 
	with parameters tuned for displaying NMF results.
	
	The function used to draw the heatmap is a mixture of the function \code{heatmap.2} 
	from the \code{gplots} package, and the function \code{heatmap.plus}
	 from the \code{heatmap.plus} package. It allows to add extra annotation rows 
	 using the \code{ColSideColor} argument.
	 See \code{\link[gplots]{heatmap.2}} and \code{\link[heatmap.plus]{heatmap.plus}}.
	}

	\item{nmfApply}{ \code{apply}-like method for objects of class \code{NMF}.
	
	When argument \code{MARGIN=1}, it calls the base method \code{apply} to apply 
	function \code{FUN} to the \emph{rows} of the basis component matrix.
	
	When \code{MARGIN=2}, it calls the base method \code{apply} to apply 
	function \code{FUN} on the \emph{columns} of the mixture coefficient matrix.
	
	See \code{\link[base]{apply}} for more details on the output format.
	}

	\item{plot}{ plots the residuals track of the run that computed object \code{x}. 
    See function \code{\link{nmf}} for details on how to enable the tracking of residuals.}
    
	\item{predict}{ Computes the dominant basis component for each sample (resp. feature) 
	based on its associated entries in the mixture coefficient matrix (i.e in \eqn{H}) 
	(resp. basis component matrix (i.e in \eqn{W})).
	
	When \code{what='samples'} the computation is performed on the mixture 
	coefficient matrix, or on the transposed basis matrix when \code{what='features'}.
	
	For each column, the dominant basis component is computed as the row index 
	 for which the entry is the maximum within the column.
	 
	If argument \code{prob=FALSE} (default), the result is a \code{factor}. 
	Otherwise it returns a list with two elements: element \code{predict} contains 
	the computed indexes ( as a \code{factor}) and element \code{prob} contains 
	the vector of the associated probabilities, that is the relative contribution 
	of the maximum entry within each column.
	
	}

	
	\item{purity}{
	Computes the purity of a clustering given a known factor.

	The purity is a measure of performance of a clustering method, in 
	recovering the classes defined by a factor a priori known (i.e. one knows 
	the true class labels).
 	Suppose we are given \eqn{l} categories, while the clustering method 
 	generates \eqn{k} clusters. Purity is given by:
 	\deqn{Purity = \frac{1}{n} \sum_{q=1}^k \max_{1 \leq j \leq l} n_q^j} 	
	, where:
	
 	- \eqn{n} is the total number of samples;
 	
 	- \eqn{n_q^j} is the number of samples in cluster \eqn{q} that belongs to 
 	original class \eqn{j} (\eqn{1 \leq j \leq l}).

	The purity is therefore a real number in \eqn{[0,1]}. 
	The larger the purity, the better the clustering performance.
	
	See \emph{Kim and Park (2007)}.
	}
	
	\item{randomize}{ permute the row entries within each column of \code{x}, 
	using a different permutation for each column.  
	
	The extra arguments in \code{...} are passed to the \code{sample} function, 
	and will be used for each column.
	
	The result is a \code{matrix} of the same dimension as \code{x} (or \code{exprs(x)} 
	in the case \code{x} is an \code{ExpressionSet} object.).
	}
	
	\item{residuals}{ returns the -- final -- residuals between the target matrix and the 
	NMF result \code{object}. They are computed using the objective function 
	associated to the NMF algorithm that returned \code{object}.
	When called with \code{track=TRUE}, the whole residuals track is returned, 
	if available. Note that method \code{nmf} does not compute the residuals track, 
	unless explicitly required.
	
	It is a S4 methods defined for the associated generic functions from package 
	\code{stats} (See \link[stats]{residuals})
	
	See \code{\link{nmf}} and \code{\linkS4class{NMFfit}}.
	}
	
	
	\item{rss}{
	returns the Residual Sum of Squares (RSS) between the target object \code{target}
	 and its estimation by the \code{object}. \emph{Hutchins et al. (2008)} used 
	 the variation of the RSS in combination with \emph{Lee and Seung}'s algorithm
    to estimate the correct number of basis vectors. The optimal rank is chosen 
    where the graph of the RSS first shows an inflexion 
    point. See references.
    
    Note that this way of estimation may not be suitable for all models. Indeed, 
    if the NMF optimization problem is not based on the Frobenius norm, the RSS is 
    not directly linked to the quality of approximation of the NMF model.
    
	}
	
	\item{sparseness}{
	Generic mathod that computes the sparseness of an object as defined in \emph{Hoyer (2004)}.
	
	 This sparseness measure quantifies how much energy of a vector is packed into only few components.
	 It is defined by:
	 \deqn{Sparseness(x) = \frac{\sqrt{n} - \frac{\sum |x_i|}{\sqrt{\sum x_i^2}}}{\sqrt{n}-1}}
	, where \eqn{n} is the length of \code{x}.
	
	 The sparseness is a real number in \eqn{[0,1]}. It is equal to 1 if and only if \code{x} contains 
	 a single nonzero component, and is equal to 0 if and only if all components of \code{x} are equal.
	 It interpolates smoothly between these two extreme values. 
	 The closer to 1 is the sparseness the sparser is the vector.
	 
	 The basic definition is for a \code{numeric} vector. The sparseness of a \code{matrix} is the mean 
	 sparseness of its column vectors. The sparseness of an object of class \code{NMF}, is the a 2-length
	 vector that contains the sparseness of the basis and mixture coefficient matrices.
	}
	
	\item{syntheticNMF}{
		Generate a synthetic matrix according to an underlying NMF model. 
		It can be used to quickly test NMF algorithms.	
	}

}

}

\seealso{	
	\code{\linkS4class{NMF}, \link[=summary,NMF-method]{summary}}  
}

\references{

 	\emph{Metagenes and molecular pattern discovery using matrix factorization}
	Brunet, J.~P., Tamayo, P., Golub, T.~R., and Mesirov, J.~P. (2004)
	Proc Natl Acad Sci U S A
	101(12), 4164--4169.

	\emph{Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis}
	Kim, H. & Park, H. (2007)
	Bioinformatics. 
	\url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
	
	\emph{Non-negative Matrix Factorization with Sparseness Constraints}
	Hoyer, P. O. (2004)
	\emph{Journal of Machine Learning Research} 5 (2004) 1457--1469
	
	\emph{Biclustering of gene expression data by non-smooth non-negative matrix factorization}
	Carmona-Saez, Pedro and Pascual-Marqui, Roberto and Tirado, F and Carazo, Jose and Pascual-Montano, Alberto (2006)
	BMC Bioinformatics
	7(1), 78
}

\author{ Renaud Gaujoux \email{renaud@cbio.uct.ac.za} }

\examples{


# generate a synthetic dataset with known classes: 50 features, 18 samples (5+5+8)
n <- 50; counts <- c(5, 5, 8);
V <- syntheticNMF(n, counts, noise=TRUE)
\dontrun{metaHeatmap(V)}

# build the class factor
groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))

# perform default NMF
res <- nmf(V, 2)
res

\dontrun{metaHeatmap(res, class=groups)}
\dontrun{metaHeatmap(res, 'features')}
# see the predicted clusters of samples
predict(res)
# compute entropy and purity
entropy(res, class=groups)
purity(res, class=groups)

# perform NMF with the right number of basis components
res <- nmf(V, 3)

\dontrun{metaHeatmap(res)}
\dontrun{metaHeatmap(res, 'features')}
entropy(res, class=groups)
purity(res, class=groups)

}

