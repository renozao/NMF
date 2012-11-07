\name{NMF-package}
\Rdversion{1.1}
\alias{NMF-package}
\alias{NMF}
\docType{package}
\title{NMF Package Overview}
\description{
The NMF package provides methods to perform Nonnegative Matrix Factorization (NMF)
, as well as a framework to develop and test new NMF algorithms.

A number of standard algorithms and seeding methods are implemented.
Tuned visualisation and post-analysis methods help in the evaluation of the 
algorithms' performances or in the interpretation of the results.

}

\author{ Renaud Gaujoux \email{renaud@cbio.uct.ac.za} }

\references{

	If you use the 'NMF' package in publications cite:

	Renaud Gaujoux, Cathal Seoighe (2010).\cr
	\emph{A flexible R package for nonnegative matrix factorization}.\cr
	BMC Bioinformatics 2010, \bold{11}:367\cr
	\url{http://www.biomedcentral.com/1471-2105/11/367}

	For pointers on further readings, please see references therein.


	Definition of Nonnegative Matrix Factorization in its modern formulation:
	
	Lee D.D. and Seung H.S. (1999).
	Learning the parts of objects by non-negative matrix factorization. 
	\emph{Nature}, \bold{401}, 788--791.

	Historical first definition and algorithms: 
		
	Paatero, P., Tapper, U. (1994).
	Positive matrix factorization: A non-negative factor model with optimal utilization of error estimates of data values.
	\emph{Environmetrics}, \bold{2}, 111--126
	, doi:10.1002/env.3170050203.
	
}

\keyword{ package }
\seealso{
\code{\link{NMF-class}}, \code{\link{nmf}}, \code{\link[Biobase]{Biobase}}
}
\examples{

# run default NMF algorithm on a random matrix
V <- rmatrix(500, 20)
res <- nmf(V, 3)
res

# compute some quality measures
summary(res)

# Visualize the results as heatmaps
\dontrun{coefmap(res)} # mixture coefficients
\dontrun{basismap(res)} # basis vectors

# run default NMF algorithm on a random matrix with actual patterns
set.seed(123456)
V <- syntheticNMF(500, 3, 20, noise=TRUE)
res <- nmf(V, 3)
res

# compute some quality measures
summary(res)

# Visualize the results as heatmaps
\dontrun{coefmap(res)} # mixture coefficients
\dontrun{basismap(res)} # basis vectors


}
\keyword{package}
