\name{nmfEstimateRank: Rank estimation}
\Rdversion{1.1}
\alias{nmfEstimateRank}
\alias{plot.NMF.rank}

\title{Estimate optimal rank for Nonnegative Matrix Factorization (NMF) models }
\description{
A critical parameter in NMF algorithms is the factorization rank \eqn{r}. 
It defines the number of basis effects used to approximate the target matrix.
Function \code{nmfEstimateRank} helps in choosing an optimal rank by implementing 
simple approaches proposed in the literature.

}
\usage{
nmfEstimateRank(x, range, method = nmf.getOption("default.algorithm"), nrun = 30, verbose=FALSE, stop=FALSE, ...)

plot.NMF.rank(x, what = c('all', 'cophenetic', 'rss', 'residuals'
						, 'dispersion', 'evar', 'sparseness'
						, 'sparseness.basis', 'sparseness.coef')
						, ref=NULL, na.rm=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%	\item{conf.interval}{
%	  a single \code{logical} specifying if confidence intervals should be estimated 
%	  for all the computed consensus measures.
%	  For each rank in \code{range}, the confidence intervals are estimated by 
%	  bootstrap, resampling 5\code{nrun} times with replacement across the \code{nrun}
%	  runs.
%	}
	
	\item{method}{
	  A single NMF algorithm, in one of the format accepted by interface \code{\link{nmf}}.
	}
	
	\item{na.rm}{
	  single logical that specifies if the rank for which the measures are NA values should 
	  be removed from the graph or not (default to \code{FALSE}).
	  This is useful when plotting results which include NAs due to error during the estimation 
	  process. See argument \code{stop} for \code{nmfEstimateRank} below.
	}
	
	\item{nrun}{
	  a \code{numeric} giving the number of run to perform for each value in \code{range}.
	}
	
	\item{range}{ a \code{numeric} vector containing the ranks of factorization to try.
	}
	
	\item{ref}{ reference object of class \code{NMF.rank}, as returned by 
	 function \code{nmfEstimateRank}.
	 The measures contained in \code{ref} are used and plotted as a reference.
	 The associated curves are drawn in \emph{red}, while those from \code{x} are 
	 drawn in \emph{blue}.
	 }
	
	\item{verbose}{ toggle verbosity.
	This parameter only affects the verbosity of the outer loop over the values in 
	\code{rank}. To print verbose (resp. debug) messages from each NMF run, one can use 
	\code{.options='v'} (resp. \code{.options='d'}) that will be passed to the 
	\code{\link{nmf}} method. 
	}
	
	\item{stop}{ logical flag for running the estimation process with fault tolerance.
	When \code{TRUE}, the whole execution will stop if any error is raised.
	When \code{FALSE} (default), the runs that raise an error will be skipped,  
	and the execution will carry on. The summary measures for the runs with errors are  
	set to NA values, and a warning is thrown.
	}
	
	\item{what}{ a \code{character} string that partially matches one of the 
	following item: \code{'all'}, \code{'cophenetic'}, \code{'rss'}, \code{'residuals'}
	, \code{'dispersion'}. 
	It specifies which measure must be plotted (\code{what='all'} plots all 
	the measures).
	}
	
	\item{x}{ 
	For \code{nmfEstimateRank} a target object to be estimated, in one of the format accepted
	  by interface \code{\link{nmf}}.
	  
	 For \code{plot.NMF.rank} an object of class \code{NMF.rank} as returned by 
	 function \code{nmfEstimateRank}.
	 }
	
	\item{\dots}{
	  For \code{nmfEstimateRank}, these are extra parameters passed to interface
	  \code{nmf}. Note that the same parameters are used for each value of the rank.
	  See \code{\link{nmf}}.
	  
	  For \code{plot.NMF.rank}, these are extra graphical parameter passed to 
	  the standard function \code{plot}. See \code{\link{plot}}.
	}
}
\details{
Given a NMF algorithm and the target matrix, a common way of estimating \eqn{r} 
is to try different values, compute some quality measures of the results, 
and choose the best value according to this quality criteria.
See \emph{Brunet et al. (2004)} and \emph{Hutchins et al. (2008)}.

The function \code{nmfEstimateRank} allows to launch this estimation procedure.
It performs multiple NMF runs for a range of rank of factorization and, for each, 
returns a set of quality measures together with the associated consensus matrix.

}
\value{
A S3 object (i.e. a list) of class \code{NMF.rank} with the following slots:
\item{measures }{a \code{data.frame} containing the quality measures for each 
rank of factorizations in \code{range}. Each row corresponds to a measure, 
each column to a rank.
}
\item{consensus }{ a \code{list} of consensus matrices, indexed by the rank of 
factorization (as a character string).}

\item{fit }{ a \code{list} of the fits, indexed by the rank of 
factorization (as a character string).}

}
\references{

 	\emph{Metagenes and molecular pattern discovery using matrix factorization}
	Brunet, J.~P., Tamayo, P., Golub, T.~R., and Mesirov, J.~P. (2004)
	Proc Natl Acad Sci U S A
	101(12), 4164--4169.

}
\author{ Renaud Gaujoux  }

\seealso{ nmf }
\examples{

set.seed(123456)
n <- 50; r <- 3; m <- 20
V <- syntheticNMF(n, r, m, noise=TRUE)

# Use a seed that will be set before each first run
\dontrun{res.estimate <- nmfEstimateRank(V, seq(2,5), method='brunet', nrun=10, seed=123456)}

# plot all the measures
\dontrun{plot(res.estimate)}
# or only one: e.g. the cophenetic correlation coefficient
\dontrun{plot(res.estimate, 'cophenetic')}

}
