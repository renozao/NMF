\name{nmfEstimateRank}
\alias{nmfEstimateRank}
\alias{plot.NMF.rank}
\title{Estimate Rank for NMF Models}
\usage{
  nmfEstimateRank(x, range,
    method = nmf.getOption("default.algorithm"), nrun = 30,
    model = NULL, ..., verbose = FALSE, stop = FALSE)

  \method{plot}{NMF.rank} (x, y = NULL,
    what = c("all", "cophenetic", "rss", "residuals", "dispersion", "evar", "sparseness", "sparseness.basis", "sparseness.coef"),
    na.rm = FALSE, ...)
}
\arguments{
  \item{x}{For \code{nmfEstimateRank} a target object to be
  estimated, in one of the format accepted by interface
  \code{\link{nmf}}.

  For \code{plot.NMF.rank} an object of class
  \code{NMF.rank} as returned by function
  \code{nmfEstimateRank}.}

  \item{range}{a \code{numeric} vector containing the ranks
  of factorization to try. Note that duplicates are removed
  and values are sorted in increasing order. The results
  are notably returned in this order.}

  \item{method}{A single NMF algorithm, in one of the
  format accepted by the function \code{\link{nmf}}.}

  \item{nrun}{a \code{numeric} giving the number of run to
  perform for each value in \code{range}.}

  \item{verbose}{toggle verbosity.  This parameter only
  affects the verbosity of the outer loop over the values
  in \code{range}. To print verbose (resp. debug) messages
  from each NMF run, one can use \code{.options='v'} (resp.
  \code{.options='d'}) that will be passed to the function
  \code{\link{nmf}}.}

  \item{stop}{logical flag for running the estimation
  process with fault tolerance.  When \code{TRUE}, the
  whole execution will stop if any error is raised.  When
  \code{FALSE} (default), the runs that raise an error will
  be skipped, and the execution will carry on. The summary
  measures for the runs with errors are set to NA values,
  and a warning is thrown.}

  \item{...}{For \code{nmfEstimateRank}, these are extra
  parameters passed to interface \code{nmf}. Note that the
  same parameters are used for each value of the rank.  See
  \code{\link{nmf}}.

  For \code{plot.NMF.rank}, these are extra graphical
  parameter passed to the standard function \code{plot}.
  See \code{\link{plot}}.}

  \item{na.rm}{single logical that specifies if the rank
  for which the measures are NA values should be removed
  from the graph or not (default to \code{FALSE}).  This is
  useful when plotting results which include NAs due to
  error during the estimation process. See argument
  \code{stop} for \code{nmfEstimateRank}.}

  \item{y}{reference object of class \code{NMF.rank}, as
  returned by function \code{nmfEstimateRank}. The measures
  contained in \code{ref} are used and plotted as a
  reference. It is typically used to plot results obtained
  on randomized data. The associated curves are drawn in
  \emph{red}, while those from \code{x} are drawn in
  \emph{blue}.}

  \item{what}{a \code{character} string that partially
  matches one of the following item: \sQuote{all},
  \sQuote{cophenetic}, \sQuote{rss}, \sQuote{residuals} ,
  \sQuote{dispersion}. It specifies which measure must be
  plotted (\code{what='all'} plots all the measures).}
}
\value{
  \code{nmfEstimateRank} returns a S3 object (i.e. a list)
  of class \code{NMF.rank} with the following elements:

  \item{measures }{a \code{data.frame} containing the
  quality measures for each rank of factorizations in
  \code{range}. Each row corresponds to a measure, each
  column to a rank. } \item{consensus }{ a \code{list} of
  consensus matrices, indexed by the rank of factorization
  (as a character string).} \item{fit }{ a \code{list} of
  the fits, indexed by the rank of factorization (as a
  character string).}
}
\description{
  A critical parameter in NMF algorithms is the
  factorization rank \eqn{r}. It defines the number of
  basis effects used to approximate the target matrix.
  Function \code{nmfEstimateRank} helps in choosing an
  optimal rank by implementing simple approaches proposed
  in the literature.
}
\details{
  Note that from version \emph{0.7}, one can equivalently
  call the function \code{\link{nmf}} with a range of
  ranks.

  Given a NMF algorithm and the target matrix, a common way
  of estimating \eqn{r} is to try different values, compute
  some quality measures of the results, and choose the best
  value according to this quality criteria. See
  \cite{Brunet et al. (2004)} and \cite{Hutchins et al.
  (2008)}.

  The function \code{nmfEstimateRank} allows to perform
  this estimation procedure. It performs multiple NMF runs
  for a range of rank of factorization and, for each,
  returns a set of quality measures together with the
  associated consensus matrix.

  In order to avoid overfitting, it is recommended to run
  the same procedure on randomized data. The results on the
  original and the randomised data may be plotted on the
  same plots, using argument \code{y}.
}
\examples{
set.seed(123456)
n <- 50; r <- 3; m <- 20
V <- syntheticNMF(n, r, m)

# Use a seed that will be set before each first run
\dontrun{
res <- nmfEstimateRank(V, seq(2,5), method='brunet', nrun=10, seed=123456)
# or equivalently
res <- nmf(V, seq(2,5), method='brunet', nrun=10, seed=123456)
}

# plot all the measures
\dontrun{plot(res)}
# or only one: e.g. the cophenetic correlation coefficient
\dontrun{plot(res, 'cophenetic')}

# run same estimation on randomized data
\dontrun{
rV <- randomize(V)
rand <- nmfEstimateRank(rV, seq(2,5), method='brunet', nrun=10, seed=123456)
plot(res, rand)
}
}
\references{
  Brunet J, Tamayo P, Golub TR and Mesirov JP (2004).
  "Metagenes and molecular pattern discovery using matrix
  factorization." _Proceedings of the National Academy of
  Sciences of the United States of America_, *101*(12), pp.
  4164-9. ISSN 0027-8424, <URL:
  http://dx.doi.org/10.1073/pnas.0308531101>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/15016911>.

  Hutchins LN, Murphy SM, Singh P and Graber JH (2008).
  "Position-dependent motif characterization using
  non-negative matrix factorization." _Bioinformatics
  (Oxford, England)_, *24*(23), pp. 2684-90. ISSN
  1367-4811, <URL:
  http://dx.doi.org/10.1093/bioinformatics/btn526>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/18852176>.
}

