\name{syntheticNMF}
\alias{syntheticNMF}
\title{Simulating Datasets}
\usage{
  syntheticNMF(n, r, p, offset = NULL, noise = TRUE,
    factors = FALSE, seed = NULL)
}
\arguments{
  \item{n}{number of rows of the target matrix.}

  \item{r}{specification of the factorization rank.
It may be a single \code{numeric}, in which case argument \code{p} is required
and \code{r} groups of samples are generated from a draw from a multinomial
distribution with equal probabilities, that provides their sizes.

It may also be a numerical vector, which contains the number of samples in
each class (i.e integers). In this case argument \code{p} is discarded
and forced to be the sum of \code{r}.}

  \item{p}{number of columns of the synthetic target matrix.
Not used if parameter \code{r} is a vector (see description of argument \code{r}).}

  \item{offset}{specification of a common offset to be added to the synthetic target
matrix, before noisification.
Its may be a numeric vector of length \code{n}, or a single numeric value that
is used as the standard deviation of a centred normal distribution from which
the actual offset values are drawn.}

  \item{noise}{a logical that indicate if noise should be added to the
matrix.}

  \item{factors}{a logical that indicates if the NMF factors should be return
together with the matrix.}

  \item{seed}{a single numeric value used to seed the random number generator
before generating the matrix.
The state of the RNG is restored on exit.}
}
\value{
  a matrix, or a list if argument \code{factors=TRUE}.

When \code{factors=FALSE}, the result is a matrix object, with the following attributes set:
\describe{
\item{coefficients}{the true underlying coefficient matrix (i.e. \code{H});}
\item{basis}{the true underlying coefficient matrix (i.e. \code{H});}
\item{offset}{the offset if any;}
\item{pData}{a \code{list} with one element \code{'Group'} that contains a factor
that indicates the true groups of samples, i.e. the most contributing basis component for each sample;}
\item{fData}{a \code{list} with one element \code{'Group'} that contains a factor
that indicates the true groups of features, i.e. the basis component
to which each feature contributes the most.}
}

Moreover, the result object is an \code{\link{ExposeAttribute}} object, which means that
relevant attributes are accessible via \code{$}, e.g., \code{res$coefficients}.
In particular, methods \code{\link{coef}} and \code{\link{basis}} will work as expected
and return the true underlying coefficient and basis matrices respectively.
}
\description{
  The function \code{syntheticNMF} generates random target
  matrices that follow some defined NMF model, and may be
  used to test NMF algorithms. It is designed to designed
  to produce data with known or clear classes of samples.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# generate a synthetic dataset with known classes: 50 features, 18 samples (5+5+8)
n <- 50
counts <- c(5, 5, 8)

# no noise
V <- syntheticNMF(n, counts, noise=FALSE)
\dontrun{aheatmap(V)}

# with noise
V <- syntheticNMF(n, counts)
\dontrun{aheatmap(V)}
}

