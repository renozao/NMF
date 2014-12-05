\name{scale.NMF}
\alias{scale.NMF}
\title{Rescaling NMF Models}
\usage{
  \method{scale}{NMF} (x, center = c("basis", "coef"),
    scale = 1)
}
\arguments{
  \item{x}{an NMF object}

  \item{center}{either a numeric normalising vector \eqn{\delta}{delta}, or either
\code{'basis'} or \code{'coef'}, which respectively correspond to using the
column sums of the basis matrix or the inverse of the row sums of the
coefficient matrix as a normalising vector.
If numeric, \code{center} should be a single value or a vector of length the
rank of the NMF model, i.e. the number of columns in the basis matrix.}

  \item{scale}{scaling coefficient applied to \eqn{D}, i.e. the value of \eqn{\alpha}{alpha},
or, if \code{center='coef'}, the value of \eqn{1/\alpha}{1/alpha} (see section \emph{Details}).}
}
\value{
  an NMF object
}
\description{
  Rescales an NMF model keeping the fitted target matrix identical.
}
\details{
  Standard NMF models are identifiable modulo a scaling factor, meaning that the
basis components and basis profiles can be rescaled without changing the fitted
values:

\deqn{X = W_1 H_1 = (W_1 D) (D^{-1} H_1) = W_2 H_2}{X = W H = (W D) (D^-1 H)}
with \eqn{D= \alpha diag(1/\delta_1, \ldots, 1\delta_r)}{D= alpha * diag(1/delta_1, ..., 1/delta_r)}

The default call \code{scale(object)} rescales the basis NMF object so that each
column of the basis matrix sums up to one.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# random 3-rank 10x5 NMF model
x <- rnmf(3, 10, 5)

# rescale based on basis
colSums(basis(x))
colSums(basis(scale(x)))

rx <- scale(x, 'basis', 10)
colSums(basis(rx))
rowSums(coef(rx))

# rescale based on coef
rowSums(coef(x))
rowSums(coef(scale(x, 'coef')))
rx <- scale(x, 'coef', 10)
rowSums(coef(rx))
colSums(basis(rx))

# fitted target matrix is identical but the factors have been rescaled
rx <- scale(x, 'basis')
all.equal(fitted(x), fitted(rx))
all.equal(basis(x), basis(rx))
}

