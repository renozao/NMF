\name{nmfReport}
\alias{nmfReport}
\title{Run NMF Methods and Generate a Report}
\usage{
  nmfReport(x, rank, method, colClass = NULL, ...,
    output = NULL, template = NULL)
}
\arguments{
  \item{x}{target matrix}

  \item{rank}{factorization rank}

  \item{method}{list of methods to apply}

  \item{colClass}{reference class to assess accuracy}

  \item{...}{extra paramters passed to \code{\link{nmf}}}

  \item{output}{output HTML file}

  \item{template}{template Rmd file}
}
\value{
  a list with the following elements:
\item{fits}{the fit(s) for each method and each value of the rank.}
\item{accuracy}{a data.frame that contains the summary assessment measures,
for each fit.}
}
\description{
  Generates an HTML report from running a set of method on a given
target matrix, for a set of factorization ranks.
}
\details{
  The report is based on an .Rmd document \code{'report.Rmd'} stored in
the package installation sub-directory \code{scripts/}, and is compiled
using \pkg{knitr}.

At the beginning of the document, a file named \code{'functions.R'} is
looked for in the current directory, and sourced if present.
This enables the definition of custom NMF methods (see \code{\link{setNMFMethod}})
or setting global options.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

\dontrun{

x <- rmatrix(20, 10)
gr <- gl(2, 5)
nmfReport(x, 2:4, method = list('br', 'lee'), colClass = gr, nrun = 5)

}
}

