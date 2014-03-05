\docType{data}
\name{esGolub}
\alias{esGolub}
\title{Golub ExpressionSet}
\format{There are 3 covariates listed.

\itemize{

\item Samples: The original sample labels.  \item ALL.AML: Whether the
patient had AML or ALL. It is a \code{\link{factor}} with levels
\code{c('ALL', 'AML')}.  \item Cell: ALL arises from two different types of
lymphocytes (T-cell and B-cell).  This specifies which for the ALL patients;
There is no such information for the AML samples. It is a
\code{\link{factor}} with levels \code{c('T-cell', 'B-cell', NA)}.

}}
\source{
  Web page for \cite{Brunet2004}:\cr
  \url{http://www.broadinstitute.org/publications/broad872}

  Original data from Golub et al.:\cr
  \url{http://www-genome.wi.mit.edu/mpr/data_set_ALL_AML.html}
}
\description{
  This data comes originally from the gene expression data
  from \cite{Golub et al. (1999)}. The version included in
  the package is the one used and referenced in
  \cite{Brunet et al. (2004)}. The samples are from 27
  patients with acute lymphoblastic leukemia (ALL) and 11
  patients with acute myeloid leukemia (AML).
}
\details{
  The samples were assayed using Affymetrix Hgu6800 chips
  and the original data on the expression of 7129 genes
  (Affymetrix probes) are available on the Broad Institute
  web site (see references below).

  The data in \code{esGolub} were obtained from the web
  page related to the paper from \cite{Brunet et al.
  (2004)}, which describes an application of Nonnegative
  Matrix Factorization to gene expression clustering. (see
  link in section \emph{Source}).

  They contain the 5,000 most highly varying genes
  according to their coefficient of variation, and were
  installed in an object of class
  \code{\link[Biobase]{ExpressionSet-class}}.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

data(esGolub)
esGolub
\dontrun{pData(esGolub)}
}
\references{
  Golub TR, Slonim DK, Tamayo P, Huard C, Gaasenbeek M,
  Mesirov JP, Coller H, Loh ML, Downing JR, Caligiuri Ma,
  Bloomfield CD and Lander ES (1999). "Molecular
  classification of cancer: class discovery and class
  prediction by gene expression monitoring." _Science (New
  York, N.Y.)_, *286*(5439), pp. 531-7. ISSN 0036-8075,
  <URL: http://www.ncbi.nlm.nih.gov/pubmed/10521349>.

  Brunet J, Tamayo P, Golub TR and Mesirov JP (2004).
  "Metagenes and molecular pattern discovery using matrix
  factorization." _Proceedings of the National Academy of
  Sciences of the United States of America_, *101*(12), pp.
  4164-9. ISSN 0027-8424, <URL:
  http://dx.doi.org/10.1073/pnas.0308531101>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/15016911>.
}
\keyword{datasets}

