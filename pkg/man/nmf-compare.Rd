\docType{methods}
\name{compare-NMF}
\alias{compare,list-method}
\alias{compare-NMF}
\alias{compare,NMFfit-method}
\alias{consensusmap,list-method}
\alias{consensusmap,NMF.rank-method}
\alias{plot,NMFList,missing-method}
\alias{summary,NMFList-method}
\title{Comparing Results from Different NMF Runs}
\usage{
  \S4method{compare}{NMFfit}(object, ...)

  \S4method{compare}{list}(object, ...)

  \S4method{summary}{NMFList}(object, sort.by = NULL,
    select = NULL, ...)

  \S4method{plot}{NMFList,missing}(x, y, skip = -1, ...)

  \S4method{consensusmap}{NMF.rank}(object, ...)

  \S4method{consensusmap}{list}(object, layout, ...,
    Rowv = FALSE, main = names(object))
}
\arguments{
  \item{...}{extra arguments passed by \code{compare} to
  \code{summary,NMFList} or to the \code{summary} method of
  each fit.}

  \item{select}{the columns to be output in the result
  \code{data.frame}.  The column are given by their names
  (partially matched).  The column names are the names of
  the summary measures returned by the \code{summary}
  methods of the corresponding NMF results.}

  \item{sort.by}{the sorting criteria, i.e. a partial match
  of a column name, by which the result \code{data.frame}
  is sorted.  The sorting direction (increasing or
  decreasing) is computed internally depending on the
  chosen criteria (e.g. decreasing for the cophenetic
  coefficient, increasing for the residuals).}

  \item{x}{an \code{NMFList} object that contains fits from
  separate NMF runs.}

  \item{y}{missing}

  \item{layout}{specification of the layout. It may be a
  single numeric or a numeric couple, to indicate a square
  or rectangular layout respectively, that is filled row by
  row. It may also be a matrix that is directly passed to
  the function \code{\link[graphics]{layout}} from the
  package \code{graphics}.}

  \item{object}{an object computed using some algorithm, or
  that describes an algorithm itself.}

  \item{skip}{an integer that indicates the number of
  points to skip/remove from the beginning of the curve. If
  \code{skip=1L} (default) only the initial residual --
  that is computed before any iteration, is skipped, if
  present in the track (it associated with iteration 0).}

  \item{Rowv}{clustering specification(s) for the rows. It
  allows to specify the
  distance/clustering/ordering/display parameters to be
  used for the \emph{rows only}. Possible values are:
  \itemize{ \item \code{TRUE} or \code{NULL} (to be
  consistent with \code{\link{heatmap}}): compute a
  dendrogram from hierarchical clustering using the
  distance and clustering methods \code{distfun} and
  \code{hclustfun}.

  \item \code{NA}: disable any ordering. In this case, and
  if not otherwise specified with argument
  \code{revC=FALSE}, the heatmap shows the input matrix
  with the rows in their original order, with the first row
  on top to the last row at the bottom. Note that this
  differ from the behaviour or \code{\link{heatmap}}, but
  seemed to be a more sensible choice when vizualizing a
  matrix without reordering.

  \item an integer vector of length the number of rows of
  the input matrix (\code{nrow(x)}), that specifies the row
  order. As in the case \code{Rowv=NA}, the ordered matrix
  is shown first row on top, last row at the bottom.

  \item a character vector or a list specifying values to
  use instead of arguments \code{distfun}, \code{hclustfun}
  and \code{reorderfun} when clustering the rows (see the
  respective argument descriptions for a list of accepted
  values). If \code{Rowv} has no names, then the first
  element is used for \code{distfun}, the second (if
  present) is used for \code{hclustfun}, and the third (if
  present) is used for \code{reorderfun}.

  \item a numeric vector of weights, of length the number
  of rows of the input matrix, used to reorder the
  internally computed dendrogram \code{d} by
  \code{reorderfun(d, Rowv)}.

  \item \code{FALSE}: the dendrogram \emph{is} computed
  using methods \code{distfun}, \code{hclustfun}, and
  \code{reorderfun} but is not shown.

  \item a single integer that specifies how many subtrees
  (i.e. clusters) from the computed dendrogram should have
  their root faded out. This can be used to better
  highlight the different clusters.

  \item a single double that specifies how much space is
  used by the computed dendrogram. That is that this value
  is used in place of \code{treeheight}. }}

  \item{main}{Main title as a character string or a grob.}
}
\description{
  The functions documented here allow to compare the fits
  computed in different NMF runs. The fits do not need to
  be from the same algorithm, nor have the same dimension.
}
\details{
  The methods \code{compare} enables to compare multiple
  NMF fits either passed as arguments or as a list of fits.
  These methods eventually call the method
  \code{summary,NMFList}, so that all its arguments can be
  passed \strong{named} in \code{...}.
}
\section{Methods}{
  \describe{

  \item{compare}{\code{signature(object = "NMFfit")}:
  Compare multiple NMF fits passed as arguments. }

  \item{compare}{\code{signature(object = "list")}:
  Compares multiple NMF fits passed as a standard list. }

  \item{consensusmap}{\code{signature(object =
  "NMF.rank")}: Draw a single plot with a heatmap of the
  consensus matrix obtained for each value of the rank, in
  the range tested with \code{\link{nmfEstimateRank}}. }

  \item{consensusmap}{\code{signature(object = "list")}:
  Draw a single plot with a heatmap of the consensus matrix
  of each element in the list \code{object}. }

  \item{plot}{\code{signature(x = "NMFList", y =
  "missing")}: \code{plot} plot on a single graph the
  residuals tracks for each fit in \code{x}. See function
  \code{\link{nmf}} for details on how to enable the
  tracking of residuals. }

  \item{summary}{\code{signature(object = "NMFList")}:
  \code{summary,NMFList} computes summary measures for each
  NMF result in the list and return them in rows in a
  \code{data.frame}. By default all the measures are
  included in the result, and \code{NA} values are used
  where no data is available or the measure does not apply
  to the result object (e.g. the dispersion for single' NMF
  runs is not meaningful). This method is very useful to
  compare and evaluate the performance of different
  algorithms. }

  }
}
\examples{
x <- rmatrix(20,10)
res <- nmf(x, 3)
res2 <- nmf(x, 2, 'lee')

# compare arguments
compare(res, res2, target=x)
# compare elements of a list
compare(list(res, res2), target=x)
}
\keyword{methods}

