\docType{methods}
\name{heatmap-NMF}
\alias{basismap}
\alias{basismap-methods}
\alias{basismap,NMFfitX-method}
\alias{basismap,NMF-method}
\alias{coefmap}
\alias{coefmap-methods}
\alias{coefmap,NMFfitX-method}
\alias{coefmap,NMF-method}
\alias{consensusmap}
\alias{consensusmap,matrix-method}
\alias{consensusmap-methods}
\alias{consensusmap,NMFfitX-method}
\alias{consensusmap,NMF-method}
\alias{heatmap-NMF}
\title{Heatmaps of NMF Factors}
\usage{
  basismap(object, ...)

  \S4method{basismap}{NMF}(object, color = "YlOrRd:50",
    ..., scale = "r1", Colv = NA, subsetRow = FALSE,
    annRow = "basis", main = "Basis components",
    info = FALSE)

  coefmap(object, ...)

  \S4method{coefmap}{NMF}(object, color = "YlOrRd:50", ...,
    scale = "c1", Rowv = NA, Colv = TRUE, annRow = NA,
    annCol = "basis", main = "Mixture coefficients",
    info = FALSE)

  consensusmap(object, ...)

  \S4method{consensusmap}{NMFfitX}(object, ...,
    annCol = c("basis", "consensus"),
    main = "Consensus matrix", info = FALSE)

  \S4method{consensusmap}{matrix}(object,
    color = "-RdYlBu", ...,
    distfun = function(x) as.dist(1 - x),
    hclustfun = "average", Rowv = TRUE, Colv = "Rowv",
    main = if (is.null(nr) || nr > 1) "Consensus matrix" else "Connectiviy matrix",
    info = FALSE)

  \S4method{coefmap}{NMFfitX}(object, ..., Colv = TRUE,
    annCol = c("basis", "consensus"))
}
\arguments{
  \item{object}{an object from which is extracted NMF
  factors or a consensus matrix}

  \item{...}{extra arguments passed to
  \code{\link{aheatmap}}.}

  \item{subsetRow}{Argument that specifies how to filter
  the rows that will appear in the heatmap. When
  \code{FALSE} (default), all rows are used. Besides the
  values supported by argument \code{subsetRow} of
  \code{\link{aheatmap}}, other possible values are:

  \itemize{ \item \code{TRUE}: only the rows that are
  basis-specific are used. The default selection method is
  from \cite{KimH2007}. This is equivalent to
  \code{subsetRow='kim'}.

  \item a single \code{character} string or numeric value
  that specifies the method to use to select the
  basis-specific rows, that should appear in the heatmap
  (cf. argument \code{method} for function
  \code{\link{extractFeatures}}).

  Note \code{\link{extractFeatures}} is called with
  argument \code{nodups=TRUE}, so that features that are
  selected for multiple components only appear once. }}

  \item{info}{if \code{TRUE} then the name of the algorithm
  that fitted the NMF model is displayed at the bottom of
  the plot, if available. Other wise it is passed as is to
  \code{aheatmap}.}

  \item{color}{colour specification for the heatmap.
  Default to palette '-RdYlBu2:100', i.e. reversed palette
  'RdYlBu2' (a slight modification of RColorBrewer's
  palette 'RdYlBu') with 100 colors. Possible values are:
  \itemize{ \item a character/integer vector of length
  greater than 1 that is directly used and assumed to
  contain valid R color specifications. \item a single
  color/integer (between 0 and 8)/other numeric value that
  gives the dominant colors. Numeric values are converted
  into a pallete by \code{rev(sequential_hcl(2, h = x, l =
  c(50, 95)))}. Other values are concatenated with the grey
  colour '#F1F1F1'. \item one of RColorBrewer's palette
  name, or one of 'RdYlBu2', 'rainbow', 'heat', 'topo',
  'terrain', 'cm'. } When the coluor palette is specified
  with a single value, and is negative or preceded a minus
  ('-'), the reversed palette is used. The number of breaks
  can also be specified after a colon (':'). For example,
  the default colour palette is specified as
  '-RdYlBu2:100'.}

  \item{scale}{character indicating how the values should
  scaled in either the row direction or the column
  direction. Note that the scaling is performed after
  row/column clustering, so that it has no effect on the
  row/column ordering. Possible values are: \itemize{ \item
  \code{"row"}: center and standardize each row separately
  to row Z-scores \item \code{"column"}: center and
  standardize each column separately to column Z-scores
  \item \code{"r1"}: scale each row to sum up to one \item
  \code{"c1"}: scale each column to sum up to one \item
  \code{"none"}: no scaling }}

  \item{Colv}{clustering specification(s) for the columns.
  It accepts the same values as argument \code{Rowv}
  (modulo the expected length for vector specifications),
  and allow specifying the
  distance/clustering/ordering/display parameters to be
  used for the \emph{columns only}. \code{Colv} may also be
  set to \code{"Rowv"}, in which case the dendrogram or
  ordering specifications applied to the rows are also
  applied to the columns. Note that this is allowed only
  for square input matrices, and that the row ordering is
  in this case by default reversed (\code{revC=TRUE}) to
  obtain the diagonal in the standard way (from top-left to
  bottom-right). See argument \code{Rowv} for other
  possible values.}

  \item{annRow}{specifications of row annotation tracks
  displayed as coloured columns on the left of the
  heatmaps. The annotation tracks are drawn from left to
  right. The same conversion, renaming and colouring rules
  as for argument \code{annCol} apply.}

  \item{main}{Main title as a character string or a grob.}

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

  \item{annCol}{specifications of column annotation tracks
  displayed as coloured rows on top of the heatmaps. The
  annotation tracks are drawn from bottom to top. A single
  annotation track can be specified as a single vector;
  multiple tracks are specified as a list, a data frame, or
  an
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
  object, in which case the phenotypic data is used
  (\code{pData(eset)}). Character or integer vectors are
  converted and displayed as factors. Unnamed tracks are
  internally renamed into \code{Xi}, with i being
  incremented for each unamed track, across both column and
  row annotation tracks. For each track, if no
  corresponding colour is specified in argument
  \code{annColors}, a palette or a ramp is automatically
  computed and named after the track's name.}

  \item{distfun}{default distance measure used in
  clustering rows and columns. Possible values are:
  \itemize{ \item \code{"correlation"} and all the
  distances supported by \code{\link{dist}} (e.g.
  \code{"euclidean"}). \item an object of class \code{dist}
  such as returned by \code{\link{dist}} or
  \code{\link{as.dist}}. }}

  \item{hclustfun}{default clustering method used to
  cluster rows and columns. Possible values are: \itemize{
  \item a method name (a character string) supported by
  \code{\link{hclust}} (e.g. \code{'average'}). \item an
  object of class \code{hclust} such as returned by
  \code{\link{hclust}} \item a dendrogram }}
}
\description{
  The NMF package ships an advanced heatmap engine
  implemented by the function \code{\link{aheatmap}}. Some
  convenience heatmap functions are implemented, that
  redefine default values for some of the arguments of
  \code{\link{aheatmap}}, to tune the output specifically
  for NMF models.

  \code{basimap} draws an annotated heatmap of the basis
  matrix.

  \code{coefmap} draws an annotated heatmap of the
  coefficient matrix.

  \code{consensusmap} plots heatmaps of consensus matrices.
}
\details{
  \code{basimap} default values for the following arguments
  of \code{\link{aheatmap}}: \itemize{ \item the color
  palette; \item the scaling specification, which by
  default scales each row separately so that they sum up to
  one (\code{scale='r1'}); \item the column ordering which
  is disabled; \item allowing for passing feature
  extraction methods in argument \code{subsetRow}, that are
  passed to \code{\link{extractFeatures}}. See argument
  description here and therein. \item the addition of a
  default named annotation track, that shows the dominant
  basis component for each row (i.e. each feature).

  This track, named \code{'basis'}, is specified as a
  single character string, that is looked-up in argument
  \code{annRow}, which may be of any type supported by
  \code{\link{aheatmap}}. The only difference is that
  character annotation vectors that have elements other
  than \code{'basis'}, are not are not allowed as-is, and
  must be encapsulated in a list (see examples). \item a
  suitable title and extra information like the fitting
  algorithm, when \code{object} is a fitted NMF model. }

  \code{coefmap} redefines default values for the following
  arguments of \code{\link{aheatmap}}: \itemize{ \item the
  color palette; \item the scaling specification, which by
  default scales each column separately so that they sum up
  to one (\code{scale='c1'}); \item the row ordering which
  is disabled; \item the addition of a default annotation
  track, that shows the most contributing basis component
  for each column (i.e. each sample).

  This track is specified in argument \code{annCol}, which
  behaves like argument \code{annRow} in \code{basismap}
  (see above). A matching row annotation track may be
  displayed using \code{annRow='basis'}. \item a suitable
  title and extra information like the fitting algorithm,
  when \code{object} is a fitted NMF model. }

  \code{consensusmap} redefines default values for the
  following arguments of \code{\link{aheatmap}}: \itemize{
  \item the color palette; \item the column ordering which
  is set equal to the row ordering, since a consensus
  matrix is symetric; \item the distance and linkage
  methods used to order the rows (and columns). The default
  is to use 1 minus the consensus matrix itself as
  distance, and average linkage. \item the addition of two
  default named annotation tracks, \code{'basis'} and
  \code{'consensus'}, that show, for each column (i.e. each
  sample), the dominant basis component in the best fit and
  the hierarchical clustering of the consensus matrix
  respectively (using 1-consensus as distance and average
  linkage).

  This track is specified in argument \code{annCol}, which
  behaves like argument \code{annRow} in
  \code{\link{basismap}}. \item a suitable title and extra
  information like the type of NMF model or the fitting
  algorithm, when \code{object} is a fitted NMF model. }
}
\section{Methods}{
  \describe{

  \item{basismap}{\code{signature(object = "NMF")}: Plots a
  heatmap of the basis matrix of the NMF model
  \code{object}. This method also works for fitted NMF
  models (i.e. \code{NMFfit} objects). }

  \item{basismap}{\code{signature(object = "NMFfitX")}:
  Plots a heatmap of the basis matrix of the best fit in
  \code{object}. }

  \item{coefmap}{\code{signature(object = "NMFfitX")}:
  Plots a heatmap of the coefficient matrix of the best fit
  in \code{object}. }

  \item{consensusmap}{\code{signature(object = "NMFfitX")}:
  Plots a heatmap of the consensus matrix obtained when
  fitting an NMF model with multiple runs. }

  \item{consensusmap}{\code{signature(object = "NMF")}:
  Plots a heatmap of the connectivity matrix of an NMF
  model. }

  \item{consensusmap}{\code{signature(object = "matrix")}:
  Main method that redefines default values for arguments
  of \code{\link{aheatmap}}. }

  }
}
\examples{
# random data with underlying NMF model
v <- syntheticNMF(50, 3, 20)
# estimate a model
x <- nmf(v, 3)

# show basis matrix
basismap(x)
basismap(x, annRow=NA) # no annotation track

# character annotation vector: ok if it does not contain 'basis'
basismap(x, annRow=c('alpha', 'beta')) # annotate first and second row
try( basismap(x, annRow=c('alpha', 'beta', 'basis')) ) # error
basismap(x, annRow=list(c('alpha', 'beta', 'basis'))) # ok
basismap(x, annRow=list('basis', c('alpha', 'beta', 'basis'))) # ok

# changing the name of the basis annotation
basismap(x, annRow=c(new_name='basis'))
# coefficient matrix
coefmap(x)
coefmap(x, annCol=NA)
coefmap(x, annCol=c('alpha', 'beta')) # annotate first and second sample
coefmap(x, annCol=list('basis', Greek=c('alpha', 'beta'))) # annotate first and second sample + basis annotation
coefmap(x, annCol=c(new_name='basis'))
}
\keyword{methods}

