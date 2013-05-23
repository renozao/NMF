\name{aheatmap}
\alias{aheatmap}
\title{Annotated Heatmaps}
\description{
  The function \code{aheatmap} plots high-quality heatmaps,
  with a detailed legend and unlimited annotation tracks
  for both columns and rows. The annotations are coloured
  differently according to their type (factor or numeric
  covariate). Although it uses grid graphics, the generated
  plot is compatible with base layouts such as the ones
  defined with \code{'mfrow'} or \code{\link{layout}},
  enabling the easy drawing of multiple heatmaps on a
  single a plot -- at last!.
}
\usage{
  aheatmap(x, color = "-RdYlBu2:100", breaks = NA,
  border_color = NA, cellwidth = NA, cellheight = NA,
  scale = "none", Rowv = TRUE, Colv = TRUE,
  revC = identical(Colv, "Rowv") || is_NA(Rowv) || (is.integer(Rowv) && length(Rowv) > 1),
  distfun = "euclidean", hclustfun = "complete",
  reorderfun = function(d, w) reorder(d, w),
  treeheight = 50, legend = TRUE, annCol = NA, annRow = NA,
  annColors = NA, annLegend = TRUE, labRow = NULL,
  labCol = NULL, subsetRow = NULL, subsetCol = NULL,
  fontsize = 10, cexRow = min(0.2 + 1/log10(nr), 1.2),
  cexCol = min(0.2 + 1/log10(nc), 1.2), filename = NA,
  width = NA, height = NA, main = NULL, sub = NULL,
  info = NULL, verbose = getOption("verbose"), ...)
}
\arguments{
  \item{x}{numeric matrix of the values to be plotted. An
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
  objects can also be passed, in which case the expression
  values are plotted (\code{exprs(x)}).}

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

  \item{breaks}{a sequence of numbers that covers the range
  of values in \code{x} and is one element longer than
  color vector. Used for mapping values to colors. Useful,
  if needed to map certain values to certain colors. If
  value is NA then the breaks are calculated automatically.
  If \code{breaks} is a single value, then the colour
  palette is centered on this value.}

  \item{border_color}{color of cell borders on heatmap, use
  NA if no border should be drawn.}

  \item{cellwidth}{individual cell width in points. If left
  as NA, then the values depend on the size of plotting
  window.}

  \item{cellheight}{individual cell height in points. If
  left as NA, then the values depend on the size of
  plotting window.}

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

  \item{revC}{a logical that specify if the \emph{row
  order} defined by \code{Rowv} should be reversed. This is
  mainly used to get the rows displayed from top to bottom,
  which is not the case by default. Its default value is
  computed at runtime, to suit common situations where
  natural ordering is a more sensible choice: no or fix
  ordering of the rows (\code{Rowv=NA} or an integer vector
  of indexes -- of length > 1), and when a symmetric
  ordering is requested -- so that the diagonal is shown as
  expected. An argument in favor of the "odd" default
  display (bottom to top) is that the row dendrogram is
  plotted from bottom to top, and reversing its reorder may
  take a not too long but non negligeable time.}

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

  \item{reorderfun}{default dendrogram reordering function,
  used to reorder the dendrogram, when either \code{Rowv}
  or \code{Colv} is a numeric weight vector, or provides or
  computes a dendrogram. It must take 2 parameters: a
  dendrogram, and a weight vector.}

  \item{subsetRow}{Specification of subsetting the rows
  before drawing the heatmap. Possible values are:
  \itemize{ \item an integer vector of length > 1
  specifying the indexes of the rows to keep; \item a
  character vector of length > 1 specyfing the names of the
  rows to keep. These are the original rownames, not the
  names specified in \code{labRow}. \item a logical vector
  of length > 1, whose elements are recycled if the vector
  has not as many elements as rows in \code{x}. } Note that
  in the case \code{Rowv} is a dendrogram or hclust object,
  it is first converted into an ordering vector, and cannot
  be displayed -- and a warning is thrown.}

  \item{subsetCol}{Specification of subsetting the columns
  before drawing the heatmap. It accepts the similar values
  as \code{subsetRow}. See details above.}

  \item{treeheight}{how much space (in points) should be
  used to display dendrograms. If specified as a single
  value, it is used for both dendrograms. A length-2 vector
  specifies separate values for the row and column
  dendrogram respectively. Default value: 50 points.}

  \item{legend}{boolean value that determines if a colour
  ramp for the heatmap's colour palette should be drawn or
  not. Default is \code{TRUE}.}

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

  \item{annRow}{specifications of row annotation tracks
  displayed as coloured columns on the left of the
  heatmaps. The annotation tracks are drawn from left to
  right. The same conversion, renaming and colouring rules
  as for argument \code{annCol} apply.}

  \item{annColors}{list for specifying annotation track
  colors manually. It is possible to define the colors for
  only some of the annotations. Check examples for
  details.}

  \item{annLegend}{boolean value specifying if the legend
  for the annotation tracks should be drawn or not. Default
  is \code{TRUE}.}

  \item{labRow}{labels for the rows.}

  \item{labCol}{labels for the columns. See description for
  argument \code{labRow} for a list of the possible
  values.}

  \item{fontsize}{base fontsize for the plot}

  \item{cexRow}{fontsize for the rownames, specified as a
  fraction of argument \code{fontsize}.}

  \item{cexCol}{fontsize for the colnames, specified as a
  fraction of argument \code{fontsize}.}

  \item{main}{Main title as a character string or a grob.}

  \item{sub}{Subtitle as a character string or a grob.}

  \item{info}{(experimental) Extra information as a
  character vector or a grob.  If \code{info=TRUE},
  information about the clustering methods is displayed at
  the bottom of the plot.}

  \item{filename}{file path ending where to save the
  picture. Currently following formats are supported: png,
  pdf, tiff, bmp, jpeg. Even if the plot does not fit into
  the plotting window, the file size is calculated so that
  the plot would fit there, unless specified otherwise.}

  \item{width}{manual option for determining the output
  file width in}

  \item{height}{manual option for determining the output
  file height in inches.}

  \item{verbose}{if \code{TRUE} then verbose messages are
  displayed and the borders of some viewports are
  highlighted. It is entended for debugging purposes.}

  \item{\dots}{graphical parameters for the text used in
  plot. Parameters passed to \code{\link{grid.text}}, see
  \code{\link{gpar}}.}
}
\details{
  The development of this function started as a fork of the
  function \code{\link{aheatmap}} from the package
  \code{pheatmap}, and provides several enhancements such
  as: \itemize{ \item it replicates the arguments of the
  base \code{\link{heatmap}} function \item annotation of
  both columns and rows, with simplified interface \item
  easy specification of clustering methods and colors \item
  return clustering data and grid grob object }

  Please read the associated vignette for more information
  and sample code.
}
\examples{
# Generate random data
n <- 50; p <- 20
x <- abs(rmatrix(n, p, rnorm, mean=4, sd=1))
x[1:10, seq(1, 10, 2)] <- x[1:10, seq(1, 10, 2)] + 3
x[11:20, seq(2, 10, 2)] <- x[11:20, seq(2, 10, 2)] + 2
rownames(x) <- paste("ROW", 1:n)
colnames(x) <- paste("COL", 1:p)

## Default heatmap
aheatmap(x)

## Distance methods
aheatmap(x, Rowv = "correlation")
aheatmap(x, Rowv = "man") # partially matched to 'manhattan'
aheatmap(x, Rowv = "man", Colv="binary")

## Scaling
aheatmap(x, scale = "row")
aheatmap(x, scale = "col") # partially matched to 'column'
aheatmap(x, scale = "r1") # each row sum up to 1
aheatmap(x, scale = "c1") # each colum sum up to 1

## Heatmap colors
aheatmap(x, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
# color specification as an integer: use R basic colors
aheatmap(x, color = 1L)
# color specification as a negative integer: use reverse basic palette
aheatmap(x, color = -1L)
# color specification as a numeric: use HCL color
aheatmap(x, color = 1)
# do not cluster the rows
aheatmap(x, Rowv = NA)
# no heatmap legend
aheatmap(x, legend = FALSE)
# cell and font size
aheatmap(x, cellwidth = 10, cellheight = 5)
aheatmap(x, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "aheatmap.pdf")

# Generate column annotations
annotation = data.frame(Var1 = factor(1:10 \\\%\\\% 2 == 0, labels = c("Class1", "Class2")), Var2 = 1:10)

aheatmap(x, annCol = annotation)
aheatmap(x, annCol = annotation, annLegend = FALSE)


# Specify colors
Var1 = c("navy", "darkgreen")
names(Var1) = c("Class1", "Class2")
Var2 = c("lightgreen", "navy")

ann_colors = list(Var1 = Var1, Var2 = Var2)

aheatmap(x, annCol = annotation, annColors = ann_colors)

# Specifying clustering from distance matrix
drows = dist(x, method = "minkowski")
dcols = dist(t(x), method = "minkowski")
aheatmap(x, Rowv = drows, Colv = dcols)
}
\author{
  Original version of \code{\link[pheatmap]{pheatmap}}:
  Raivo Kolde <rkolde@gmail.com>

  Enhancement into \code{aheatmap}: Renaud Gaujoux
}

