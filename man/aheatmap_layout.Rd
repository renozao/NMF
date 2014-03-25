\name{aheatmap_layout}
\alias{aheatmap_layout}
\title{Annotated Heatmap Layout}
\usage{
  aheatmap_layout(layout = "daml")
}
\arguments{
  \item{layout}{layout specification that indicates the
  relative position of the dendrograms, annotations, matrix
  and labels. Two layouts can be defined: one horizontal,
  which relates to components associated to rows, and one
  vertical, which relates to components associated with
  columns. Each layout is specified as a character strings,
  composed of characters that encode the order of each
  component:

  \describe{ \item{d}{ for dendrogram;} \item{a}{ for
  annotation tracks;} \item{m}{ for data matrix;} \item{l}{
  for labels.} } The specification must contain one
  instance of each character.

  The default horizontal and vertical layout is
  \code{"daml"}, and can also be specified as \code{"."}.
  Separate layouts can be passed as a character vector with
  2 element (e.g., \code{c("daml", "mald")}), or as a
  single string, with layouts separated by \code{"|"}
  (e.g., \code{". | almd"}). If only one layout
  specification is passed, then it is used for both
  horizontal and vertical layouts.

  \strong{Examples:} \itemize{ \item \code{layout = "dlma"}
  puts labels at the leaves of the dendrograms and
  annotation track below or at the right of the data matrix
  \item \code{layout = ". | amld"} use the default layout
  for rows, put column annotation track on top of the data
  matrix, followed by column labels and dendrogram. }}
}
\description{
  Shows a diagram of an annotated heatmap layout
  specification.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

aheatmap_layout()
aheatmap_layout('amld')
aheatmap_layout('amld | .')
}

