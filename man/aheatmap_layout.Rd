\name{aheatmap_layout}
\alias{aheatmap_layout}
\title{Annotated Heatmap Layout Preview}
\usage{
  aheatmap_layout(layout = "daml", size = NULL)
}
\arguments{
  \item{layout}{layout specification that indicates the relative position
of the heatmap's components.
Two layouts can be defined: one horizontal, which relates to components associated to rows,
and one vertical, which relates to components associated with columns.
Each layout is specified as a character strings, composed of characters
that encode the order of each component: dendrogram (d), annotation tracks (a),
data matrix (m), labels (l) and legend (L).
See section \emph{Layout syntax} for a complete specification}

  \item{size}{list defining the size of each component (mainly for internal use).}
}
\description{
  Shows a diagram of an annotated heatmap layout for given
  specification.
}
\details{
  Layout syntax:

  Layouts are specified as character strings that can
  contain the following characters, each associated with a
  given component or behaviour:

  \strong{Components} \describe{ \item{\sQuote{d}}{
  dendrogram component} \item{\sQuote{a}}{ annotation
  tracks} \item{\sQuote{m}}{ data matrix}
  \item{\sQuote{l}}{ labels} \item{\sQuote{L}}{ scale
  legend} }

  \strong{Behaviours} \describe{ \item{\sQuote{^}}{ align
  top (resp. left) for horizontal (resp. vertical) layout.}
  \item{\sQuote{-}}{ align middle (resp. center) for
  horizontal (resp. vertical) layout.} \item{\sQuote{_}}{
  align bottom (resp. right) for horizontal (resp.
  vertical) layout. If used alone (i.e. \code{layout =
  "_"}), then this is equivalent to \code{"|.L_"}, which
  places the legend horizontally on the bottom-right
  corner.} \item{\sQuote{*}}{ used either alone or after
  after \sQuote{L} to specifiy that the legend should
  expand to full height/width.} } The specification must
  contain one instance of each of these character.

  The default horizontal/vertical layout is \code{"daml"},
  and can also be specified as \code{"."}.

  Separate layouts can be passed as a character vector with
  2 element (e.g., \code{c("daml", "mald")}), or as a
  single string, with layouts separated by \code{"|"}
  (e.g., \code{"almd | L."}). When using this separator, a
  layout specification may be omitted, indicating that the
  default layout shoud be used: \code{"almd|"} is
  equivalent to \code{"almd | ."}. If only one layout
  specification is passed (i.e. a string without
  \code{"|"}), then it is used for both horizontal and
  vertical layouts.

  \strong{Shortcuts} \itemize{ \item \code{layout = "*"} is
  a shortcut for \code{layout = ".L*"}, which expands the
  legend to take up full height; \item \code{layout = "_"}
  is a shortcut for \code{layout = "|.L_"}, which puts the
  legend at bottom-right corner; \item \code{layout = "_^"}
  is a shortcut for \code{layout = "|.L^"}, which puts the
  legend at bottom-left corner; \item \code{layout = "_*"}
  is a shortcut for \code{layout = "|.L*"}, which puts the
  legend at bottom, expanded to take up full width; \item
  \code{layout = "^"} is a shortcut for \code{layout =
  "L^.|"}, which puts the legend on the top-left corner. }

  \strong{Examples:} \itemize{ \item \code{layout = "dlma"}
  puts labels at the leaves of the dendrograms and
  annotation track below or at the right of the data matrix
  \item \code{layout = ". | amld"} use the default layout
  for rows, put column annotation track on top of the data
  matrix, followed by column labels and dendrogram. }
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# default layout
aheatmap_layout()

# Common row/column layout: annotations > data > labels > dendrogram
aheatmap_layout('amld')

# Separate row/column layout: row as above / column as default
aheatmap_layout('amld | .')

## Legend
# horizontal bottom-right
aheatmap_layout('_')
# hotizontal top-left (equivalent to "|L^.")
aheatmap_layout('^')
}

