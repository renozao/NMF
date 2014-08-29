\name{cluster_mat}
\alias{cluster_mat}
\title{Cluster Matrix Rows in Annotated Heatmaps}
\usage{
  cluster_mat(mat, param, distfun, hclustfun, reorderfun,
    na.rm = TRUE, subset = NULL, verbose = FALSE)
}
\arguments{
  \item{mat}{original input matrix that has already been appropriately subset in
the caller function (\code{aheatmap})}

  \item{param}{clustering specifications}

  \item{distfun}{Default distance method/function}

  \item{hclustfun}{Default clustering (linkage) method/function}

  \item{reorderfun}{Default reordering function}

  \item{na.rm}{Logical that specifies if NA values should be removed}

  \item{subset}{index (integer) vector specifying the subset indexes used to
subset mat. This is required to be able to return the original indexes.}
}
\description{
  Cluster Matrix Rows in Annotated Heatmaps
}
\keyword{internal}

