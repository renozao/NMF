\name{silhouette.NMF}
\alias{silhouette.NMF}
\title{Silhouette of NMF Clustering}
\usage{
  \method{silhouette}{NMF} (x, what = NULL, ...)
}
\arguments{
  \item{x}{an NMF object, as returned by
  \code{\link{nmf}}.}

  \item{what}{defines the type of clustering the computed
  silhouettes are meant to assess: \code{'samples'} for the
  clustering of samples (i.e. the columns of the target
  matrix), \code{'features'} for the clustering of features
  (i.e. the rows of the target matrix), \code{'consensus'}
  for the consensus clustering of samples, and
  \code{'cmap'} for the consensus clustering of samples but
  ordered has in the default hierarchical clustering used
  by \code{\link{consensusmap}} when plotting the heatmap
  of the consensus matrix (for multi-run NMF fits).}

  \item{...}{extra arguments not used.}
}
\description{
  Silhouette of NMF Clustering
}
\examples{
res <- nmfCheck(nrun = 5)

# sample clustering from best fit
plot(silhouette(res))

# from consensus
plot(silhouette(res, what = 'consensus'))

# feature clustering
plot(silhouette(res, what = 'features'))

# average silhouette are computed in summary measures
summary(res)
}
\seealso{
  \code{\link[NMF]{predict}}
}

