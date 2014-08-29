\name{profplot}
\alias{profplot}
\alias{profplot.default}
\title{Plotting Expression Profiles}
\usage{
  profplot(x, ...)

  \method{profplot}{default} (x, y,
    scale = c("none", "max", "c1"), match.names = TRUE,
    legend = TRUE, confint = TRUE, Colv, labels,
    annotation, ..., add = FALSE)
}
\arguments{
  \item{x}{a matrix or an NMF object from which is extracted the mixture
coefficient matrix. It is extracted from the best fit if \code{x} is the
results from multiple NMF runs.}

  \item{y}{a matrix or an NMF object from which is extracted the mixture
coefficient matrix.
It is extracted from the best fit if \code{y} is the results from multiple NMF runs.}

  \item{scale}{specifies how the data should be scaled before plotting.
If \code{'none'} or \code{NA}, then no scaling is applied and the "raw" data is plotted.
If \code{TRUE} or \code{'max'} then each row of both matrices
are normalised with their respective maximum values.
If \code{'c1'}, then each column of both matrix is scaled into proportions (i.e. to sum up to one).
Default is \code{'none'}.}

  \item{match.names}{a logical that indicates if the profiles in \code{y}
should be subset and/or re-ordered to match the profile names in \code{x}
(i.e. the rownames). This is attempted only when both \code{x} and \code{y}
have names.}

  \item{legend}{a logical that specifies whether drawing the legend or not, or
coordinates specifications passed to argument \code{x} of
\code{\link{legend}}, that specifies the position of the legend.}

  \item{confint}{logical that indicates if confidence intervals for the
R-squared should be shown in legend.}

  \item{Colv}{specifies the way the columns of \code{x} are ordered before
plotting. It is used only when \code{y} is missing.  It can be: \itemize{
\item a single numeric value, specifying the index of a row of \code{x},
that is used to order the columns by \code{x[, order(x[abs(Colv),])]}.
Decreasing order is specified with a negative index.  \item an integer
vector directly specifying the order itself, in which case the columns are
ordered by \code{x[, Colv]} \item a factor used to order the columns by
\code{x[, order(Colv)]} and as argument \code{annotation} if this latter is
missing or not \code{NA}.  \item any other object with a suitable
\code{order} method. The columns are by \code{x[, order(Colv)]} }}

  \item{labels}{a character vector containing labels for each sample (i.e.
each column of \code{x}). These are used for labelling the x-axis.}

  \item{annotation}{a factor annotating each sample (i.e. each column of
\code{x}). If not missing, a coloured raw is plotted under the x-axis and
annotates each sample accordingly. If argument \code{Colv} is a factor, then
it is used to annotate the plot, unless \code{annotation=NA}.}

  \item{...}{graphical parameters passed to \code{\link{matplot}} or \code{\link{matpoints}}.}

  \item{add}{logical that indicates if the plot should be added as points to a previous plot}
}
\description{
  Plotting Expression Profiles

  When using NMF for clustering in particular, one looks
  for strong associations between the basis and a priori
  known groups of samples. Plotting the profiles may
  highlight such patterns.
}
\details{
  The function can also be used to compare the profiles
  from two NMF models or mixture coefficient matrices. In
  this case, it draws a scatter plot of the paired
  profiles.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# create a random target matrix
v <- rmatrix(50, 10)

# fit a single NMF model
res <- nmf(v, 3)
profplot(res)

# ordering according to first profile
profplot(res, Colv=1) # increasing
profplot(res, Colv=-1) # decreasing

# fit a multi-run NMF model
res2 <- nmf(v, 3, nrun=3)
profplot(res2)

# draw a profile correlation plot: this show how the basis components are
# returned in an unpredictable order
profplot(res, res2)

# looking at all the correlations allow to order the components in a "common" order
profcor(res, res2)
}
\seealso{
  \code{\link{profcor}}
}
\keyword{aplot}

