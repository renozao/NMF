\name{runNmf}
\alias{runNmf}
\title{A function to run NMF and extract important features from the results.}
\description{This function is designed to be used once the user has defined the rank factor and algorithm most suitable for his/her data. Besides running the \code{\link{nmf}} function, it also uses the \code{\link{exploreNMF}} to extract metagene features and further information.}
\usage{runNmf(original.data, ann, r, prefix, .opt = "vP", alg, maxIter = 30000, nrun = 1000)}
\arguments{
  \item{original.data}{The original data matrix submiited to the \code{nmf} analysis.}
  \item{ann}{Information about the \code{data} to be displayed on the heatmaps.}
  \item{r}{Factorization rank to be used.}
  \item{prefix}{A string used to name the generated plots and the folder to save them.}
  \item{.opt}{this argument is used to set runtime options.

It can be a \code{list} containing named options with their values, or, in
the case only boolean/integer options need to be set, a character string
that specifies which options are turned on/off or their value, in a unix-like
command line argument way.

The string must be composed of characters that correspond to a given option
(see mapping below), and modifiers '+' and '-' that toggle options on and off respectively.
E.g. \code{.options='tv'} will toggle on options \code{track} and \code{verbose},
while \code{.options='t-v'} will toggle on option \code{track} and toggle off
option \code{verbose}.

Modifiers '+' and '-' apply to all option character found after them:
\code{t-vp+k} means \code{track=TRUE}, \code{verbose=parallel=FALSE},
and \code{keep.all=TRUE}.
The default behaviour is to assume that \code{.options} starts with a '+'.

for options that accept integer values, the value may be appended to the
option's character e.g. \code{'p4'} for asking for 4 processors or \code{'v3'}
for showing verbosity message up to level 3.

The following options are available (the characters after \dQuote{-} are those
to use to encode \code{.options} as a string):
\describe{

\item{debug - d}{ Toggle debug mode (default: \code{FALSE}).
Like option \code{verbose} but with more information displayed.}

\item{keep.all - k}{ used when performing multiple runs (\code{nrun}>1): if
\code{TRUE}, all factorizations are saved and returned (default: \code{FALSE}).
Otherwise only the factorization achieving the minimum residuals is returned.}

\item{parallel - p}{ this option is useful on multicore *nix or Mac machine
only, when performing multiple runs (\code{nrun} > 1) (default: \code{TRUE}).
If \code{TRUE}, the runs are performed using the parallel foreach backend
defined in argument \code{.pbackend}.
If this is set to \code{'mc'} or \code{'par'} then \code{nmf} tries to
perform the runs using multiple cores with package
\code{link[doParallel]{doParallel}} -- which therefore needs to be installed.

If equal to an integer, then \code{nmf} tries to perform the computation on
the specified number of processors.
When passing options as a string the number is appended to the option's character
e.g. \code{'p4'} for asking for 4 processors.

If \code{FALSE}, then the computation is performed sequentially using the base
function \code{\link{sapply}}.

Unlike option 'P' (capital 'P'), if the computation cannot be performed in
parallel, then it will still be carried on sequentially.

\strong{IMPORTANT NOTE FOR MAC OS X USERS:} The parallel computation is
based on the \code{doMC} and \code{multicore} packages, so the same care
should be taken as stated in the vignette of \code{doMC}: \emph{\dQuote{it
is not safe to use doMC from R.app on Mac OS X. Instead, you should use doMC
from a terminal session, starting R from the command line.}} }

\item{parallel.required - P}{ Same as \code{p}, but an error is thrown if
the computation cannot be performed in parallel or with the specified number
of processors.}

\item{shared.memory - m}{ toggle usage of shared memory (requires the
\pkg{synchronicity} package).
Default is as defined by \code{nmf.getOption('shared.memory')}.}

\item{restore.seed - r}{ deprecated option since version 0.5.99.
Will throw a warning if used.}

\item{simplifyCB - S}{ toggle simplification of the callback results.
Default is \code{TRUE}}

\item{track - t}{ enables error tracking (default: FALSE).
If \code{TRUE}, the returned object's slot \code{residuals} contains the
trajectory of the objective values, which can be retrieved via
\code{residuals(res, track=TRUE)}
This tracking functionality is available for all built-in algorithms.
}

\item{verbose - v}{ Toggle verbosity (default: \code{FALSE}).
If \code{TRUE}, messages about the configuration and the state of the
current run(s) are displayed.
The level of verbosity may be specified with an integer value, the greater
the level the more messages are displayed.
Value \code{FALSE} means no messages are displayed, while value \code{TRUE}
is equivalent to verbosity level 1.
}

}}
  \item{alg}{One NMF algorithm to be used.}
  \item{maxIter}{Maximum of interactions to perform.}
  \item{nrun}{Number of runs to perform.}
}
\value{
An object of class \code{\linkS4class{NMFfit}} and files to be saved in the folder named after the \code{prefix}.
}
\author{
Cintia C Palu
}
\seealso{
\code{\link{nmfAlgorithm}}
\code{\link{exploreNMF}}
}
\examples{
n <- 20; counts <- c(5, 3, 2);
p <- sum(counts)
x <- syntheticNMF(n, counts)
dim(x)

# build the true cluster membership
groups <- unlist(mapply(rep, seq(counts), counts))

res.runNmf=runNmf(x, groups,r = 3,alg = "snmf/r", prefix='runNmf')
}

