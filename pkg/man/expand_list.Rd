\name{expand_list}
\alias{expand_dots}
\alias{expand_list}
\title{Expanding Lists}
\usage{
  expand_list(x, ..., .exact = TRUE, .names = !.exact)

  expand_dots(..., .exclude = NULL)
}
\arguments{
  \item{x}{input list}

  \item{...}{extra named arguments defining the default
  items. A list of default values can also be passed as a a
  single unnamed argument.}

  \item{.exact}{logical that indicates if the names in
  \code{x} should be partially matched against the
  defaults.}

  \item{.names}{logical that only used when
  \code{.exact=FALSE} and indicates that the names of items
  in \code{x} that partially match some defaults should be
  expanded in the returned list.}

  \item{.exclude}{optional character vector of argument
  names to exclude from expansion.}
}
\value{
  a list
}
\description{
  \code{expand_list} expands a named list with a given set
  of default items, if these are not already in the list,
  partially matching their names.

  \code{expand_dots} expands the \code{...} arguments of
  the function in which it is called with default values,
  using \code{expand_list}. It can \strong{only} be called
  from inside a function.
}
\examples{
expand_list(list(a=1, b=2), c=3)
expand_list(list(a=1, b=2, c=4), c=3)
# with a list
expand_list(list(a=1, b=2), list(c=3, d=10))
# no partial match
expand_list(list(a=1, b=2, c=5), cd=3)
# partial match with names expanded
expand_list(list(a=1, b=2, c=5), cd=3, .exact=FALSE)
# partial match without expanding names
expand_list(list(a=1, b=2, c=5), cd=3, .exact=FALSE, .names=FALSE)

# works also inside a function to expand a call with default arguments
f <- function(...){
	cl  <- match.call()
	expand_list(cl, list(a=3, b=4), .exact=FALSE)
}
f()
f(c=1)
f(a=2)
f(c=1, a=2)
# expanding dot arguments

f <- function(...){
	expand_dots(list(a=2, bcd='a', xxx=20), .exclude='xxx')
}

# add default value for all arguments
f()
# add default value for `bcd` only
f(a=10)
# expand names
f(a=10, b=4)
}

