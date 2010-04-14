\name{advanced}
\Rdversion{1.1}
\alias{NMF-advanced}
\alias{nmf.options}
\alias{nmf.getOption}
\alias{nmf.options.reset}
\alias{nmf.options.runtime}

\title{ Advanced usage of package NMF }
\description{
	Allow the user to get/set/define package NMF specific options in the same way as with
	base functions \code{\link[base]{options}} and \code{\link[base]{getOption}}.
	
}
\usage{

nmf.options(..., runtime = FALSE)
nmf.getOption(name)
nmf.options.reset()
nmf.options.runtime()

}
%- maybe also 'usage' for other objects documented here.
\arguments{

	\item{\dots}{ any options can be defined, using 'name = value' or by
	   passing a list of such tagged values.  However, only the ones below are 
	   used in package NMF. Further, 'nmf.options('name') == nmf.options()['name']', 
	   see the example.
	}
		
	\item{name}{ a character string holding an option name.
	}
  	
  	\item{runtime}{ a boolean used to specify if main interface function \code{\link{nmf}} 
	  should store the option into the initial \code{\linkS4class{NMF}} object before performing
	  the computation.
	}
		
}
\details{

\describe{

	\item{nmf.getOption}{ get the value of a single option.
	}
		
	\item{nmf.options}{ gets/sets/defines options in the way of base function \code{\link[base]{options}}.
	Invoking 'nmf.options()' with no arguments returns a list with the current values of the options. 
	To access the value of a single option, one should use \code{nmf.getOption("error.track")}, e.g., 
	rather than \code{nmf.options("error.track")} which is a \code{list} of length one.	
	}
		
	\item{nmf.options.reset}{ Reset all \emph{built-in options} to their default values. 
	Note that only built-in are reset. The options defined by the user during the current session 
	will keep their values.	
	}
	
}
}
\value{
	For \code{nmf.getOption}, the current value set for option \code{name}, or \code{NULL}
    if the option is unset.
	
    For \code{nmf.options()}, a list of all set options sorted by name.  For
    \code{options(name)}, a list of length one containing the set value, or
    \code{NULL} if it is unset.  For uses setting one or more options, a
    list with the previous values of the options changed (returned
    invisibly).
        
}

\section{Options set in package NMF}{

\describe{

\item{debug}{ logical. Similar to option \code{'verbose'} (see below), but reports more information.
}
\item{default.algorithm}{ \code{character}. The default NMF algorithm used by the \code{\link{nmf}} 
method when called without argument \code{method}.
}
\item{default.seed}{ \code{character}. The default seeding method used by the \code{\link{nmf}} 
method when called without argument \code{seed}.
}
\item{error.track}{ logical. Should the estimation error be tracked during the computations?
If set to \code{TRUE} then the error track can be plotted using method \code{\link{plot}}. 
The step size of the error track is set via option \code{track.interval} (see below).
}
\item{parallel.backend}{ \code{character} or \code{numeric}. 
The default parallel backend used by the \code{\link{nmf}} 
method when called with argument \code{nrun} greater than 1.

Currently it accepts the following values: \code{'mc'} or a number that specifies 
the number of cores to use, \code{'seq'} or \code{NULL} to use sequential backend for \code{foreach}, 
and the empty string \code{''} to completely disable the parallel computation. 

}
\item{track.interval}{ numeric. The number of iterations performed between two consecutive error points. 
For performance reason, this value should be too small, as the computation of the estimation error can be  
time consuming (Default value is 30).
}
\item{verbose}{ logical. Should R report extra information about the computations?
}


}

}

\author{ Renaud Gaujoux \email{renaud@cbio.uct.ac.za} }

\seealso{
	\code{\link{options}}
}

\examples{

	# save all options value 
	op <- nmf.options(); 
	utils::str(op) # op may contain functions.

    nmf.getOption("track.interval") == nmf.options()$track.interval # the latter needs more memory
    
    x <- matrix(runif(50*10), 50, 10) # create a random target matrix
    # or define a synthetic data with a hidden pattern using function syntheticNMF (see ?syntheticNMF) 
    \dontrun{x <- syntheticNMF(50, 5, 10, noise=TRUE)}
    
    # perform default NMF computation
    res <- nmf(x, 3)
    
    # Toogle on verbose mode
    nmf.options(verbose = TRUE)    
    res <- nmf(x, 3)    

	# Toogle on debug mode
    nmf.options(debug = TRUE)    
    res <- nmf(x, 3)

    # set the error track step size, and save previous value
    old.o <- nmf.options(track.interval = 5)
    old.o
    
    # check options
    utils::str(nmf.options())
    # reset to default values
    nmf.options.reset()
    utils::str(nmf.options())
	
}
