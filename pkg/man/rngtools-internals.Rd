\name{RNG-Internals}
\alias{RNG-internals}

\alias{initialize,rstream2-method}
\alias{initialize,rstream.user-method}
\alias{RNGinfo}
\alias{RNGhooks}
\alias{RNGlib}
\alias{RNGlibs}
\alias{RNGprovider}
\alias{RNGrecovery}
\alias{RNGseed}
\alias{RNGseq}
\alias{RNGwrap}

\alias{setRNG}
\alias{setRNG-methods}
\alias{setRNG,ANY-method}
\alias{setRNG,character-method}
\alias{setRNG,numeric-method}
\alias{setRNG,rstream-method}

\title{Undocumented Internal Methods and Functions for RNG Handling}

\description{
The functions documented here used internally by the package NMF to manage RNG 
settings.
}

\seealso{
	\code{\link[NMF]{getRNG-methods}} and the package's vignette
}

\author{Renaud Gaujoux}
\keyword{internal}

\examples{

# summary of the RNG settings
RNGinfo()

# get active RNG library
RNGlib()
RNGlib(full=TRUE)
# get info about a given RNG library 
RNGlib('rstream')

# get all RNG library (except the RNG wrapper)
RNGlibs()
# get all RNG library (including the RNG wrapper)
RNGlibs(all=TRUE)
RNGlibs(all=TRUE, full=TRUE)

# get current RNGprovider
RNGprovider()

r <- new('rstream.mrg32k3a')
setRNG(r)
\dontrun{
# could have done this as well
setRNG('rstream') 
}
# to see messages about the whole process
setRNG(r, verbose=TRUE)

RNGprovider()
RNGinfo()

.require <- require
.library <- library

if( .require(rlecuyer) ){
	
	# activate RNG library from rstream: use wrapper because rlecuyer hides 
	# rstream 	
	RNGwrap('rstream')
	RNGlibs(all=TRUE)
	
	# activate RNG library from rlecuyer: do not use wrapper because there is no 
	# issue in directly using the library
	RNGwrap('rlecuyer')
	RNGlibs(all=TRUE)
	
	detach('package:rlecuyer', unload=TRUE)
	
	# load another RNG library that provides other hooks than user_unif_rand	
	if( .require(randtoolbox) ){
		
		RNGhooks('randtoolbox')
		RNGhooks(who=TRUE)
	
		# rlecuyer only provides user_unif_rand
		.library(rlecuyer)
		RNGhooks('rlecuyer')
			
		# activate RNG library from rlecuyer: here one must use the wrapper 
		# library to ensure that the hooks from randtoolbox are not called
		RNGwrap('rlecuyer')
		RNGlibs(all=TRUE)
		
		detach('package:randtoolbox', unload=TRUE)
	}
	
	detach('package:rlecuyer', unload=TRUE)
	
	# activate RNG library from rstream: no need of wrapper as it is the only 
	# RNG library 	
	RNGwrap('rstream')
	RNGlibs(all=TRUE)
	
}

}
