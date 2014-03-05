# Original Matlab algorithm from Brunet et al. (2004)
# 
# Author: Renaud Gaujoux
# Created: 23 Nov 2012
###############################################################################

#' @include registry-algorithms.R
NULL

#' The algorithm \sQuote{.M#brunet} provide access to the \emph{original} Matlab 
#' algorithm from \cite{Brunet2004}, through the \pkg{RcppOctave} package.
#' The Matlab code used can be found in the \pkg{NMF} package's 'm-files/' 
#' sub-directory:
#' 
#' \samp{
#' library(RcppOctave)
#' file.show(system.mfile('brunet.m', package='NMF'))
#' }
#' 
#' @rdname KL-nmf
#' @aliases brunet_M-nmf
nmfAlgorithm.brunet_M <- setNMFMethod('.M#brunet',
	, objective= 'KL'
	, mcode = 'brunet.m'
	, algorithm = function (y, x, verbose=nmf.getOption('verbose')
                            , maxIter = nmf.getOption('maxIter') %||% 2000L){
		RcppOctave::.CallOctave('brunet', y, nbasis(x), verbose, basis(x), coef(x), maxIter);
	}
)
