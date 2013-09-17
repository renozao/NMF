# Functions produce reports   
# 
# Author: Renaud Gaujoux
# Created: 23 Jul 2013
###############################################################################


#' Run NMF Methods and Generate a Report
#' 
#' Generates an HTML report from running a set of method on a given 
#' target matrix, for a set of factorization ranks.
#' 
#' The report is based on an .Rmd document \code{'report.Rmd'} stored in 
#' the package installation sub-directory \code{scripts/}, and is compiled
#' using \pkg{knitr}.
#'  
#' At the begining of the document, a file named \code{'functions.R'} is 
#' looked for in the current directory, and sourced if present.
#' This enables the definition of custom NMF methods (see \code{\link{setNMFMethod}}) 
#' or setting global options.   
#' 
#' @param x target matrix
#' @param rank factorization rank
#' @param method list of methods to apply
#' @param colClass reference class to assess accuracy
#' @param ... extra paramters passed to \code{\link{nmf}}
#' @param output output HTML file
#' @param template template Rmd file
#' 
#' @return a list with the following elements:
#' \item{fits}{the fit(s) for each method and each value of the rank.}
#' \item{accuracy}{a data.frame that contains the summary assessment measures, 
#' for each fit.} 
#' 
#' @export
#' @examples 
#' 
#' \dontrun{
#' 
#' x <- rmatrix(20, 10)
#' gr <- gl(2, 5)
#' nmfReport(x, 2:4, method = list('br', 'lee'), colClass = gr, nrun = 5)
#' 
#' }
nmfReport <- function(x, rank, method, colClass = NULL, ..., output = NULL, template = NULL){
	
	library(knitr)
	if( is.null(template) )
		template <- system.file('scripts/report.Rmd', package = 'NMF')
	x <- force(x)
	rank <- force(rank)
	method <- force(method)
	args <- list(...)
	nmfRun <- function(x, rank, method, ...){
		args <- expand_dots(args)
		str(args)
		do.call(nmf, c(list(x, rank, method), args))
	}
    accuracy <- NA
    res <- NA
	knit2html(template)
	res <- list(fits = res, accuracy = accuracy)
	saveRDS(res, file = 'report_results.rds')
	invisible(res)
}
