# Description and generation of data
# 
# Author: Renaud Gaujoux
###############################################################################

#' Golub ExpressionSet
#' 
#' This data comes originally from the gene expression data from \cite{Golub1999}. 
#' The version included in the package is the one used and referenced in \cite{Brunet2004}.
#' The samples are from 27 patients with acute lymphoblastic leukemia (ALL) and 
#' 11 patients with acute myeloid leukemia (AML).
#' 
#' The samples were assayed using Affymetrix Hgu6800 chips and the original
#' data on the expression of 7129 genes (Affymetrix probes) are available on
#' the Broad Institute web site (see references below).
#' 
#' The data in \code{esGolub} were obtained from the web page related to 
#' the paper from \cite{Brunet2004}, which describes an application of 
#' Nonnegative Matrix Factorization to gene expression clustering.
#' (see link in section \emph{Source}).
#' 
#' They contain the 5,000 most highly varying genes according to their
#' coefficient of variation, and were installed in an object of class
#' \code{\link[Biobase]{ExpressionSet-class}}.
#' 
#' @format There are 3 covariates listed.
#' 
#' \itemize{
#' 
#' \item Samples: The original sample labels.  \item ALL.AML: Whether the
#' patient had AML or ALL. It is a \code{\link{factor}} with levels
#' \code{c('ALL', 'AML')}.  \item Cell: ALL arises from two different types of
#' lymphocytes (T-cell and B-cell).  This specifies which for the ALL patients;
#' There is no such information for the AML samples. It is a
#' \code{\link{factor}} with levels \code{c('T-cell', 'B-cell', NA)}.
#' 
#' }
#' 
#' @source
#' Web page for \cite{Brunet2004}:\cr
#' \url{http://www.broadinstitute.org/publications/broad872}
#'  
#' Original data from Golub et al.:\cr
#' \url{http://www-genome.wi.mit.edu/mpr/data_set_ALL_AML.html}
#' 
#' @name esGolub
#' @docType data
#' @keywords datasets
#' @examples
#' 
#' # requires package Biobase to be installed
#' if( require(Biobase) ){
#' 
#' 	data(esGolub)
#' 	esGolub
#' 	\dontrun{pData(esGolub)}
#' 
#' }
#'  
NULL

