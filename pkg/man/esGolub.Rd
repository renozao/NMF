 \name{esGolub}
 \docType{data}
 \alias{esGolub}
 \title{Golub ExpressionSet from Brunet et al. Paper}
 \description{
   	The original data is related to Golub et al., and this version is the one 
   	used and referenced in Brunet et al.
   	The samples are from 27 patients with acute lymphoblastic leukemia (ALL) 
   	and 11 patients with acute myeloid leukemia (AML).
   	 
   	The samples were assayed using Affymetrix Hgu6800 chips and the original data 
   	on the expression of 7129 genes (Affymetrix probes) are available 
   	on the Broad Institute web site (see references below). 
   	
   	The data in \code{esGolub} were obtained from the web site related to 
   	Brunet et al.'s publication on an application of Nonnegative Matrix 
   	Factorization (see link in section \emph{Source}). 
   	
   	They contain the 5,000 most highly varying genes according to their 
   	coefficient of variation, and were installed in an object of class 
   	\code{\link[Biobase]{ExpressionSet-class}}.
 }
 \usage{data(esGolub)}
 \format{
 There are 3 covariates listed. 
 
\itemize{

\item Samples: The original sample labels. 
\item ALL.AML: Whether the patient had AML or ALL. It is a \code{\link{factor}} 
with levels \code{c('ALL', 'AML')}. 
\item Cell: ALL arises from two different types of lymphocytes (T-cell and B-cell). 
This specifies which for the ALL patients; There is no such information for 
the AML samples. It is a \code{\link{factor}} with levels 
\code{c('T-cell', 'B-cell', NA)}.
 
 }
 
 }
 \source{\url{http://www.broadinstitute.org/publications/broad872}}
 \references{
 
 	Brunet, J.P., Tamayo, P., Golub, T.R., and Mesirov, J.P. (2004)
 	\emph{Metagenes and molecular pattern discovery using matrix factorization}.	
	Proc Natl Acad Sci USA 101(12), 4164--4169.
	
	T. R. Golub et al. (1999)
	\emph{Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression Monitoring}.
 	Science, 531-537, 1999
 	 	
 	Original data from Golub et al.:\cr
 	\url{http://www-genome.wi.mit.edu/mpr/data_set_ALL_AML.html}  
 	 	   
 }
 \examples{
 	data(esGolub)
 	esGolub
 	\dontrun{pData(esGolub)}
 }
 \keyword{datasets}

