# Run all unit tests in installed directory unitTests
# 
# Author: Renaud Gaujoux
# Creation: 26 Oct 2011
###############################################################################

tests <- try( pkgmaker::utest('package:NMF', quiet=FALSE) )

td <- utestPath('package:NMF')
resfile <- list.files(td, pattern=".+\\.html", full.names=TRUE)

library(mail)
sendmail('renaud@cbio.uct.ac.za',
		paste("NMF: unit test results",
			"[", if( is(tests, 'try-error') ) 'ERROR' else "OK", "]"),
		paste(readLines(resfile), collapse="\n")
)
