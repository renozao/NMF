# Run all unit tests in installed directory unitTests
# 
# Author: Renaud Gaujoux
# Creation: 26 Oct 2011
###############################################################################

tests <- try( pkgmaker::utest('package:NMF', quiet=FALSE) )

td <- pkgmaker:::utestPath(package='package:NMF')
resfile <- list.files(td, pattern=".+\\.txt", full.names=TRUE)
cat("Result files:\n")
print(resfile)

if( length(resfile) ){
	# send
	library(mail)
	sapply(resfile, function(f){
		sendmail('renaud@cbio.uct.ac.za',
				paste("NMF: unit test results",
					"-", basename(f),
					"-", if( is(tests, 'try-error') ) 'ERROR' else "OK", "]"),
				paste(readLines(f), collapse="\n")
		)
	})
}
