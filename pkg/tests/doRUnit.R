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
				
		# build message
		msg <- c("**************\nR.version Info\n**************\n", capture.output(R.version))
		sys <- Sys.info()
		msg <- c(msg, "**************\nSystem Info\n**************\n"
				, sapply(names(sys), function(n){ paste(n, ': ', sys[n], sep='')}))
		msg <- c(msg, "**************\nRESULTS:\n**************\n", readLines(f))
		# collapse
		msg <- paste(msg, collapse="\n")
		# send email
		sendmail('renaud@cbio.uct.ac.za',
				, paste("NMF: unit test results",
					"-", basename(f),
					"-", if( is(tests, 'try-error') ) 'ERROR' else "OK", "]")
				, msg 
		)
	})
}
