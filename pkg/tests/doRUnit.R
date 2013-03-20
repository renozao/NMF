# Run all unit tests in installed directory unitTests
# 
# Author: Renaud Gaujoux
# Creation: 26 Oct 2011
###############################################################################

library(pkgmaker)
tests <- try( utest('package:NMF', quiet=FALSE) )

testdir <- pkgmaker:::utestPath(package='package:NMF')

resfile <- list.files(testdir, pattern=".+\\.txt", full.names=TRUE)
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
		# subject
		subject <- paste("Package NMF: unit test results"
						, "-", basename(f), "-"
						, "[", if( is(tests, 'try-error') ) 'ERROR' else "OK", "]"
						, sep='')
		if( isCRANcheck() ){
			subject <- paste('CRAN check -', subject)
		}
		# try send email
		if( !userIs('renaud') ) try( sendmail('renaud@cbio.uct.ac.za', subject, msg) )
		else write(msg, file=file.path(testdir, paste("check_", basename(f), sep='')))
	})
}

# stop if testing locally
if( userIs('renaud') && is(tests, 'try-error') ){
	stop(tests)
}
