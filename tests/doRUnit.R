# Run all unit tests in installed directory unitTests
# 
# Author: Renaud Gaujoux
# Creation: 26 Oct 2011
###############################################################################

library(pkgmaker)

# Skip checks except if run locally    
if( !isFALSE(Sys.getenv_value('_R_LOCAL_CHECK_')) ){
    
# skip tests on CRAN checks
#if( !isCRANcheck() ){
library(NMF)
library(RUnit)
nmf.options(maxIter=100L)
    
tests <- try( utest('package:NMF', quiet=FALSE) )

if( FALSE ){
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

} # end if NOT CRAN check

# stop if error
if( is(tests, 'try-error') ){
	stop(tests)
}

}
