# Utility functions to work with octave
# 
# Author: Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###############################################################################

# Load library 'foreign for octave file reading
library(foreign)

save.ascii <- function(X, file, digits=11, quote=FALSE, ...){
	
	X <- format(X, digits=digits)
	write.table(X, file=file, row.names=FALSE, col.names=FALSE, quote=quote, ...)
	
	return(invisible(X))
}

