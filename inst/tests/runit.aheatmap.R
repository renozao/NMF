# Unit tests for the heatmap drawing function
# 
# Author: Renaud Gaujoux
# Creation: 18 Nov 2011
###############################################################################

checkPlot <- function(...){
	if( isCHECK() ) return()
	pkgmaker::checkPlot(...)
}

test.mfrow <- function(){
	
	x <- rmatrix(20,10)
	
	checkPlot({
	op <- par(mfrow=c(1,2))
	on.exit(par(op))	
	aheatmap(x)
	aheatmap(x*100)
	}, "Using mfrow correctly generates two heatmaps side by side.")
}
