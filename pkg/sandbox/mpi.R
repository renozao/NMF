# Running NMF in parallel using MPI and doMPI
# 
# Author: Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###############################################################################

## 0. Create and register an MPI cluster
library(doMPI)
cl <- startMPIcluster()
registerDoMPI(cl)

lib.dir <- 'lib'
error <- try({ #START_TRY

## 1. Setup MPI options
initEnvir <- function(envir, lib.dir) {

	.libPaths(c(lib.dir, .libPaths()))
	# init log file	
	fname <- tempfile(paste('log', mpi.comm.rank(), Sys.info()['nodename'], '', sep='_'), tmpdir=getwd())
	# open logfile and redirect output
	envir$f <- file(fname, "a+")	
	sink(envir$f, append=TRUE)
	sink(envir$f, append=TRUE, type='message')

}

initArgs <- list(lib.dir=lib.dir)

finalEnvir <- function(envir) {
  	# clean-up and close log file
	sink(type="message")
	sink()	
	close(envir$f)
}


mpiopts <- list(initEnvir=initEnvir, initArgs=initArgs, finalEnvir=finalEnvir)

## 2. Schedule the runs accross the workers
nrun <- 100;
nworker <- getDoParWorkers();
ntasks <- rep(round(nrun/nworker), nworker)
if( (remain <- nrun %% nworker) > 0 ) 
	ntasks[1:remain] <- ntasks[1:remain] + 1

message("# :: MPI-nmf :: Use ", getDoParWorkers(), " workers for ", nrun , " runs")
t <- system.time({

	res <- foreach(i=1:getDoParWorkers(), n=ntasks, .verbose=TRUE
			, .options.mpi=mpiopts
			, .packages = c('NMF', 'doMC', 'Biobase')) %dopar% {
	
		cat('# :: MPI-dopar :: Task set number ', i, " -> ", n, " tasks\n")	
		# load Golub dataset
		data(esGolub)
		# run NMF
		nmf(esGolub, 3, 'brunet', nrun=n, .opt='v')
	}

})

## 3. save the result in a file
save(res, t, file='result.RData')

message("Done")

}) #END_TRY

if( is(error, 'try-error') ) message("Error in MPI run:", error)

## 4. Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
