# Test parallel computations
# 
# Author: Renaud Gaujoux
###############################################################################

# make the internal functions/objects visible
if( isNamespaceLoaded('NMF') ){
	setupLibPaths <- NMF:::setupLibPaths
	setupBackend <- NMF:::setupBackend
}

library(stringr)
library(foreach)

shared_DEACTIVATED <- function(...){
	msg <- NULL
	if( .Platform$OS.type == 'windows' ) msg <- str_c(..., ' [OS: Windows]')
	else if( !require.quiet(bigmemory) ) msg <- str_c(..., ' [missing: bigmemory]')
	else if( !require.quiet(synchronicity) ) msg <- str_c(..., ' [missing: synchronicity]')
	
	if( !is.null(msg) ) DEACTIVATED(msg)
}

check_shared_memory <- function(.msg, libs=TRUE, seq=FALSE){
	
	.test <- function(.msg, mutex, libs, seq){
		
		mess <- function(...){
			paste(.msg
				, if( mutex ) "With mutex" else "No mutex"
				, ":", ...)
		}
		
		mtx <- if( mutex ) ts_eval() else force
		
		if( libs ) setupLibPaths()
		
		alpha <- 5
		res <- foreach(i=1:4) %dopar% {
			t <- Sys.time()
			if( i==1 ) mtx(Sys.sleep(3))
			else if( i== 2) Sys.sleep(0.2)
			mtx({a <- runif(i); c <- 10 * i; d <- alpha + i})
			b <- c
			list(i, Sys.getpid(), t, Sys.time(), a, b, c, d)
		}
		
		pids <- sapply(res, '[[', 2)
		wtime <- sapply(res, function(x) round(as.numeric(x[[4]] - x[[3]]), 2))
		pid <- unique(pids)
		stopifnot( length(pid) == if( seq ) 1L else 2L )
		
		# check evaluation
		checkEquals( length(unlist(lapply(res,'[[', 5))), 4 *5 /2, mess("Evaluation of random draws is OK"))
		checkIdentical( sapply(res,'[[', 6), 1:4 * 10, mess("Evaluation outside eval call is OK"))
		checkIdentical( sapply(res,'[[', 7), 1:4 * 10, mess("Evaluation inside eval call is OK"))
		checkIdentical( sapply(res,'[[', 8), alpha + 1:4, mess("Evaluation inside eval call with exported variable is OK"))
		
		# return time differences
		ipid <- if( seq ) 1:2
				else c(which(pids == pid[1])[1L], which(pids == pid[2])[1L])
		wt <- wtime[ipid]
		
#		message(mess())
#		message( str_out(wt))
		wt
	}
	
	mess <- function(...) paste(.msg, ":", ...)
	
	# restore doSEQ backend on.exit
	on.exit( registerDoSEQ() )
		
	# no mutex
	wtime <- .test(mess(), mutex=FALSE, libs, seq)
	checkTrue( wtime[1] >= 2 , mess("No mutex: Thread 1 waits 2 second (", wtime[1], ')'))
	checkTrue( wtime[2] <  1 , mess("No mutex: Thread 2 does not wait at all (", wtime[2], ')'))
	
	# check mutex lock
	shared_DEACTIVATED("NMF shared memory feature not available.")
	
	wtime <- .test(mess(), mutex=TRUE, libs, seq)
	checkTrue( wtime[1] >= 2 , mess("With mutex : Thread 1 waits 2 seconds (", wtime[1], ')'))
	if( !seq )
		checkTrue( wtime[2] > 2 , mess("With mutex: Thread 2 also waits at least 2 seconds (", wtime[2], ')'))
	
}

test.shared_memory_doSEQ <- function(){
	
	# doSEQ
	registerDoSEQ()
	check_shared_memory('doSEQ', libs=FALSE, seq=TRUE)
	
}

test.shared_memory_doMC <- function(){
	# doParallel (doMC)
	library(doParallel)
	registerDoParallel(2)
	check_shared_memory('doParallel - Multicore', libs=FALSE)
	
}

test.shared_memory_doParallel <- function(){
	# doParallel (doParallel)
	cl <- makeCluster(2, outfile='wout.log')
	on.exit( stopCluster(cl), add=TRUE)
	registerDoParallel(cl)
	check_shared_memory('doParallel')
	
}

test.shared_memory_doMPI <- function(){
	DEACTIVATED("NMF shared memory feature does not currently work with doMPI.")	
	if( !require(doMPI) ) DEACTIVATED("Package doMPI not available.")
	# doMPI
	cl_MPI <- startMPIcluster(2)
	on.exit( closeCluster(cl_MPI), add=TRUE)
	registerDoMPI(cl_MPI)
	check_shared_memory('doMPI')
}

test.setupBackend <- function(){
	
	# restore backend on.exit
	on.exit( registerDoSEQ() )
	
	checkException( setupBackend(-1, 'par', TRUE), "Invalid number of cores (optional)")
	checkException( setupBackend(-1, 'par', FALSE), "Invalid number of cores (required)")
	checkException( setupBackend(10, 'par', FALSE), "Required too many cores")
	checkException( setupBackend(1, 'toto', FALSE), "Required unknown backend")
	
}


test.gVariable <- function(){
	
	# restore backend on.exit
	on.exit( registerDoSEQ() )
	
	.check <- function(.msg, libs=TRUE, seq=FALSE){
		on.exit( registerDoSEQ() )
		
		.test <- function(shared){
			mess <- function(...) paste(.msg, ' + shared=', shared, ":", ...)
            cat(mess(), "\n")
			
			# run foreach loop
			v <- gVariable(123, shared=shared)
			if( libs ) setupLibPaths(verbose = TRUE)
			res <- foreach(i=1:20) %dopar% { 
				if(i==1) v(456) else if( i== 2) Sys.sleep(0.2); c(Sys.getpid(), v())
			}
			
			# extract result data
			pids <- sapply(res, '[', 1)
			vals <- sapply(res, '[', 2)
			pid <- unique(pids)
			stopifnot( length(pid) == if( seq ) 1L else 2L )
			
			# when not shared: only the iterations run by the first process see changes
			if( !shared && !seq ){
				checkIdentical( unique(vals[pids==pid[1]]), 456, mess("Value change in first process affects first process"))
				checkIdentical( unique(vals[pids==pid[2]]), 123, mess("Value change ins first process does not affect second process"))
			}
			else{
				checkIdentical( unique(vals), 456
					, mess("Value change affects all processes"))
			}
		}
		.test(FALSE)
		
		shared_DEACTIVATED("NMF global shared variables not available.")
		.test(TRUE)
	}
	
	# doSEQ
	registerDoSEQ()
	.check('doSEQ', libs=FALSE, seq=TRUE)
	
	# doParallel (Multicore)
	library(doParallel)
	registerDoParallel(2)
	.check('doParallel - Multicore')
	
	# doParallel (doSNOW)
	cl <- makeCluster(2, outfile='wout.log')
	on.exit( stopCluster(cl), add=TRUE)
	registerDoParallel(cl)
	.check('doParallel')
		
	# doMPI
	if( !require(doMPI) ) DEACTIVATED("Package doMPI not available.")
	cl_MPI <- startMPIcluster(2)
	on.exit( closeCluster(cl_MPI), add=TRUE)
	registerDoMPI(cl_MPI)
	.check('doMPI')
	
}

test.ForeachBackend <- function(){
	
	.check <- function(type, n, ...){
		b <- ForeachBackend(...)
		checkIdentical(class(b), c(str_c(type, '_backend'), 'foreach_backend'), str_c(type, ": Class is ok"))
		b
	}
	
	# doParallel (Multicore)
	library(doParallel)
	.check('doParallel', 3, 'PAR', 3)
	
	# doParallel (SNOW)
	cl <- makeCluster(2)
	on.exit( stopCluster(cl), add=TRUE)
	b <- .check('doParallel', 2, cl)
	
	# doMPI
	if( !require(doMPI) ) DEACTIVATED("Package doMPI not available.")
	b <- .check('doMPI', 2, 'MPI', 2)
	cl_MPI <- startMPIcluster(2)
	on.exit( closeCluster(cl_MPI), add=TRUE)
	b <- .check('doMPI', 2, cl_MPI)
	
}

test.nmf <- function(){
	
	on.exit( registerDoSEQ() )
	
	set.seed(123456)
	a <- rmatrix(20,10)
	nmf.options(cores=2)
	
	checkTrue( isNMFfit(resREF <- nmf(a, 2, seed=123, nrun=3, .opt='v3')), "Default works")
	
	cl_loadedNamespaces <- function(cl=NULL){
		if( is_NA(cl) ) return()
		if( is.null(cl) ) unique(unlist(foreach(i=1:2) %dopar% { loadedNamespaces() }))
		else unique(unlist(clusterApplyLB(cl, 1:2, function(i){ loadedNamespaces() })))
	}
	
	.check <- function(msg, .options=NULL, ..., LOADED_NAMESPACES=NA){
		
		be <- getDoBackend()
		checkTrue( isNMFfit(res2 <- nmf(a, 2, seed=123, nrun=3, .opt=str_c('v3', .options), ...)), str_c(msg, " works"))
		# retrieve namespaces on dopar processes if necessary
		ns <- if( !is_NA(LOADED_NAMESPACES) ) cl_loadedNamespaces(LOADED_NAMESPACES)
		
		checkTrue( nmf.equal(resREF, res2), str_c(msg, ": result is identical to default") )
		checkIdentical( consensus(resREF, no.attrib=TRUE), consensus(res2, no.attrib=TRUE)
					, str_c(msg, ": consensus matrice (no.attrib) is identical to default") )
		checkIdentical( consensus(resREF), consensus(res2), str_c(msg, ": consensus matrice is identical to default") )
		checkTrue( identical(be, getDoBackend()), str_c(msg, ": backend is restored") )
		
		# check restoration on error
		checkException( nmf(a, 2, method=function(...) 1L, seed=123, nrun=3, .opt=str_c('v3', .options), ...), str_c(msg, " throw error if bad method"))
		checkTrue( identical(be, getDoBackend()), str_c(msg, ": backend is restored after error") )
		
		# return loaded namespaces
		ns
	}
	
	library(parallel)
	
	.check('SEQ', .pbackend='SEQ')
	# Multicore
	if( parallel::detectCores() > 1 ) 
		.check('P2', .options='P2')
	
	# SNOW-type from char spec
	.check('.pbackend="psock"', .options='P2', .pbackend='PSOCK')
	
	# SNOW-type
	cl <- makeCluster(2)
	on.exit( stopCluster(cl), add=TRUE)
	lpkg0 <- cl_loadedNamespaces(cl)
	lpkg1 <- .check('.pbackend=cl + SNOW-like cluster', .pbackend=cl, LOADED_NAMESPACES=cl)
#	if( !isCHECK() ){ #TODO: investigate why this test does not work in R CMD check
#		checkTrue(length(setdiff(lpkg1, lpkg0)) > 0, "Provided cluster was used if .pbackend=cl")
#	}
	stopCluster(cl)
	
	library(doParallel)
	cl <- makeCluster(2)
	registerDoParallel(cl)
	lpkg0 <- cl_loadedNamespaces(NULL)
	# no .pbackend => use registered cluster
	lpkg1 <- .check('doParallel registered cluster + P2 [should not use registered cluster]', .opt='P2'
					, LOADED_NAMESPACES=NULL)
	checkIdentical(lpkg0, lpkg1, "Registered cluster was not used if .opt='P2'")
	# .pbackend=NULL
	lpkg1 <- .check('.pbackend=NULL + doParallel registered cluster', .pbackend=NULL
					, LOADED_NAMESPACES=NULL)
	checkTrue(length(setdiff(lpkg1, lpkg0))>0, "Registered cluster was used if .pbackend=NULL")
	
	# MPI
	if( !require(doMPI) ) DEACTIVATED("Package doMPI not available.")
	cl_MPI <- startMPIcluster(2)
	on.exit( closeCluster(cl_MPI), add=TRUE)
	.check('.pbackend=cl_MPI + MPI cluster', .pbackend=cl_MPI)
	registerDoMPI(cl_MPI)
	.check('.pbackend=NULL + doMPI registered MPI cluster', .pbackend=NULL)	
}
