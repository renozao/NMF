# Test parallel computations
# 
# Author: Renaud Gaujoux
# Converted from RUnit on 22 Feb 2020
###############################################################################

# Setup ----
nmf.options(maxIter = 50L)
on.exit(nmf.options(maxIter = NULL), add = TRUE)

# Utils ----
str_c <- paste0
shared_DEACTIVATED <- function(...){
    msg <- NULL
    if( .Platform$OS.type == 'windows' ) msg <- str_c(..., ' [OS: Windows]')
    else if( !require.quiet(bigmemory) ) msg <- str_c(..., ' [missing: bigmemory]')
    else if( !require.quiet(synchronicity) ) msg <- str_c(..., ' [missing: synchronicity]')

    if( !is.null(msg) ) skip(msg)
}

check_shared_memory <- function(.msg, libs=TRUE, seq=FALSE){

    .test <- function(.msg, mutex, libs, seq){

        mess <- function(...){
            paste(.msg
                  , if( mutex ) "With mutex" else "No mutex"
                  , ":", ...)
        }

        mtx <- if( mutex ) ts_eval() else force

        if( libs ) NMF:::setupLibPaths()

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
        expect_equal( length(unlist(lapply(res,'[[', 5))), 4 *5 /2, info = mess("Evaluation of random draws is OK"))
        expect_identical( sapply(res,'[[', 6), 1:4 * 10, mess("Evaluation outside eval call is OK"))
        expect_identical( sapply(res,'[[', 7), 1:4 * 10, mess("Evaluation inside eval call is OK"))
        expect_identical( sapply(res,'[[', 8), alpha + 1:4, mess("Evaluation inside eval call with exported variable is OK"))

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
    expect_true( wtime[1] >= 2 , mess("No mutex: Thread 1 waits 2 second (", wtime[1], ')'))
    expect_true( wtime[2] <  1 , mess("No mutex: Thread 2 does not wait at all (", wtime[2], ')'))

    # check mutex lock
    shared_DEACTIVATED("NMF shared memory feature not available.")

    wtime <- .test(mess(), mutex=TRUE, libs, seq)
    expect_true( wtime[1] >= 2 , mess("With mutex : Thread 1 waits 2 seconds (", wtime[1], ')'))
    if( !seq )
        expect_true( wtime[2] > 2 , mess("With mutex: Thread 2 also waits at least 2 seconds (", wtime[2], ')'))

}


# Tests ----
test_that("test.ForeachBackend", {
    .check <- function(type, n, ...) {
        b <- ForeachBackend(...)
        expect_identical(class(b), c(str_c(type, "_backend"), "foreach_backend"),
            str_c(type, ": Class is ok"))
        b
    }
    library(doParallel)
    .check("doParallel", 3, "PAR", 3)
    cl <- makeCluster(2)
    on.exit(stopCluster(cl), add = TRUE)
    b <- .check("doParallel", 2, cl)
    if (!require.quiet("doMPI"))
        skip("Package doMPI not available.")
    skip("doMPI checks are disabled due to issues in doMPI::closeCluster")
    b <- .check("doMPI", 2, "MPI", 2)
    cl_MPI <- startMPIcluster(2)
    on.exit(closeCluster(cl_MPI), add = TRUE)
    b <- .check("doMPI", 2, cl_MPI)
})

test_that("test.gVariable", {
    on.exit(registerDoSEQ())
    .check <- function(.msg, libs = TRUE, seq = FALSE) {
        on.exit(registerDoSEQ())
        .test <- function(shared) {
            mess <- function(...) paste(.msg, " + shared=", shared,
                ":", ...)
            cat(mess(), "\n")
            v <- gVariable(123, shared = shared)
            if (libs)
                NMF:::setupLibPaths(verbose = TRUE)
            res <- foreach(i = 1:20) %dopar% {
                if (i == 1)
                  v(456)
                else if (i == 2)
                  Sys.sleep(0.2)
                c(Sys.getpid(), v())
            }
            pids <- sapply(res, "[", 1)
            vals <- sapply(res, "[", 2)
            pid <- unique(pids)
            stopifnot(length(pid) == if (seq)
                1L
            else 2L)
            if (!shared && !seq) {
                expect_identical(unique(vals[pids == pid[1]]),
                  456, mess("Value change in first process affects first process"))
                expect_identical(unique(vals[pids == pid[2]]),
                  123, mess("Value change ins first process does not affect second process"))
            }
            else {
                expect_identical(unique(vals), 456, mess("Value change affects all processes"))
            }
        }
        .test(FALSE)
        shared_DEACTIVATED("NMF global shared variables not available.")
        .test(TRUE)
    }
    registerDoSEQ()
    .check("doSEQ", libs = FALSE, seq = TRUE)
    library(doParallel)
    registerDoParallel(2)
    .check("doParallel - Multicore")
    cl <- makeCluster(2, outfile = "wout.log")
    on.exit(stopCluster(cl), add = TRUE)
    registerDoParallel(cl)
    .check("doParallel")
    if (!require.quiet("doMPI"))
        skip("Package doMPI not available.")
    skip("doMPI checks are disabled due to issues in doMPI::closeCluster")
    cl_MPI <- startMPIcluster(2)
    on.exit(closeCluster(cl_MPI), add = TRUE)
    registerDoMPI(cl_MPI)
    .check("doMPI")
})

test_that("test.nmf", {
    on.exit(registerDoSEQ())
    set.seed(123456)
    a <- rmatrix(20, 10)
    nmf.options(cores = 2)
    
    # default run (mutlticore)
    expect_true(isNMFfit(resREF <- nmf(a, 2, seed = 123, nrun = 2, .opt = "v3")), info = "Default works")
    # cl_loadedNamespaces <- function(cl = NULL) {
    #     if (is_NA(cl)) 
    #         return()
    #     if (is.null(cl)) 
    #         unique(unlist(foreach(i = 1:2) %dopar% {
    #             loadedNamespaces()
    #         }))
    #     else unique(unlist(clusterApply(cl, seq(length(cl)), function(i) {
    #         loadedNamespaces()
    #     })))
    # }
    .check <- function(msg, .options = NULL, method = nmf.getOption("default.algorithm"), ...) {
        be <- getDoBackend()
        expect_true(isNMFfit(res2 <- nmf(a, 2, seed = 123, nrun = 2, method = method,
                                         .opt = str_c("v3", .options), ...)), str_c(msg, " works"))
        # ns <- if (!is_NA(LOADED_NAMESPACES)) cl_loadedNamespaces(LOADED_NAMESPACES)
        expect_true(nmf.equal(resREF, res2), str_c(msg, ": result is identical to default"))
        expect_identical(consensus(resREF, no.attrib = TRUE), consensus(res2, 
                                                                        no.attrib = TRUE), str_c(msg, ": consensus matrice (no.attrib) is identical to default"))
        expect_identical(consensus(resREF), consensus(res2), str_c(msg, 
                                                                   ": consensus matrice is identical to default"))
        expect_true(identical(be, getDoBackend()), str_c(msg, ": backend is restored"))
        expect_error(nmf(a, 2, method = function(...) 1L, seed = 123, 
                         nrun = 2, .opt = str_c("v3", .options), ...), "should return an instance .* 'NMFfit'")
        expect_true(identical(be, getDoBackend()), str_c(msg, ": backend is restored after error"))
        # ns
    }
    #
    # create a method that just returns the result and set some value in the worker Global environment:
    # the value of variable is checked to see if the cluster was used.
    .fake_method <- function(y, x, abcd, res_object){
        assign("ABCD", abcd, envir = .GlobalEnv)
        res_object
        
    }
    .check_cluster_variable <- function(msg, cl, expected_value){
        expect_true(is.null(get0("ABCD", envir = .GlobalEnv)), paste0(msg, ": Local global environment is not affected"))
        n <- length(cl)
        val <- unlist(clusterApply(cl, seq(n), function(i) get0("ABCD", .GlobalEnv)), use.names = FALSE)
        expect_identical(val, rep(expected_value, n), paste0(msg, ": cluster value is set as expected"))
        
    }
    #
    library(parallel)
    cl <- makeCluster(2)
    on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
    
    .check("SEQ", .pbackend = "SEQ")
    if (NMF:::getAllCores() > 1){
        .check("P2", .options = "P2")
        .check("P2", .options = "P2", method = .fake_method, abcd = 1, res_object = resREF)
        .check_cluster_variable("P2", cl, NULL)
        
    }
    # PSOCK: independent processes
    .check(".pbackend=\"psock\"", .options = "P2", .pbackend = "PSOCK")
    .check(".pbackend=\"psock\"", .options = "P2", .pbackend = "PSOCK", method = .fake_method, abcd = 2, res_object = resREF)
    .check_cluster_variable(".pbackend=\"psock\"", cl, NULL)
    # Pre-defined cluster
    .check(".pbackend=cl + SNOW-like cluster", .pbackend = cl)
    .check(".pbackend=cl + SNOW-like cluster", .pbackend = cl, method = .fake_method, abcd = 3, res_object = resREF)
    .check_cluster_variable(".pbackend=cl + SNOW-like cluster", cl, 3)
    # 
    
    # Pre-registered cluster but run with P2 option: cluster is NOT used
    library(doParallel)
    registerDoParallel(cl)
    .check("doParallel registered cluster + P2 [should not use registered cluster]", .opt = "P2")
    .check("doParallel registered cluster + P2 [should not use registered cluster]", .opt = "P2", 
           method = .fake_method, abcd = 4, res_object = resREF)
    .check_cluster_variable("doParallel registered cluster + P2 [should not use registered cluster]", cl, 3)
    #
    
    # Pre-registered cluster with .pbackend = NULL: cluster used
    .check(".pbackend=NULL + doParallel registered cluster", .pbackend = NULL)
    .check(".pbackend=NULL + doParallel registered cluster", .pbackend = NULL, method = .fake_method, abcd = 5, res_object = resREF)
    .check_cluster_variable(".pbackend=NULL + doParallel registered cluster", cl, 5)
    
    if (!require.quiet("doMPI")) 
        skip("Package doMPI not available.")
    skip("doMPI checks are disabled due to issues in doMPI::closeCluster")
    cl_MPI <- startMPIcluster(2)
    on.exit(closeCluster(cl_MPI), add = TRUE)
    .check(".pbackend=cl_MPI + MPI cluster", .pbackend = cl_MPI)
    registerDoMPI(cl_MPI)
    .check(".pbackend=NULL + doMPI registered MPI cluster", .pbackend = NULL)
})

test_that("test.setupBackend", {
    on.exit(registerDoSEQ())
    expect_error(setupBackend(-1, "par", TRUE), info = "Invalid number of cores (optional)")
    expect_error(setupBackend(-1, "par", FALSE), info = "Invalid number of cores (required)")
    expect_error(setupBackend(10, "par", FALSE), info = "Required too many cores")
    expect_error(setupBackend(1, "toto", FALSE), info = "Required unknown backend")
})

test_that("test.shared_memory_doMC", {
    library(doParallel)
    registerDoParallel(2)
    check_shared_memory("doParallel - Multicore", libs = FALSE)
})

test_that("test.shared_memory_doMPI", {
    skip("NMF shared memory feature does not currently work with doMPI.")
    if (!require.quiet("doMPI"))
        skip("Package doMPI not available.")
    cl_MPI <- startMPIcluster(2)
    on.exit(closeCluster(cl_MPI), add = TRUE)
    registerDoMPI(cl_MPI)
    check_shared_memory("doMPI")
})

test_that("test.shared_memory_doParallel", {
    cl <- makeCluster(2, outfile = "wout.log")
    on.exit(stopCluster(cl), add = TRUE)
    registerDoParallel(cl)
    check_shared_memory("doParallel")
})

test_that("test.shared_memory_doSEQ", {
    registerDoSEQ()
    check_shared_memory("doSEQ", libs = FALSE, seq = TRUE)
})

