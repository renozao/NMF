#' Unit Testing script for NMF package: seeding methods.
#'
#' @author Renaud Gaujoux
#' @creation 17 Jul 2009
#' Converted from RUnit on 22 Feb 2020

.testData <- function(n=20, r=3, m=10, ...){
    syntheticNMF(n, r, m, ...)
}

nmf.options(maxIter = 10L)
on.exit(nmf.options(maxIter = NULL), add = TRUE)

#' Utility function for \code{test.seed}: performs a set of test on a seeded object
check.seed <- function(title, obj, V, r, seeding.meth, expect.class, exact.class=TRUE){
    
    expect_true( inherits(obj, 'NMFfit'), paste(title, ": result class inherits from 'NMFfit'") )
    expect_true( inherits(fit(obj), 'NMF'), paste(title, ": model class inherits from 'NMF'") )
    
    if( exact.class )
        expect_true( is(fit(obj), expect.class), paste(title, ": default class returned is '", expect.class, "'") )
    else
        expect_true( inherits(fit(obj), expect.class), paste(title, ": default class returned inherits from '", expect.class, "'") )
    
    expect_true( !is.empty.nmf(obj) , paste(title, ": Seeded object is not empty"))
    expect_equal( nbasis(obj), r , info = paste(title, ": Seeded object has correct rank"))
    expect_equal( nrow(obj), nrow(V) , info = paste(title, ": Seeded object has correct number of rows"), checkNames=FALSE)
    expect_equal( ncol(obj), ncol(V) , info = paste(title, ": Seeded object has correct number of columns"), checkNames=FALSE)
    expect_equal( seeding(obj), seeding.meth, info = "Seeding method's name is correctly set")
    
}

#' Utility check function: checks the range of value in a NMF object
check.range <- function(title, nmf.fit, max){
    obj <- fit(nmf.fit)
    expect_true( all(basis(obj) <= max & basis(obj) >= 0), paste(title, ': All entries of W are between 0 and', max) )
    expect_true( all(coef(obj) <= max & coef(obj) >= 0), paste(title, ': All entries of H are between 0 and', max) )
}

check.seed.change <- function(msg, expr, base){
    
    bs <- .Random.seed
    e <- parent.frame()
    eval(expr, env=e)
    
    if( base ) expect_true( any(bs != .Random.seed), paste(msg, ": .Random.seed IS changed")) 
    else expect_identical(bs, .Random.seed, paste(msg, ": .Random.seed is NOT changed"))
    
}


# Tests ----
test_that("test.nndsvd", {
    .seedTest <- 123456
    set.seed(.seedTest)
    r <- 3
    V <- .testData(r = r)
    check.seed.change("seeding with \"nndsvd\"", obj <- seed(V, 
        r, "nndsvd"), FALSE)
    check.seed("With matrix", obj, V, r, "nndsvd", "NMF")
    obj.bis <- seed(V, r, "nndsvd")
    expect_true(nmf.equal(obj.bis, obj, identical = FALSE), info = "Seeded NMF models are identical for every run")
})

test_that("test.none", {
    r <- 3
    V <- .testData(r = r)
    n <- nrow(V)
    m <- ncol(V)
    obj <- seed(V, r, "none")
    expect_true(is(obj, "NMFfit"), info = "Seeded object is an instance of class 'NMFfit'")
    expect_true(is(fit(obj), "NMFstd"), info = "Seeded model is an instance of class 'NMFstd'")
    expect_true(is.empty.nmf(obj), info = "Should not initialize the NMF object")
    expect_equal(r, nbasis(obj), info = "Seeded object have the correct rank")
    obj.init <- nmfModel(r)
    obj <- seed(V, obj.init, "none")
    expect_true(identical(fit(obj), obj.init), info = "Empty object: seeded object is identical to the initial one")
    obj.init <- nmfModel(r, model = "NMFstd", W = matrix(seq(n * 
        r), n, r), H = matrix(seq(r * m), r, m))
    obj <- seed(V, obj.init, "none")
    expect_true(identical(fit(obj), obj.init), info = "Dummy object: seeded object is identical to the initial one")
})

test_that("test.random", {
    .seedTest <- 123456
    V.na <- matrix(NA, 50, 20)
    r <- 3
    expect_warning(obj <- seed(V.na, r, "random"), "returning -Inf")
    check.range("NA matrix", obj, 1)
    max <- 0.05
    obj <- seed(matrix(max, 50, 20), r, "random")
    check.range(paste("Matrix", max), obj, max)
    set.seed(.seedTest)
    expect_warning(obj <- seed(V.na, r, "random"), "returning -Inf")
    check.seed("With matrix", obj, V.na, r, "random", "NMFstd")
    expect_warning(obj2 <- seed(V.na, r, "random"), "returning -Inf")
    expect_false(identical(obj2, obj), info = "Seeded objects are different if seed has not been fixed before seeding")
    set.seed(.seedTest)
    expect_warning(obj.bis <- seed(V.na, r, "random"), "returning -Inf")
    expect_true(nmf.equal(obj.bis, obj), info = "Seeded NMF models are identical if seed has been reset to the same value before seeding")
    set.seed(.seedTest)
    nmfOff <- nmfModel(r, model = "NMFOffset")
    expect_warning(obj <- seed(V.na, nmfOff, "random"), "returning -Inf")
    max <- 1
    check.seed("With NMFOffset object", obj, V.na, r, "random", 
        "NMFOffset")
    check.range(paste("Object NMFOffset", max), obj, max)
    expect_true(all(offset(fit(obj)) <= max & offset(fit(obj)) >= 
        0), info = paste("Object NMFOffset: All entries of Offset are between 0 and", 
        max))
    expect_warning(res <- seed(V.na, r, .seedTest), "returning -Inf")
    set.seed(.seedTest)
    expect_warning(obj <- seed(V.na, r, "random"), "returning -Inf")
    expect_true(nmf.equal(obj, res), info = "Seeded NMF models are identical when setting random generator seed and call method 'seed' with the same numerical seed (except for name of seeding method)")
})

test_that("test.restore", {
    skip("The option 'restore.seed' is deprecated. Related tests are now in test.seed.effect")
    r <- 3
    V <- .testData(r = r)
    os <- .Random.seed
    nmf(V, r)
    expect_false((all.equal(os, .Random.seed) == TRUE), info = "call with no seed: seed is correctly NOT restored")
    os <- .Random.seed
    nmf(V, r, seed = 1)
    expect_identical(.Random.seed, os, info = "Default behaviour is to restore the seed: seed is correctly restored")
    os <- .Random.seed
    nmf(V, r, .opt = "r", seed = 12)
    expect_identical(.Random.seed, os, info = "force seed restoration with 'r': seed correctly restored")
    os <- .Random.seed
    nmf(V, r, .opt = list(restore.seed = TRUE), seed = 123)
    expect_identical(.Random.seed, os, info = "force seed restoration with 'restore.seed=TRUE': seed correctly restored")
    os <- .Random.seed
    nmf(V, r, .opt = "-r", seed = 1234)
    expect_false((all.equal(os, .Random.seed) == TRUE), info = "Disable seed restoration with '-r': seed correctly NOT restored")
    os <- .Random.seed
    nmf(V, r, .opt = list(restore.seed = FALSE), seed = 12345)
    expect_false((all.equal(os, .Random.seed) == TRUE), info = "force seed restoration with 'restore.seed=FALSE': seed correctly NOT restored")
})

test_that("test.seed.effect", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    check.seed.change("After single run without seed", nmf(V, 
        r), TRUE)
    check.seed.change("After single run without seed (NO-REPRO has no effect)", 
        nmf(V, r), TRUE)
    check.seed.change("After single run with seed", nmf(V, r, 
        seed = 123), FALSE)
    check.seed.change("After multiple runs without seed (sapply)", 
        nmf(V, r, nrun = 3, .opt = "-p"), TRUE)
    check.seed.change("NO-REPRO: After multiple runs without seed (sapply)", 
        nmf(V, r, nrun = 3, .opt = "-pR"), TRUE)
    check.seed.change("After multiple runs without seed (foreach-MC)", 
        nmf(V, r, nrun = 3, .opt = "P", .pbackend = "par"), TRUE)
    check.seed.change("NO-REPRO: After multiple runs without seed (foreach-MC)", 
        nmf(V, r, nrun = 3, .opt = "P-R", .pbackend = "par"), 
        TRUE)
    check.seed.change("After multiple runs without seed (foreach-SEQ)", 
        nmf(V, r, nrun = 3, .opt = "P", .pbackend = "seq"), TRUE)
    check.seed.change("NO-REPRO: After multiple runs without seed (foreach-SEQ)", 
        nmf(V, r, nrun = 3, .opt = "P-R", .pbackend = "seq"), 
        TRUE)
    check.seed.change("After multiple runs with seed (sapply)", 
        nmf(V, r, nrun = 3, seed = 1234, .opt = "-p"), FALSE)
    check.seed.change("NO-REPRO: After multiple runs with seed (sapply)", 
        nmf(V, r, nrun = 3, seed = 1234, .opt = "-pR"), FALSE)
    check.seed.change("After multiple runs with seed (foreach-MC)", 
        nmf(V, r, nrun = 3, seed = 1234, .opt = "P", .pback = "par"), 
        FALSE)
    check.seed.change("NO-REPRO: After multiple runs with seed (foreach-MC)", 
        nmf(V, r, nrun = 3, seed = 1234, .opt = "P-R", .pback = "par"), 
        FALSE)
    check.seed.change("After multiple runs with seed (foreach-SEQ)", 
        nmf(V, r, nrun = 3, seed = 1234, .opt = "P", .pback = "seq"), 
        FALSE)
    check.seed.change("NO-REPRO: After multiple runs with seed (foreach-SEQ)", 
        nmf(V, r, nrun = 3, seed = 1234, .opt = "P-R", .pback = "seq"), 
        FALSE)
})

test_that("test.zzz.all", {
    set.seed(123)
    r <- 3
    V <- .testData(r = r)
    n <- nrow(V)
    m <- ncol(V)
    algorithms <- nmfAlgorithm()
    algorithms <- algorithms[!algorithms %in% c("ls-nmf", "pe-nmf", 
        "siNMF")]
    seed.methods <- nmfSeed()
    seed.methods <- seed.methods[which(seed.methods != "none")]
    test_algo <- function(name.algo, ..., target_rank = r) {
        sapply(seed.methods, function(name.seed, ...) {
            message("\n###########\n# ", name.algo, " + ", name.seed, 
                "\n#############")
            err <- try(obj <- nmf(..., method = name.algo, seed = name.seed))
            expect_true(!is(err, "try-error"), paste("Run OK - Algo:", 
                name.algo, "+ Seed:", name.seed, if (is(err, 
                  "try-error")) 
                  paste("[Error: ", err, "]")
                else NULL))
            check.seed(paste("Algo:", name.algo, "+ Seed:", name.seed), 
                obj, V, target_rank, name.seed, "NMF", exact.class = FALSE)
        }, ...)
    }
    sapply(algorithms, test_algo, V, r)
    test_algo("pe-nmf", V, r, alpha = 1, beta = 0.1)
    test_algo("ls-nmf", V, r, weight = rmatrix(V))
    g <- gl(2, m/2)
    test_algo("siNMF", V ~ g, r, target_rank = r + 2)
})

