#' Unit Testing script for NMF package: NMF interface for algorithms.
#'
#' @author Renaud Gaujoux
#' @creation 14 May 2009
# Converted from RUnit: 16 Feb 2020

library(rngtools)
checkIdenticalRNG <- function(x, y = getRNG(), msg){
    expect_true(checkRNG(x, y), msg)
    
}

.testData <- function(n=20, r=3, m=10, ...){
	syntheticNMF(n, r, m, ...)
}

#' Utility function for \code{test.seed}: performs a set of test on a seeded object
check.seed <- function(title, obj, V, r, expect.class){
		
	expect_true( isNMFfit(obj), paste(title, ": class returned is a valid NMF fit object") )
	expect_true( is(fit(obj), expect.class), paste(title, ": default class returned is ", expect.class) )
	expect_true( !is.empty.nmf(fit(obj)) , paste(title, ": Seeded object is not empty"))
	expect_equal( nbasis(obj), r , info = paste(title, ": Seeded object has correct rank"))
	expect_equivalent( nrow(obj), nrow(V) , info = paste(title, ": Seeded object has correct number of rows"), checkNames=FALSE)
	expect_equivalent( ncol(obj), ncol(V) , info = paste(title, ": Seeded object has correct number of columns"), checkNames=FALSE)
	
}

check.res <- function(title, obj, V, r, expect.class, algo=NULL, seed=NULL
					, rng=NULL, rngref=NULL){
	
	# check the same thing as in the seed
	check.seed(title, obj, V, r, expect.class)
	
	# check if some slots are correctly set
	if( is.null(algo) ) algo <- nmf.getOption('default.algorithm')
	expect_equal( algorithm(obj), algo, info = paste(title, ": Slot 'method' is correctly set"))
	
	if( is.null(seed) ) seed <- nmf.getOption('default.seed')
	expect_equal( seeding(obj), seed, info = paste(title, ": Slot 'seed' is correctly set"))
	
	# rng
	if( !is.null(rng) ){
				
		# check RNG
		expect_true( !rng.equal(rng), paste(title, ": The RNG after the fit is different from the one used to sed the computation"))
		if( nrun(obj) == 1 ){
			expect_true(rng.equal(obj, rng), paste(title, ": The fit's RNG seed is correctly set"))
		}else{			
			if( is.list(rng) ){				
				if( is(obj, 'NMFfitXn') )
					expect_true( all(mapply(rng.equal, obj, rng))
							, paste(title, ": The RNGs used in the multi-run computation are from the correct sequence"))
			}			
		}
		
		# check RNG_1
		rng1 <- if( is.list(rng) ) rng[[1]] else rng
		if( !is.null(rngref) )
			expect_true( !rng.equal(rng1, rngref), paste(title, ": The initial current RNG is different from the first RNG used in computation"))		
		expect_true(rng1.equal(obj, rng1), paste(title, ": The first fit's RNG seed is correctly set"))
		checkIdenticalRNG( getRNG1(obj), rng1, paste(title, ": The first RNG used in the computation is given by getRNG1"))
	}
	# ref rng 
	if( !is.null(rngref) ){
		expect_true( rng.equal(rngref), paste(title, ": The current RNG was not affected by the computation"))
		if( is.null(rng) )
			expect_true( rng.equal(obj, rngref), paste(title, ": The fit's RNG is the same as the reference RNG"))
	}
}


test_that("test.compare", {
   compare <- NMF::compare
    r <- 3
    V <- .testData(r = r)
    m <- ncol(V)
    res <- list()
    res$brunet <- nmf(V, r, "brunet")
    res$ns <- nmf(V, r, "ns")
    res$offset <- nmf(V, r, "off")
    res$snmfr <- nmf(V, r, "snmf/r")
    res$snmfl <- nmf(V, r, "snmf/l")
    classes <- as.factor(sample(seq(r), m, replace = TRUE))
    expect_true(is.data.frame(compare(res, class = classes)), 
        info = "Result of method 'compare' (list) is a data.frame")
    expect_true(is.data.frame(compare(res, class = classes, sort.by = "method")), 
        info = "Result of method 'compare' (list) is a data.frame")
    res$brunet <- nmf(V, r, "brunet", nrun = 10)
    expect_true(is.data.frame(compare(res, class = classes)), 
        info = "Result of method 'compare' (list with multiple runs) is a data.frame")
    res <- nmf(V, r, list("brunet", "lee", "offset"))
    expect_true(is.data.frame(compare(res, class = classes)), 
        info = "Result of method 'compare' (list with multiple methods) is a data.frame")
    res <- nmf(V, r, list("brunet", "lee", "ns"), nrun = 3)
    expect_true(is.data.frame(compare(res, class = classes)), 
        info = "Result of method 'compare' (list with multiple runs + multiple methods) is a data.frame")
})

test_that("test.nmf.callback", {
    r <- 3
    V <- .testData(r = r)
    cb <- function(object, i) {
        1
    }
    res <- nmf(V, r, nrun = 3, .callback = cb, .opt = "P")
    expect_equal(rep(1, 3), res$.callback, info = "Result of callback is in: res$.callback (PAR)")
    res <- nmf(V, r, nrun = 3, .callback = cb, .opt = "-p")
    expect_equal(rep(1, 3), res$.callback, info = "Result of callback is in: res$.callback (SEQ)")
    cb <- algorithm
    res <- nmf(V, r, nrun = 3, .callback = cb, .opt = "P")
    expect_equal(rep("brunet", 3), res$.callback, info = "Result of callback can use the result object of each run (PAR)")
    res <- nmf(V, r, nrun = 3, .callback = cb, .opt = "-p")
    expect_equal(rep("brunet", 3), res$.callback, info = "Result of callback can use the result object of each run (SEQ)")
    cb <- function(object, i) {
        stop()
    }
    expect_warning(res <- nmf(V, r, "br", nrun = 3, .callback = cb, 
        .opt = "Pk"), "discarding argument .*\\.callback")
    expect_true(is(res, "NMFfitXn"), info = "Callback function is not used with option 'keep.all=TRUE' (PAR)")
    expect_warning(res <- nmf(V, r, "br", nrun = 3, .callback = cb, 
        .opt = "k-p"), "discarding argument .*\\.callback")
    expect_true(is(res, "NMFfitXn"), info = "Callback function is not used with option 'keep.all=TRUE' (SEQ)")
    cb <- function(object, i) {
        stop("BIG ERROR")
    }
    expect_warning(res <- nmf(V, r, "br", nrun = 3, .callback = cb, .opt = "P"),
                   "All NMF fits .* callback .* threw an error.* run #1: BIG ERROR")
    expect_true(isNMFfit(res), info = "Error in callback function does not stop the copmutation (PAR)")
    expect_true(is.list(res$.callback), info = "res$.callback is a list when there is an error (PAR)")
    expect_equal(rep(TRUE, 3), sapply(res$.callback, function(x) is(x, 
        "error")), info = "Error in callback function returns errors in res$.callback (PAR)")
    
    expect_warning(res <- nmf(V, r, "br", nrun = 3, .callback = cb, .opt = "-p"),
                   "All NMF fits .* callback .* threw an error.* run #1: BIG ERROR")
    expect_true(isNMFfit(res), info = "Error in callback function does not stop the copmutation (SEQ)")
    expect_true(is.list(res$.callback), info = "res$.callback is a list when there is an error (SEQ)")
    expect_equal(rep(TRUE, 3), sapply(res$.callback, function(x) is(x, 
        "error")), info = "Error in callback function returns errors in res$.callback (SEQ)")

    res <- nmf(V, r, "br", nrun = 3, .callback = summary, .opt = "P")
    expect_true(is.matrix(res$.callback), info = "res$.callback is a list when there is NO error (PAR)")
    res <- nmf(V, r, "br", nrun = 3, .callback = summary, .opt = "P-S")
    expect_true(is.list(res$.callback), info = "res$.callback is a list when there is NO error (PAR) and simplifyCB=FALSE")
    res <- nmf(V, r, "br", nrun = 3, .callback = summary, .opt = "-p")
    expect_true(is.matrix(res$.callback), info = "res$.callback is a list when there is NO error (SEQ)")
    res <- nmf(V, r, "br", nrun = 3, .callback = summary, .opt = "-pS")
    expect_true(is.list(res$.callback), info = "res$.callback is a list when there is NO error (SEQ) and simplifyCB=FALSE")
    cb <- function(object, i) {
        if (i == 1) 
            stop("BIG ERROR ", i)
        summary(object)
    }
    expect_warning(res <- nmf(V, r, "br", nrun = 3, .callback = cb, .opt = "P"),
                   "All NMF fits .* callback .* threw an error.* run #1: BIG ERROR")
    expect_true(is.list(res$.callback), info = "res$.callback is a list when there is at least one error (PAR)")
    expect_equal(c(TRUE, FALSE, FALSE), sapply(res$.callback, 
        function(x) is(x, "error")), info = "Error in callback function returns errors mixed with values in res$.callback (PAR)")
    expect_warning(res <- nmf(V, r, "br", nrun = 3, .callback = cb, .opt = "-p"),
                   "All NMF fits .* callback .* threw an error.* run #1: BIG ERROR")
    expect_true(is.list(res$.callback), info = "res$.callback is a list when there is at least one error (SEQ)")
    expect_equal(c(TRUE, FALSE, FALSE), sapply(res$.callback, 
        function(x) is(x, "error")), info = "Error in callback function returns errors mixed with values in res$.callback (SEQ)")
})

test_that("test.nmf.custom", {
    r <- 3
    V <- .testData(r = r)
    my.algo <- function(x, seed, theta = 0) {
        seed$extra.param <- theta
        seed
    }
    res <- nmf(V, r, my.algo, name = "dummy")
    check.res("No argument (only name)", res, V, r, "NMFstd", 
        "dummy")
    expect_equal(0, res$extra.param, info = "No argument: NO argument is provided to the algorithm")
    res <- nmf(V, r, my.algo, name = "dummy", model = list())
    check.res("Argument 'model' an empty list", res, V, r, "NMFstd", 
        "dummy")
    expect_equal(0, res$extra.param, info = "No argument: NO argument is provided to the algorithm")
    expect_error(nmf(V, r, my.algo, name = "dummy", model = list(theta = 1)), 
        info = "Error when passing non model parameter in argument 'model', with same name as an algo parameter")
    res <- nmf(V, r, my.algo, name = "dummy", theta = 10)
    check.res("No argument (only name)", res, V, r, "NMFstd", 
        "dummy")
    expect_equal(10, res$extra.param, info = "NMFstd + Extra argument in '...': extra argument IS provided to the algorithm")
    expect_error(nmf(V, r, my.algo, name = "dummy", toto = 1), 
        info = "Error if NMFstd + unused parameter in argument '...'")
    res <- nmf(V, r, my.algo, model = "NMFns")
    default.theta <- fit(res)@theta
    custom.theta <- 1
    expect_equal(0, res$extra.param, info = "NMFns + no args: extra argument is NOT provided to the algorithm")
    res <- nmf(V, r, my.algo, name = "dummy", model = list("NMFns", 
        theta = custom.theta))
    check.res("NMFns + With extra argument in argument 'model'", 
        res, V, r, "NMFns", "dummy")
    expect_equal(0, res$extra.param, info = "NMFns + Argument in 'model': extra argument is NOT provided to the algorithm")
    expect_equal(custom.theta, fit(res)@theta, info = "NMFns + Argument in 'model': extra argument IS used in model")
    res <- nmf(V, r, my.algo, name = "dummy", model = "NMFns", 
        theta = custom.theta)
    check.res("NMFns + With extra argument in argument '...' and 'model' a model name", 
        res, V, r, "NMFns", "dummy")
    expect_equal(0, res$extra.param, info = "NMFns + Argument in '...', 'model' a model name: extra argument is NOT provided to the algorithm")
    expect_equal(custom.theta, fit(res)@theta, info = "NMFns + Argument in '...' and 'model' a model name: extra argument used in model")
    res <- nmf(V, r, my.algo, name = "dummy", model = list("NMFns"), 
        theta = 1)
    check.res("NMFns + With extra argument in argument '...' and 'model' a list with model name", 
        res, V, r, "NMFns", "dummy")
    expect_equal(1, res$extra.param, info = "NMFns + Argument in '...' and 'model' is a list with model name: extra argument IS provided to the algorithm")
    expect_equal(default.theta, fit(res)@theta, info = "NMFns + Argument in '...' and 'model' is list with model name: extra argument NOT used in model")
    res <- nmf(V, r, my.algo, name = "dummy", model = list("NMFns", 
        theta = custom.theta), theta = 1)
    check.res("NMFns + With extra argument in argument 'model'", 
        res, V, r, "NMFns", "dummy")
    expect_equal(1, res$extra.param, info = "NMFns + Different values in argument in '...' and 'model': correct extra argument IS provided to the algorithm")
    expect_equal(custom.theta, fit(res)@theta, info = "NMFns + Different values in '...' and 'model': correct extra argument IS used in model")
    V.neg <- V
    V.neg[1, 1] <- -1
    expect_error(nmf(V.neg, r, my.algo, name = "dummy"), info = "Throw an error if some input entries are negative and algoritham is declared NOT MIXED")
    res <- nmf(V.neg, r, my.algo, name = "dummy", mixed = TRUE)
    check.res("Call with dummy MIXED algorithm on input with negative entries", 
        res, V.neg, r, "NMFstd", "dummy")
    my.algo <- function(target, start, param1, param2) {
        basis(start)[1, ] <- -1
        start
    }
    res <- nmf(V, r, my.algo, name = "dummy")
    check.res("Call with dummy algorithm and MIXED output", res, 
        V, r, "NMFstd", "dummy")
})

test_that("test.nmf.default", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    check.res("Call with rank (only)", nmf(V, r), V, r, "NMFstd")
})

test_that("test.nmf.dots", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    expect_error(nmf(V, r, toto = 5), info = "Error if there is an unused parameter in '...'")
    expect_error(nmf(V, r, "nsNMF", model = list(), theta = 0.6), 
        info = "Forcing argument to go to algo: error if there is an unused parameter in '...'")
    res <- nmf(V, r, "nsNMF", theta = 0.6)
    check.res("Call with argument '...' to specify extra initialization parameter", 
        res, V, r, "NMFns", "nsNMF")
    expect_equal(0.6, fit(res)@theta, info = "Call algo:nsNMF with theta in argument '...': argument correctly passed to model")
})

test_that("test.nmf.method", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    expect_error(nmf(V, r, "zzz"), info = "Throw an error when: inexistent algorithm name")
    expect_error(nmf(V, r, toto = 3), info = "Throw an error when: unused argument is passed to algorithm")
    check.res("Call with algorithm name", nmf(V, r, "nsNMF"), 
        V, r, "NMFns", "nsNMF")
    check.res("Call with algorithm name and partial match", nmf(V, 
        r, "ns"), V, r, "NMFns", "nsNMF")
    old.rseed <- getRNG()
    res <- nmf(V, r, list("ns", "br", "lee"))
    expect_identical(c("nsNMF", "brunet", "lee"), names(res), 
        info = "Argument list(): names are set correctly to the complete method names")
    expect_true(all(sapply(res, function(x) identical(getRNG(x), 
        getRNG(res[[1]])))), info = "Initial RNG settings are the same for each method")
    expect_false(identical(old.rseed, getRNG()), info = "RNG setting is different after the run")
    new.rseed <- getRNG()
    setRNG(old.rseed)
    nmf(V, r, "lee")
    expect_identical(getRNG(), new.rseed, info = "RNG setting after the run is the same as if one has run only the last method")
    res <- nmf(V, r, list("ns", "br", "lee"), nrun = 3)
    expect_identical(c("nsNMF", "brunet", "lee"), names(res), 
        info = "Argument list() + multiple run: names are set correctly to the complete method names")
    expect_true(all(sapply(res, function(x) identical(getRNG1(x), 
        getRNG1(res[[1]])))), info = "Argument list() + multiple runs: Initial RNG settings are the same for each method")
    ml <- list("ns", "brunet", "lee")
    expect_error(nmf(V, r, ml, .parameters = 2:3), info = "Error if argument .parameters not a list")
    expect_error(nmf(V, r, ml, .parameters = list(list(copy = TRUE))), 
        info = "Error if argument .parameters has no names")
    expect_error(nmf(V, r, ml, .parameters = list(br = list(), 
        list(copy = TRUE))), info = "Error if argument .parameters has missing names")
    expect_error(nmf(V, r, ml, .parameters = list(br = list(), 
        brun = list(copy = TRUE))), info = "Error if argument .parameters has multiple matching names")
    expect_warning(nmf(V, r, ml, .parameters = list(br = list(aaa = 1))), 
        info = "Error if unused argument in selected method-specific parameters")
    expect_warning(nmf(V, r, ml, .parameters = list(br = list(), 
        toto = list())), info = "Warning if unused elements in .parameters")
    expect_true(all(isNMFfit(res <- nmf(V, r, ml, seed = 123, 
        .parameters = list(br = list(maxIter = 10), ns = list(maxIter = 2))))), 
        info = "List of methods working if called with correct .parameters")
    expect_identical(2L, niter(res$nsNMF))
    expect_identical(10L, niter(res$brunet))
    res_lee <- nmf(V, r, "lee", seed = 123)
    expect_true(niter(res$lee) > 10L, info = "Method without method-specific parameter specification correctly runs without them: niter > 10L")
    expect_identical(niter(res$lee), niter(res$lee), info = "Method without method-specific parameter specification correctly runs without them: niter equal")
    expect_true(nmf.equal(res$lee, res_lee), info = "Method without method-specific parameter specification correctly runs without them: identical result")
})

test_that("test.nmf.model", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    check.res("Call with empty argument 'model'", nmf(V, r, model = list()), 
        V, r, "NMFstd")
    expect_error(nmf(V, r, model = NA), info = "Error if there argument 'model' is of a bad type: NA")
    expect_error(nmf(V, r, model = "toto"), info = "Error if there argument 'model' is of a bad type: character")
    expect_error(nmf(V, r, model = 12), info = "Error if there argument 'model' is of a bad type: numeric")
    expect_error(nmf(V, r, model = list("toto")), info = "Error if there is a not named element in argument 'model'")
    expect_error(nmf(V, r, model = list(toto = 5)), info = "Error if there is a bad slot name in argument 'model'")
    res <- nmf(V, r, "nsNMF", model = list(theta = 0.6))
    check.res("Call with argument 'model' to specify extra initialization parameter", 
        res, V, r, "NMFns", "nsNMF")
    expect_equal(0.6, fit(res)@theta, info = "Call algo:nsNMF with theta in argument 'model': argument correctly passed to model")
})

test_that("test.nmf.multirank", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    old.rseed <- getRNG()
    ranks <- 2:4
    old.rseed <- getRNG()
    res <- nmf(V, ranks, nrun = 1)
    expect_true(is(res, "NMF.rank"), info = "Result is a NMF.rank object")
    expect_identical(c("measures", "consensus", "fit"), names(res), 
        info = "result names are corrects")
    expect_identical(as.character(ranks), names(res$fit), info = "names of fits are the ranks as character strings")
    fits <- res$fit
    expect_true(all(sapply(fits, function(x) identical(getRNG(x), 
        getRNG(fits[[1]])))), info = "Initial RNG settings are the same for each rank")
    expect_false(identical(old.rseed, getRNG()), info = "RNG setting is different after the run")
    new.rseed <- getRNG()
    setRNG(old.rseed)
    nmf(V, tail(ranks, 1))
    expect_identical(getRNG(), new.rseed, info = "RNG setting after the run is the same as if one has run only the last rank")
    res <- nmf(V, ranks, nrun = 3)
    expect_true(is(res, "NMF.rank"), info = "multiple runs: Result is a NMF.rank object")
    expect_identical(c("measures", "consensus", "fit"), names(res), 
        info = "multiple runs: result names are corrects")
    expect_identical(as.character(ranks), names(res$fit), info = "multiple runs: names of fits are the ranks as character strings")
    fits <- res$fit
    expect_true(all(sapply(fits, function(x) identical(getRNG1(x), 
        getRNG1(fits[[1]])))), info = "multiple runs: Initial RNG settings are the same for each rank")
})

test_that("test.nmf.options", {
    x <- rmatrix(20, 10)
    .check <- function(msg, it, ...) {
        .msg <- function(...) paste(msg, ":", ...)
        res <- nmf(x, 2, ...)
        t <- residuals(res, track = TRUE)
        expect_true(!is.null(names(t)), .msg("Track has names"))
        expect_true(0 == names(t)[1], .msg("First value in track is for iteration 0"))
        t <- t[-1]
        lags <- head(diff(as.numeric(names(t))), length(t) - 
            2)
        expect_identical(lags, rep(it, length(t) - 2), .msg("Track interval is correct"))
    }
    .test <- function(msg, ...) {
        .check(paste(msg, "- single run"), ..., nrun = 1)
        .check(paste(msg, "- multiple runs"), ..., nrun = 3)
    }
    .test("Default call -> use option", nmf.getOption("track.interval"), 
        .options = "t")
    .test("Specified in .options=\"t5\" -> use value from \"t\"", 
        5, .options = "t5")
    nmf.options(track.interval = 7)
    .test("Default call after changing option -> use new option", 
        7, .options = "t")
})

test_that("test.nmf.parameters", {
    my.algo <- function(target, start, param1, param2) {
        start
    }
    r <- 3
    V <- .testData(r = r)
    check.res("Call with custom algo, model NMF", nmf(V, r, my.algo, 
        name = "dummy", model = "NMFstd"), V, r, "NMFstd", "dummy", 
        "random")
    check.res("Call with custom algo, model NMFns", nmf(V, r, 
        my.algo, name = "dummy", model = "NMFns"), V, r, "NMFns", 
        "dummy", "random")
    res <- nmf(V, r, my.algo, name = "dummy", model = list("NMFns", 
        theta = 0.3))
    check.res("Call with custom algo, model NMFns, theta in argument model: TEST RES", 
        res, V, r, "NMFns", "dummy", "random")
    expect_equal(0.3, fit(res)@theta, info = "Call with custom algo, model NMFns, theta in argument model: argument correctly passed to model")
    res <- nmf(V, r, my.algo, name = "dummy", model = "NMFns", 
        theta = 0.3)
    check.res("Call with custom algo, model NMFns, theta in extra argument: TEST RES", 
        res, V, r, "NMFns", "dummy", "random")
    expect_equal(0.3, fit(res)@theta, info = "Call with custom algo, model NMFns, theta in extra argument: argument correctly passed to model")
    res <- nmf(V, r, my.algo, name = "dummy", model = "NMFstd", 
        param1 = 0.6)
    check.res("Call with custom algo, model NMFns, plus an extra argument: TEST RES", 
        res, V, r, "NMFstd", "dummy", "random")
    expect_equal(list(param1 = 0.6), res@parameters, info = "Call with custom algo, model NMFns, plus an extra argument: argument is passed correctly to algorithm")
    my.algo2 <- function(target, start, theta) {
        start
    }
    res <- nmf(V, r, my.algo2, name = "dummy", model = list("NMFns", 
        theta = 0.3), theta = 0.6)
    check.res("Call with custom algo, model NMFns, theta in argument model AND extra argument: TEST RES", 
        res, V, r, "NMFns", "dummy", "random")
    expect_equal(0.3, fit(res)@theta, info = "Call with custom algo, model NMFns, theta in argument model AND extra argument: argument model passed to model")
    expect_equal(list(theta = 0.6), res@parameters, info = "Call with custom algo, model NMFns, theta in argument model AND extra argument: extra argument passed to algorithm")
})

test_that("test.nmf.seed.argument", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    expect_error(nmf(V, r, seed = "zzz"), info = "Throw an error when: inexistent seeding method name")
    expect_error(nmf(V, r, seed = matrix(1, 2, 2)), info = "Throw an error when: invalid object as seeding method")
    expect_error(nmf(V, r, seed = list("zzz")), info = "Throw an error when: inexistent seeding method name (passed as list)")
    expect_error(nmf(V, r, seed = list(method = "zzz")), info = "Throw an error when: inexistent seeding method name (passed as named list)")
    expect_error(nmf(V, r, seed = list(toto = 1, method = "random")), 
        info = "Throw an error when: unused argument is passed to seeding method")
    expect_error(nmf(V, r, seed = numeric()), info = "Throw an error when: seed argument is an empty numeric")
    expect_error(nmf(V, r, seed = c(1, 2)), info = "Throw an error when: seed argument is a numeric of invalid length (2)")
    expect_error(nmf(V, r, seed = rep(5, 7)), info = "Throw an error when: seed argument is an invalid numeric value for .Random.seed (7)")
    set.seed(123)
    rngRef <- getRNG()
    check.res("Call with only a NON random seeding method name", 
        nmf(V, r, seed = "nndsvd"), V, r, "NMFstd", seed = "nndsvd", 
        rngref = rngRef)
    check.res("Call with only a NON random seeding method name partially match", 
        nmf(V, r, seed = "nnd"), V, r, "NMFstd", seed = "nndsvd")
    set.seed(1234)
    rngRef <- getRNG()
    check.res("Call with \"random\" seeding method name", res <- nmf(V, 
        r, seed = "random"), V, r, "NMFstd", seed = "random", 
        rng = rngRef)
    msg <- function(...) paste("Call with only a numerical seed:", 
        ...)
    rngRef <- getRNG()
    s <- nextRNG(123456)
    check.res(msg(), res <- nmf(V, r, seed = 123456), V, r, "NMFstd", 
        seed = "random", rng = s, rngref = rngRef)
    msg <- function(...) paste("Call with a 6-length numerical seed:", 
        ...)
    runif(10)
    rngRef <- getRNG()
    nseed <- c(1, 2, 3, 4, 5, 6)
    s <- RNGseq(1, nseed)
    check.res(msg(), res <- nmf(V, r, seed = nseed), V, r, "NMFstd", 
        seed = "random", rng = s, rngref = rngRef)
    msg <- function(...) paste("Multirun - parallel + numeric seed (keep all):", 
        ...)
    runif(10)
    rngRef <- getRNG()
    sRNG <- RNGseq(3, seed = 5698)
    check.res(msg(), res <- nmf(V, r, nrun = 3, seed = 5698, 
        .opt = "kP"), V, r, "NMFstd", seed = "random", rng = sRNG, 
        rngref = rngRef)
    msg <- function(...) paste("Multirun - parallel + list of seeds (keep all):", 
        ...)
    runif(10)
    rngRef <- getRNG()
    check.res(msg(), res2 <- nmf(V, r, nrun = 3, seed = sRNG, 
        .opt = "kP"), V, r, "NMFstd", seed = "random", rng = sRNG, 
        rngref = rngRef)
    checkIdenticalRNG(res2, res, msg("The best fit's RNG is the same as when seeding with corresponding single numeric seed"))
    msg <- function(...) paste("Multirun - parallel + numeric seed (keep best):", 
        ...)
    runif(10)
    res2 <- nmf(V, r, nrun = 3, seed = 5698, .opt = "P")
    checkIdenticalRNG(res2, res, msg("The best fit's RNG is the same as when keeping all the fits"))
    checkIdenticalRNG(getRNG1(res2), sRNG[[1]], msg("The first RNG used in the computation of the NMFfitX1 object is given by getRNG1"))
    expect_true(rng1.equal(res2, sRNG[[1]]), info = msg("The first RNG used in the computation is correct"))
    msg <- function(...) paste("Multirun - parallel + single 7-length numeric seed (keep best):", 
        ...)
    runif(10)
    res2 <- nmf(V, r, nrun = 3, seed = sRNG[[1]], .opt = "P")
    checkIdenticalRNG(res2, res, msg("The best fit's RNG is the same as when keeping all the fits"))
    checkIdenticalRNG(getRNG1(res2), sRNG[[1]], msg("The first RNG used in the computation of the NMFfitX1 object is given by getRNG1"))
    expect_true(rng1.equal(res2, sRNG[[1]]), info = msg("The first RNG used in the computation is correct"))
    msg <- function(...) paste("Multirun - parallel + rstream seed (keep best):", 
        ...)
    runif(10)
    res2 <- nmf(V, r, nrun = 3, seed = sRNG[[1]], .opt = "P")
    checkIdenticalRNG(res2, res, msg("The best fit's RNG is the same as when keeping all the fits"))
    checkIdenticalRNG(getRNG1(res2), sRNG[[1]], msg("The first RNG used in the computation of the NMFfitX1 object is given by getRNG1"))
    obj.s <- rnmf(r, V)
    rngRef <- getRNG()
    res <- nmf(V, obj.s)
    check.res("Call with rank = <NMF object>", res, V, r, "NMFstd", 
        "brunet", "NMF", rngref = rngRef)
    expect_true(nmf.equal(res, nmf(V, obj.s)), info = "Run with rank=<NMF object> is deterministic")
    res.s <- nmf(V, seed = obj.s)
    check.res("Call with seed = <NMF object>", res, V, r, "NMFstd", 
        "brunet", "NMF", rngref = rngRef)
    expect_true(nmf.equal(res, res.s), info = "Run with rank=<NMF object> returns identical result as with seed=<NMF object>")
    check.res("Call with only a seeding method name and some extra parameters (element method first and not named)", 
        nmf(V, r, seed = list("nndsvd", densify = "average")), 
        V, r, "NMFstd", seed = "nndsvd")
    check.res("Call with only a seeding method name and some extra parameters (element method second and named)", 
        nmf(V, r, seed = list(densify = "average", method = "nndsvd")), 
        V, r, "NMFstd", seed = "nndsvd")
    check.res("Call with both algorithm and seeding method", 
        nmf(V, r, "lee", seed = "nndsvd"), V, r, "NMFstd", "lee", 
        "nndsvd")
})

test_that("test.nmf.seed.equivalent", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    set.seed(1234)
    ss_SEQ <- nmf(V, r, nrun = 3, .opt = "k-p")
    runif(10)
    SEQ <- nmf(V, r, nrun = 3, seed = 1234, .opt = "k-p")
    expect_true(nmf.equal(ss_SEQ, SEQ, all = TRUE), info = "Multiple run using sapply with a numeric seed is equivalent to set.seed + run (sapply)")
    set.seed(1234)
    ss_PAR <- nmf(V, r, nrun = 3, .opt = "kP2")
    runif(10)
    PAR <- nmf(V, r, nrun = 3, seed = 1234, .opt = "kP2")
    expect_true(nmf.equal(ss_SEQ, PAR, all = TRUE), info = "Multiple run using foreach with a numeric seed is equivalent to set.seed + run (sapply)")
    expect_true(nmf.equal(ss_PAR, PAR, all = TRUE), info = "Multiple run using foreach with a numeric seed is equivalent to set.seed + run (foreach)")
    PAR_SEQ <- nmf(V, r, nrun = 3, seed = 1234, .opt = "kP", 
        .pbackend = "seq")
    expect_true(nmf.equal(PAR_SEQ, PAR, all = TRUE), info = "Multiple run using foreach with a numeric seed is equivalent to foreach sequential with numeric seed")
    set.seed(1234)
    ss_PAR_SEQ <- nmf(V, r, nrun = 3, .opt = "kP", .pbackend = "seq")
    expect_true(nmf.equal(ss_PAR_SEQ, PAR, all = TRUE), info = "Multiple run using foreach with a numeric seed is equivalent to set.seed + foreach sequential")
    set.seed(1234)
    ss_SEQ_noR <- nmf(V, r, nrun = 3, .opt = "k-pR")
    runif(10)
    SEQ_noR <- nmf(V, r, nrun = 3, seed = 1234, .opt = "k-pR")
    expect_true(nmf.equal(ss_SEQ_noR, SEQ_noR, all = TRUE), info = "Multiple run using sapply with a numeric seed WITHOUT option R is equivalent to set.seed + run (sapply) WITHOUT option R")
    expect_false(nmf.equal(SEQ_noR, SEQ, all = TRUE), info = "Multiple run using sapply with a numeric seed WITHOUT option R is NOT equivalent to set.seed + run (sapply)")
    set.seed(1234)
    ss_SEPA <- replicate(3, nmf(V, r))
    expect_false(nmf.equal(ss_SEPA, ss_SEQ, all = TRUE), info = "Fits of multiple runs with sapply and a seed are NOT the same as the ones obtained from set.seed + separate fits")
    expect_true(nmf.equal(ss_SEPA, ss_SEQ_noR, all = TRUE), info = "Fits of multiple runs WITHOUT option R with sapply and a seed are the same as the ones obtained from set.seed + separate fits")
    expect_false(nmf.equal(ss_SEPA, PAR, all = TRUE), info = "Fits of multiple runs with foreach and a seed are NOT the same as the ones obtained from set.seed + separate fits")
    PAR_noR <- nmf(V, r, nrun = 3, seed = 1234, .opt = "kP2-R")
    expect_true(nmf.equal(PAR_noR, PAR, all = TRUE), info = "Option -R has no effect on true parallel computations")
})

test_that("test.nmf.seed.fault", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    os <- .Random.seed
    try(res <- nmf(V, r, method = function(...) {
    }))
    expect_identical(.Random.seed, os, info = ".Random.seed is NOT changed after error in single run without seed")
    os <- .Random.seed
    try(res <- nmf(V, r, nrun = 3, method = function(...) {
    }, .opt = "-p"))
    expect_identical(.Random.seed, os, info = ".Random.seed is NOT changed after error in multiple runs without seed (sapply)")
    os <- .Random.seed
    try(res <- nmf(V, r, nrun = 3, method = function(...) {
    }, .opt = "P2"))
    expect_identical(.Random.seed, os, info = ".Random.seed is NOT changed after error in multiple runs without seed (foreach-MC)")
    os <- .Random.seed
    try(res <- nmf(V, r, nrun = 3, method = function(...) {
    }, .opt = "P1"))
    expect_identical(.Random.seed, os, info = ".Random.seed is NOT changed after error in multiple runs without seed (foreach-SEQ)")
})

test_that("test.nmf.seed.repro", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    res <- replicate(2, nmf(V, r, seed = 123))
    expect_true(nmf.equal(res[[1]], res[[2]]), info = "Results are reproducible: single run")
    res <- replicate(2, nmf(V, r))
    expect_false(nmf.equal(res[[1]], res[[2]]), info = "Results are NOT the same if not seeded: single run")
    res <- replicate(2, nmf(V, r, nrun = 3, seed = 123, .opt = "kP"), 
        simplify = FALSE)
    expect_true(nmf.equal(res[[1]], res[[2]], all = TRUE), info = "Results are reproducible: multiple run - Parallel")
    res <- replicate(2, nmf(V, r, nrun = 3, .opt = "kP"), simplify = FALSE)
    expect_false(nmf.equal(res[[1]], res[[2]]), info = "Results are NOT the same if not seeded: multiple run - Parallel")
    res <- replicate(2, nmf(V, r, nrun = 3, seed = 123, .opt = "k-p"), 
        simplify = FALSE)
    expect_true(nmf.equal(res[[1]], res[[2]], all = TRUE), info = "Results are reproducible: multiple run - sapply")
    res <- replicate(2, nmf(V, r, nrun = 3, .opt = "k-p"), simplify = FALSE)
    expect_false(nmf.equal(res[[1]], res[[2]]), info = "Results are NOT the same if not seeded: multiple run - sapply")
    res <- list(nmf(V, r, nrun = 3, seed = 123, .opt = "kP", 
        .pbackend = "seq"), nmf(V, r, nrun = 3, seed = 123, .opt = "kP", 
        .pbackend = "mc"))
    expect_true(nmf.equal(res[[1]], res[[2]], all = TRUE), info = "Identical results from seeded foreach MC and SEQ")
    res <- list(nmf(V, r, nrun = 3, .opt = "kP", .pbackend = "seq"), 
        nmf(V, r, nrun = 3, .opt = "kP", .pbackend = "mc"))
    expect_false(nmf.equal(res[[1]], res[[2]], all = TRUE), info = "NON-identical results from non-seeded foreach MC and SEQ")
})

test_that("test.nmf.stop", {
    r <- 3
    V <- .testData(r = r)
    if (is.null(maxIter <- nmf.getOption("maxIter"))) 
        maxIter <- 2000L
    maxIter <- maxIter + 100L
    expect_identical(37L, niter(nmf(V, r, .stop = 37L)), info = "Integer stopping criterium: fixed number of iterations")
    expect_identical(maxIter, niter(nmf(V, r, .stop = maxIter)), 
        info = "Integer stopping criterium greater than default maxIter: fixed number of iterations is honoured")
    expect_identical(67L, niter(nmf(V, r, .stop = 200L, maxIter = 67L)), 
        info = "Integer stopping criterium greater than provided maxIter: maxIter is honoured")
    expect_true(niter(nmf(V, r, .stop = 37)) != 37, info = "Numeric stopping criterium: stationarity threshold")
    expect_true(niter(nmf(V, r, .stop = "nmf.stop.stationary")) != 
        37, info = "stopping criterium 'stationary' as full name is passed correctly")
    expect_true(niter(nmf(V, r, .stop = "stationary")) != 37, 
        info = "stopping criterium 'stationary' as short name is passed correctly")
})

test_that("test.nmfModel.formula", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    w <- rmatrix(nrow(V), r)
    h <- rmatrix(r, ncol(V))
    cx <- runif(ncol(V))
    bx <- runif(nrow(V))
    .check <- function(res, dims, msg, cterm = NULL, bterm = NULL) {
        .msg <- function(...) paste0(msg, ": ", ...)
        expect_true(isNMFfit(res), .msg("Result is an NMFfit object"))
        expect_equal(dim(res), dims, info = .msg("Dimensions [", str_dim(res), 
            "] are as expected [", str_dim(dims = dims), "]"))
        if (!is.null(cterm)) {
            if (is.null(dim(cterm))) 
                cterm <- matrix(cterm, 1L)
            else if (is.data.frame(cterm)) 
                t(as.matrix(cterm))
            else if (!is.matrix(cterm)) 
                stop("Unexpected error: invalid data type [", 
                  class(cterm), "]")
            n <- nrow(cterm)
            ft <- coef(res)[tail(1:nbasis(res), n), , drop = FALSE]
            dimnames(ft) <- NULL
            expect_identical(cterm, ft, "Fixed coef term don't change")
        }
        if (!is.null(bterm)) {
            if (is.null(dim(bterm))) 
                bterm <- matrix(bterm, ncol = 1L)
            else if (is.data.frame(bterm)) 
                as.matrix(bterm)
            else if (!is.matrix(cterm)) 
                stop("Unexpected error: invalid data type [", 
                  class(bterm), "]")
            n <- ncol(bterm)
            ft <- basis(res)[, tail(1:nbasis(res), n), drop = FALSE]
            dimnames(ft) <- NULL
            expect_identical(bterm, ft, "Fixed basis term don't change")
        }
    }
    .check(nmf(V ~ cx), c(dim(V), 1L), cterm = cx, "Single coef term")
    .check(nmf(V ~ h), c(dim(V), nrow(h)), cterm = h, "Matrix coef term")
    .check(nmf(t(V) ~ t(w)), c(dim(t(V)), ncol(w)), cterm = t(w), 
        "Matrix coef term (transpose)")
    .check(nmf(V ~ data.frame(t(h))), c(dim(V), nrow(h)), cterm = h, 
        "Data frame coef term")
})

test_that("test.parallel", {
    r <- 3
    V <- .testData(r = r)
    ref <- nmf(V, r, nrun = 3, .opt = "k-p", seed = 123456)
    res <- nmf(V, r, nrun = 3, .opt = "kP", seed = 123456, .pbackend = "seq")
    expect_true(nmf.equal(res, ref), info = "Identical results with seeded parallel .pbackend='seq' and standard sequential '-p'")
    res <- nmf(V, r, nrun = 3, .opt = "kP", seed = 123456, .pbackend = "par")
    expect_true(nmf.equal(res, ref), info = "Identical results with seeded parallel .pbackend='par' and standard sequential '-p'")
    res <- nmf(V, r, nrun = 3, .opt = "kP", seed = 123456, .pbackend = "mc")
    expect_true(nmf.equal(res, ref), info = "Identical results with seeded parallel .pbackend='mc' and standard sequential '-p'")
    res <- nmf(V, r, nrun = 3, .pbackend = NA, seed = 123456)
    expect_true(nmf.equal(res, ref), info = "Identical results with seeded .pbackend=NA and standard sequential")
    cl <- makeCluster(2)
    on.exit(stopCluster(cl), add = TRUE)
    res <- nmf(V, r, nrun = 3, .pbackend = cl, seed = 123456)
    expect_true(nmf.equal(res, ref), info = "Identical results with seeded .pbackend=cl and standard sequential")
    registerDoParallel(cl)
    on.exit(registerDoSEQ(), add = TRUE)
    res <- nmf(V, r, nrun = 3, .pbackend = NULL, seed = 123456)
    expect_true(nmf.equal(res, ref), info = "Identical results with seeded .pbackend=NULL + registered backend and standard sequential")
})

test_that("test.registry", {
    checkNotNull <- function(x, ...) expect_true(!is.null(x), ...)
    dummy.method <- function() {
    }
    expect_identical(FALSE, removeNMFMethod("algo.tata"), info = "removeNMFMethod a method that does not exist: should not generate an error")
    expect_error(nmfAlgorithm("algo.toto"), info = "Try to access a method that does not exist: should generate an error")
    expect_true(is.null(nmfAlgorithm("algo.toto", error = FALSE)), 
        info = "Try to access a method that does not exist with error=FALSE: should NOT generate an error and return NULL")
    on.exit({
        removeNMFMethod("dummy")
    }, add = TRUE)
    checkNotNull(setNMFMethod("dummy", dummy.method), "Register works on dummy -- empty -- method")
    expect_error(setNMFMethod("dummy", dummy.method), info = "Try to register an algorithm with an existing name")
    checkNotNull(setNMFMethod("dummy", dummy.method, overwrite = TRUE), 
        "Overwrite an existing algorithm ")
    expect_true(is(nmfAlgorithm("dummy"), "NMFStrategyFunction"), 
        info = "Get method by exact match")
    expect_true(is(nmfAlgorithm("dum"), "NMFStrategyFunction"), 
        info = "Get method by partial match")
    expect_equal("dummy", name(nmfAlgorithm("dum")), info = "The method's full name is set in slot 'name'")
})

test_that("test.seed", {
    r <- 3
    V <- .testData(r = r)
    n <- nrow(V)
    m <- ncol(V)
    obj <- seed(V, r)
    check.seed("Call with rank", obj, V, r, "NMFstd")
    obj <- seed(V, r, 123456)
    check.seed("Call with rank and numeric seed", obj, V, r, 
        "NMFstd")
    obj <- seed(V, r, "nndsvd", densify = "average")
    check.seed("Call with name and extra parameters", obj, V, 
        r, "NMFstd")
    expect_error(seed(V, r, "random", toto = 1), info = "Throw an error when: unused parameter is passed to seeding method")
    class.in <- "NMFOffset"
    obj <- seed(V, list(class.in, r))
    check.seed("Call with class", obj, V, r, class.in)
    class.in <- "NMFOffset"
    obj <- nmfModel(r, model = class.in)
    obj <- seed(V, obj)
    check.seed("Call with object", obj, V, r, class.in)
    expect_error({
        obj <- seed(V, r, "seed.toto")
    }, info = "Error when calling with an undefined seeding method")
    expect_true(inherits(seed(V, r, "random"), "NMFfit"), info = "No error when calling with a defined seeding method")
})

test_that("test.summary", {
    r <- 3
    V <- .testData(r = r)
    expect_true(is.numeric(summary(nmf(V, r))), info = "Single run: result is numeric")
    expect_true(is.numeric(summary(nmf(V, r, nrun = 3))), info = "Multiple run (no keep): result is numeric")
    expect_true(is.numeric(summary(nmf(V, r, nrun = 3, .options = "k"))), 
        info = "Multiple run with keep: result is numeric")
})

