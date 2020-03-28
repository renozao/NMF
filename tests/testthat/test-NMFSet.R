# Unit testing for multiple NMF runs
# 
# Author: Renaud Gaujoux
# Converted from RUnit on 22 Feb 2020
###############################################################################

# make the internal functions visible
if( isNamespaceLoaded('NMF') ){
    join <- NMF:::NMFfitX
}

.testData <- function(n=20, r=3, m=10, ...){
    syntheticNMF(n, r, m, ...)
}

check.result <- function(obj, V, r, size, msg=NULL){
    expect_true( is(obj, 'NMFfitX'), paste(msg, ' -> result is an object of class NMFfitX', sep=''))
    expect_equal( nrun(obj), size, info = paste(msg, ' -> number of run is correctly set', sep=''))
}

test_that("test.fit", {
    set.seed(123456)
    n <- 20
    r <- 3
    m <- 30
    V <- rmatrix(n, m)
    check.result <- function(res) {
        cl <- class(res)
        expect_true(is(fit(res), "NMF"), paste(cl, "- fit: returns an object of class 'NMF'"))
        expect_true(!isNMFfit(fit(res)), paste(cl, "- fit: does not return the complete fit"))
        expect_true(is(minfit(res), "NMFfit"), paste(cl, "- minfit: returns a 'NMFfit' object"))
    }
    check.result(nmf(V, r))
    check.result(nmf(V, r, nrun = 3, .opt = "-k"))
    check.result(nmf(V, r, nrun = 3, .opt = "k"))
})

test_that("test.interface", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    rownames(V) <- seq(nrow(V))
    colnames(V) <- seq(ncol(V))
    res <- nmf(V, r, nrun = 5)
    expect_true(is.factor(predict(res)), info = "Clusters are computed without error")
    expect_equal(ncol(V), length(predict(res)), info = "Clusters for samples have correct length")
    expect_equal(nrow(V), length(predict(res, "features")), info = "Clusters for features have correct length")
    expect_true(is.list(predict(res, prob = TRUE)), info = "Clusters are computed without error")
    expect_true(is.character(featureNames(res)), info = "Feature names are accessible without error")
    expect_true(is.character(sampleNames(res)), info = "Sample names are accessible without error")
})

test_that("test.join.multipleAndSingleRunsMethods", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    resL <- list()
    nruns <- c(1, 3, 4, 1)
    resL$a <- nmf(V, r)
    resL$b <- nmf(V, r, nrun = nruns[2])
    resL$c <- nmf(V, r, nrun = nruns[3])
    resL$d <- nmf(V, r)
    res <- join(resL)
    check.result(res, V, r, sum(nruns), "Join multiple runs + single runs")
})

test_that("test.join.multipleRuns", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    resL <- list()
    nruns <- c(2, 3, 4)
    resL$brunet <- nmf(V, r, nrun = nruns[1])
    resL$brunet2 <- nmf(V, r, nrun = nruns[2])
    resL$brunet3 <- nmf(V, r, nrun = nruns[3])
    res <- join(resL)
    check.result(res, V, r, sum(nruns), "Join multiple runs")
})

test_that("test.join.singleRuns", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    resL <- list()
    resL$brunet <- nmf(V, r)
    resL$brunet2 <- nmf(V, r)
    resL$brunet3 <- nmf(V, r)
    res <- join(resL)
    check.result(res, V, r, length(resL), "Simple join")
    res <- join(resL, runtime.all = {
        tt <- runif(5)
    })
    check.result(res, V, r, length(resL), "Simple join + set time")
    expect_true(all.equal(tt, as.numeric(runtime.all(res))), 
        info = "Simple join + set time: time is correctly set")
    res <- join(resL, .merge = TRUE)
    check.result(res, V, r, length(resL), "Merging join")
})

test_that("test.multipleruns", {
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    res <- nmf(V, r, nrun = 5)
    check.result(res, V, r, 5, "Multiple runs")
})

