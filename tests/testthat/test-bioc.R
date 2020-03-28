# Unit tests for the Bioconductor layer
# 
# Author: Renaud Gaujoux
# Created: Mar 6, 2013
# Converted from RUnit: 16 Feb 2020
###############################################################################


# check extension of Bioc access methods
# => this also serves to check that they are correctly exported
test_that("test.access", {
    x <- nmfModel(3, 20, 10)
    .check <- function(fun, val, newval) {
        msg <- function(...) paste(fun, ": ", ..., sep = "")
        f <- match.fun(fun)
        expect_identical(f(x), val, msg("on fresh object returns NULL"))
        if (isNumber(newval)) 
            newval <- paste("aaa", 1:newval)
        res <- try(eval(parse(text = paste(fun, "(x) <- newval", 
            sep = ""))))
        expect_true(!is(res, "try-error"), msg("setting value works"))
        expect_identical(f(x), newval, msg("new value is correct"))
        res <- try(eval(parse(text = paste(fun, "(x) <- val", 
            sep = ""))))
        expect_true(!is(res, "try-error"), msg("resetting value works"))
        expect_identical(f(x), val, msg("reset value is correct"))
    }
    expect_identical(nmeta(x), nbasis(x), info = "nmeta is defined")
    .check("featureNames", NULL, 20)
    .check("sampleNames", NULL, 10)
    .check("basisnames", NULL, 3)
    .check("metaprofiles", coef(x), rmatrix(coef(x)))
    .check("metagenes", basis(x), rmatrix(basis(x)))
})

