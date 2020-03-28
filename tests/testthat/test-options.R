#' Unit Testing script for NMF package: package specific options
#'
#' @author Renaud Gaujoux
#' @creation 25 July 2009
#' Converted from RUnit on 22 Feb 2020

test_that("test.nmf.getOption", {
    OLD <- nmf.options()
    on.exit(nmf.options(OLD))
    expect_true(is.null(nmf.getOption("toto")), info = "Unknow in nmf.getOption returns NULL")
    nmf.options(toto = 5)
    expect_equal(5, nmf.getOption("toto"), info = "nmf.getOption returns correct value")
    nmf.options(toto = NULL)
})

test_that("test.nmf.options", {
    OLD <- nmf.options()
    on.exit(nmf.options(OLD))
    expect_true(is.list(nmf.options()), info = "Options are returned as a list")
    expect_true(is.list(nmf.options("error.track", "debug")), 
        info = "Options are returned as a list")
    expect_true(is.list(nmf.options("error.track")), info = "Single option is returned as a list")
    expect_equal(list(toto = NULL), nmf.options("toto"), info = "Unknow in nmf.options returns NULL correctly named")
    expect_false(is.null(nmf.options(toto = 6)), info = "Can add new option")
    nmf.options(tata = 10)
    opt <- nmf.options()
    expect_equal(list(toto = 6, tata = 10, titi = NULL), o <- nmf.options(toto = 5, 
        tata = 9, titi = 25), info = "Changing options return old values")
    expect_identical(opt, {
        nmf.options(o)
        nmf.options()
    }, info = "Restoring options works")
    nmf.options(toto = 4)
    expect_equal(list(toto = 4), nmf.options(toto = NULL), info = "Removing an option returns old value")
    expect_false(is.element("toto", names(nmf.options())), info = "Removing an option actually removes it from the option list")
})

