#' Unit Testing script for NMF package: virtual class NMFStrategy
#'
#' @author Renaud Gaujoux
#' @creation 14 Aug 2009
#' Converted from Runit 22 Feb 2020

# make the internal functions/objects visible
# if( isNamespaceLoaded('NMF') ){
#     name <- NMF:::name
#     `name<-` <- NMF:::`name<-`
#     is.mixed <- NMF:::is.mixed
# }

.TestSeed <- 123456

check.slots.methods <- function(obj, title=''){
    
    expect_equal( name(obj), obj@name, info = paste(title, ": slot methods for 'name' give same result"))
    expect_equal( objective(obj), obj@objective, info = paste(title, ": slot methods for 'objective' give same result"))
    expect_equal( modelname(obj), obj@model, info = paste(title, ": slot methods for 'model' give same result"))
    expect_equal( is.mixed(obj), obj@mixed, info = paste(title, ": slot methods for 'mixed' give same result"))
    
}

checkClass <- function(x, cl, msg){
    expect_equal(class(x)[1L], cl, info = paste(msg, "[object of class '", class(x)[1L],"']"))
}


test_that("test.accessors", {
    setClass("A", contains = "NMFStrategy")
    on.exit(removeClass("A"), add = TRUE)
    a <- new("A")
    expect_error(validObject(a), info = "Prototype object is not valid")
    check.slots.methods(a, "Prototype object")
    expect_error(name(a) <- "", info = "Method name<-: Error if setting slot 'name' to ''")
    expect_error(name(a) <- character(), info = "Method name<-: Error if setting slot 'name' to character()")
    expect_error(name(a) <- 4, info = "Method name<-: Error if setting slot 'name' to not character")
    name(a) <- "toto"
    expect_equal("toto", a@name, info = "Method name<-: set slot 'name' correctly")
    check.slots.methods(a, "Object after name set by name<-")
    expect_error(objective(a) <- "", info = "Method objective<-: Error if setting slot 'objective' to ''")
    expect_error(objective(a) <- character(), info = "Method objective<-: Error if setting slot 'objective' to character()")
    expect_error(objective(a) <- 4, info = "Method objective<-: Error if setting slot 'objective' to not character")
    objective(a) <- "toto"
    expect_equal("toto", a@objective, info = "Method objective<-: set slot 'objective' correctly")
    check.slots.methods(a, "Object after name set by objective<-")
    objective(a) <- function(...) {
    }
    expect_equal(function(...) {
    }, a@objective, info = "Method objective<-: set slot 'objective' correctly")
    check.slots.methods(a, "Object after name set by objective<-")
})

test_that("test.constructor", {
    expect_error(new("NMFStrategy"), info = "Class is virtual")
    setClass("A", contains = "NMFStrategy")
    on.exit(removeClass("A"), add = TRUE)
    expect_error(new("A", name = ""), info = "Error if slot name is an empty character string")
    expect_true(validObject(new("A", name = "toto")), info = "No error if object with non empty name")
    expect_error(new("A", name = "toto", objective = 4), info = "Error if slot objective is NOT a character or function (numeric)")
    expect_error(new("A", name = "toto", objective = ""), info = "Error if slot objective is an empty character string")
    expect_true(validObject(new("A", name = "toto", objective = "tata")), 
        info = "No error if slot objective is a non empty character string")
    expect_true(validObject(new("A", name = "toto", objective = function() {
    })), info = "No error if slot objective is a function")
    expect_error(new("A", name = "toto", model = ""), info = "Error if slot model is an empty character string")
    expect_error(new("A", name = "toto", model = "toto"), info = "Error if slot model is not an sub-class of class NMF")
    expect_true(validObject(new("A", name = "toto", model = "NMFstd")), 
        info = "Valid object if slot model is 'NMFstd'")
    expect_true(validObject(new("A", name = "toto", model = "NMFns")), 
        info = "Valid object if slot model is a subclass of class 'NMF'")
    expect_true(validObject(new("A", name = "toto", model = c("NMFns", 
        "NMFOffset"))), info = "Valid object if slot model is vector of subclasses of class 'NMF'")
    expect_error(new("A", name = "toto", mixed = "toto"), info = "Error if slot mixed is not logical")
    expect_error(new("A", name = "toto", mixed = c(TRUE, FALSE)), 
        info = "Error if slot mixed is not length 1")
    expect_true(validObject(new("A", name = "toto", mixed = TRUE)), 
        info = "Valid object if slot mixed is TRUE")
    expect_true(validObject(new("A", name = "toto", mixed = TRUE)), 
        info = "Valid object if slot mixed is FALSE")
})

test_that("test.constructorMethod", {
    checkClass(NMFStrategy("a", Update = function(i, y, x, ...) {
    }), "NMFStrategyIterative", "With argument `Update`")
    checkClass(NMFStrategy("a", function(i, y, x, ...) {
    }), "NMFStrategyFunction", "With method=function")
    checkClass(NMFStrategy("a", algorithm = function(i, y, x, 
        ...) {
    }), "NMFStrategyFunction", "With argument `algorithm`")
})

