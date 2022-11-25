#' Unit Testing script for NMF package: NMF models and utilities.
#'
#' @author Renaud Gaujoux
#' @creation 22 April 2009

# make the internal functions/objects visible
if( isNamespaceLoaded('NMF') ){
    .predict.nmf <- NMF:::.predict.nmf
    is.same <- NMF:::is.same
}

.TestSeed <- 123456

checkOnlyNAs <- function(x, msg){
    expect_true( all(is.na(basis(x))) , paste(msg, ': Only NAs in W'))	
    expect_true( all(is.na(coef(x))), paste(msg, ': Only NAs in H'))	
}

basicHM <- aheatmap

checkPlot <- function(...){
    if( isCHECK() ) return()
    pkgmaker::checkPlot(...)
}

#' checks the validity of a NMF object: dimensions
check.object <- function(obj, n, m, r, title='Check NMF object', class='NMFstd'){
    
    # convert dimensions to integers
    n <- as.integer(n)
    m <- as.integer(m)
    r <- as.integer(r)
    
    # check class
    expect_true(is(obj, class), paste(title, ': object is of class "', class, "'"));
    # check virtual interface accessors
    expect_true(is.same(obj@W, basis(obj)), msg('method `basis` returns slot W'))
    expect_true(is.same(obj@H, coef(obj)), msg('method `coef` returns slot H'))
    
    expect_identical(nrow(basis(obj)), n, paste(title, ': number of rows in basis matrix is correct'))
    expect_identical(nrow(obj), n, paste(title, ': nrow returns correct value'))
    expect_identical(ncol(basis(obj)), r, paste(title, ': number of columns in basis matrix is correct'))
    expect_identical(nbasis(obj), r, paste(title, ': nbasis returns correct value'))
    expect_identical(nrow(coef(obj)), r, paste(title, ': number of rows in coef matrix correct'))
    expect_identical(ncol(coef(obj)), m, paste(title, ': number of columns in coef matrix is correct'))
    expect_identical(ncol(obj), m, paste(title, ': ncol returns correct value'))
}


check.NMF.class <- function(class, ...){
    
    set.seed(.TestSeed)
    
    msg <- function(...) paste('New -', class,':', ...) 
    
    # base constructor: should build an empty model
    a <- new(class, ...)
    
    # check if there is only NAs
    checkOnlyNAs(a, msg('with no parameters'))	
    check.object(a, 0, 0, 0, msg('with no parameters'), class)
    
    # define some dimensions to use as template
    n <- 50; m <- 10; r <- 3;
    w <- rmatrix(n, r); h <- rmatrix(r, m)
    
    # base constructor with one of the two matrices: exceptions	
    expect_error(new(class, W=w, ...), "W and H are not compatible", info = msg('error if only basis matrix'))
    expect_error(new(class, H=h, ...), "W and H are not compatible", info = msg('error if only coef matrix'))	
    
    # base constructor with two matrices		
    a <- new(class, W=w, H=h, ...)
    
    # check the dimensions of the internal matrices
    check.object(a, n, m, r, msg('with two matrices - '), class)
    expect_identical(basis(a), w, msg('with two matrices - basis matrix is correctly set')); 
    expect_identical(coef(a), h, msg('with two matrices - coef matrix is correctly set'));
    
    # check error with wrong dimension
    expect_error(a <- new(class, W=matrix(1,n,r+1), H=matrix(1,r,m), ...), 
                 "W and H are not compatible",
                 info = msg('Error if incompatible dimensions (W)'))
    expect_error(a <- new(class, W=matrix(1,n,r), H=matrix(1,r+1,m), ...), 
                 "W and H are not compatible", 
                 info = msg('Error if incompatible dimensions (H)'))	
}

check.dimnames <- function(x, dn, msg){
    
    expect_identical(dimnames(x), dn, paste(msg, '-',"dimnames returns correct value"))
    expect_identical(dimnames(x)[c(1,3)], dimnames(basis(x)), paste(msg, '-', "dimnames returns value consistent with basis"))
    expect_identical(dimnames(x)[c(3,2)], dimnames(coef(x)), paste(msg, '-', "dimnames returns value consistent with coef"))
    
    expect_identical(rownames(x), dn[[1]], paste(msg, '-', "rownames returns correct value"))
    expect_identical(rownames(basis(x)), rownames(x), paste(msg, '-', "rownames returns same value as rownames(basis)"))
    
    expect_identical(colnames(x), dn[[2]], paste(msg, '-', "colnames returns correct value"))
    expect_identical(colnames(coef(x)), colnames(x), paste(msg, '-', "colnames returns same value as colnames(basis)"))
    
    expect_identical(basisnames(x), dn[[3]], paste(msg, '-', "basisnames returns correct value"))
    expect_identical(colnames(basis(x)), basisnames(x), paste(msg, '-', "basisnames returns same value as colnames(basis)"))
    expect_identical(rownames(coef(x)), basisnames(x), paste(msg, '-', "basisnames returns same value as rownames(coef)"))
    
}

test_that("test.basis", {
    n <- 50
    m <- 10
    r <- 3
    a <- nmfModel(r)
    expect_true(is.same(basis(a), a@W), info = "Method 'basis' correctly returns slot W")
    W.ext <- matrix(1, n, r)
    basis(a) <- W.ext
    expect_identical(W.ext, a@W, info = "Method 'basis<-' correctly sets slot W")
})

test_that("test.class.NMFns", {
    check.NMF.class("NMFns")
    t <- 0.8
    check.NMF.class("NMFns", theta = t)
    a <- new("NMFns", theta = t)
    expect_identical(t, a@theta, info = "Slot theta is correctly set")
    expect_error(a <- new("NMFns", theta = -1), info = "Negative value of theta throws an exception")
    expect_error(a <- new("NMFns", theta = 1.2), info = "Value of theta > 1 throws an exception")
    set.seed(.TestSeed)
    n <- 50
    m <- 10
    r <- 3
    W <- rmatrix(n, r)
    H <- rmatrix(r, m)
    a.ns <- new("NMFns", W = W, H = H, theta = 0)
    a.std <- new("NMFstd", W = W, H = H)
    expect_identical(fitted(a.std), fitted(a.ns), info = "Values fitted correspond to standard model if theta=0")
    a.ns <- new("NMFns", W = W, H = H, theta = 0.4)
    s <- smoothing(a.ns)
    expect_equal(c(r, r), dim(s), info = "Smoothing matrix: dimension are correct")
    expect_true(all(s >= 0), info = "Smoothing matrix: all entries are nonnegative")
    expect_equal(rep(1, r), rowSums(s), info = "Smoothing matrix: sum of rows are ones")
    expect_equal(rep(1, r), colSums(s), info = "Smoothing matrix: sum of columns are ones")
    expect_equal(basis(a.ns) %*% s %*% coef(a.ns), fitted(a.ns), 
        info = "Fitted values are correct (product of basis, smoothing and coef).")
})

test_that("test.class.NMFstd", {
    check.NMF.class("NMFstd")
})

test_that("test.coef", {
    n <- 50
    m <- 10
    r <- 3
    a <- nmfModel(r)
    expect_true(is.same(coef(a), a@H), info = "Method 'coef' correctly returns slot H")
    ext <- matrix(1, r, m)
    coef(a) <- ext
    expect_identical(ext, a@H, info = "Method 'coef<-' correctly sets slot H")
})

test_that("test.connectivity", {
    n <- 50
    m <- 10
    r <- 3
    set.seed(.TestSeed)
    a <- nmfModel(r, c(n, m))
    a <- rnmf(a)
    con <- connectivity(a)
    expect_true(is.matrix(con), info = "The result is a matrix")
    expect_true(all(con %in% c(0, 1)), info = "All entries are 0 or 1")
    expect_true(all(t(con) == con), info = "The connectivity matrix is symmetric")
})

test_that("test.deviance", {
    n <- 10
    m <- 5
    set.seed(.TestSeed)
    y <- rmatrix(n, m)
    x <- rnmf(nmfModel(3, y))
    expect_warning(dev <- deviance(x, y), "Undefined distance method: .* [returned NA]")
    expect_equal(NA_real_, dev, info = "With no method: NA")
    expect_error(deviance(x, y, "toto"), info = "Error if undefined method")
    norm1 <- function(x, y) {
        sum(abs(fitted(x) - y))
    }
    funname <- paste("test", paste(sample(1:9, 20, replace = TRUE), 
        collapse = ""), sep = ".")
    assign(funname, norm1, envir = .GlobalEnv)
    on.exit(rm(list = funname, envir = .GlobalEnv), add = TRUE)
    expect_true(deviance(x, y, norm1) > 0, info = "Works with a user-defined function in .GlobalEnv")
    expect_true(deviance(x, y, funname) > 0, info = "Works with name of a user-defined function .GlobalEnv")
    meth <- "euclidean"
    expect_true(deviance(x, y, meth) > 0, info = "Euclidean: positive")
    expect_equal(sum((fitted(x) - y)^2)/2, deviance(x, y, meth), 
        info = "Euclidean: OK")
    meth <- "KL"
    expect_true(deviance(x, y, meth) > 0, info = "KL: positive")
    z <- y
    z[1, ] <- 0
    expect_true(deviance(x, z, meth) != Inf, info = "Ok if some zeros in the target")
    z <- x
    basis(z)[1, ] <- 0
    expect_identical(Inf, deviance(z, y, meth), info = "Infinite if some zeros in the estimate")
})

test_that("test.dimensions", {
    n <- as.integer(50)
    m <- as.integer(10)
    r <- as.integer(3)
    a <- nmfModel(r, n, m)
    expect_identical(c(n, m, r), dim(a), info = "Function 'dim' is OK")
    expect_identical(n, nrow(a), info = "Function 'nrow' is OK")
    expect_identical(m, ncol(a), info = "Function 'ncol' is OK")
    expect_identical(r, nbasis(a), info = "Function 'nbasis' is OK")
})

test_that("test.dimnames", {
    set.seed(.TestSeed)
    n <- 20
    m <- 10
    r <- 3
    w <- rmatrix(n, r)
    h <- rmatrix(r, m)
    M <- nmfModel(r, n, m)
    a <- M
    expect_error({
        dimnames(a) <- 1:n
    }, info = "set to vector")
    expect_error({
        dimnames(a) <- list(seq(n - 1))
    }, info = "Error of wrong dimension (nrow-1)")
    expect_error({
        dimnames(a) <- list(seq(n + 1))
    }, info = "Error of wrong dimension (nrow+1)")
    msg <- function(...) paste("Dimnames with 2 elements -", 
        ...)
    a <- M
    check.dimnames(a, NULL, "No dimnames => NULL")
    check.dimnames({
        dimnames(a) <- NULL
        a
    }, NULL, msg("set to NULL"))
    check.dimnames({
        dimnames(a) <- list()
        a
    }, NULL, msg("set to list()"))
    msg <- function(...) paste("Dimnames with 1 element -", ...)
    a <- M
    dn <- list(letters[1:nrow(a)])
    dn.name <- setNames(dn, "rows")
    check.dimnames({
        a <- M
        dimnames(a) <- dn
        a
    }, c(dn, list(NULL, NULL)), msg("Set dimnames"))
    check.dimnames({
        a <- M
        dimnames(a) <- dn.name
        a
    }, c(dn.name, list(NULL, NULL)), msg("Set with names"))
    check.dimnames({
        a <- M
        rownames(a) <- dn[[1]]
        a
    }, c(dn, list(NULL, NULL)), msg("Set rownames"))
    check.dimnames({
        a <- M
        colnames(a) <- letters[1:ncol(a)]
        a
    }, list(NULL, letters[1:ncol(a)], NULL), msg("Set colnames"))
    check.dimnames({
        a <- M
        basisnames(a) <- letters[1:nbasis(a)]
        a
    }, list(NULL, NULL, letters[1:nbasis(a)]), msg("Set basisnames"))
    check.dimnames({
        dimnames(a) <- NULL
        a
    }, NULL, msg("Reset to NULL"))
    msg <- function(...) paste("Dimnames with 2 elements -", 
        ...)
    a <- M
    dn <- list(letters[1:nrow(a)], letters[seq(nrow(a) + 1, nrow(a) + 
        ncol(a))])
    dn.name <- setNames(dn, c("rows", "cols"))
    check.dimnames({
        dimnames(a) <- dn
        a
    }, c(dn, list(NULL)), msg("Set dimnames"))
    check.dimnames({
        dimnames(a) <- NULL
        a
    }, NULL, msg("Reset to NULL"))
    check.dimnames({
        dimnames(a) <- dn.name
        a
    }, c(dn.name, list(NULL)), msg("Set with names"))
    check.dimnames({
        dimnames(a) <- NULL
        a
    }, NULL, msg("Reset to NULL (2)"))
    msg <- function(...) paste("Dimnames with 3 elements -", 
        ...)
    a <- M
    dn <- list(letters[1:nrow(a)], letters[seq(nrow(a) + 1, nrow(a) + 
        ncol(a))], letters[seq(nrow(a) + ncol(a) + 1, nrow(a) + 
        ncol(a) + nbasis(a))])
    dn.name <- setNames(dn, c("rows", "cols", "basis"))
    check.dimnames({
        dimnames(a) <- dn
        a
    }, dn, msg("Set dimnames"))
    check.dimnames({
        dimnames(a) <- NULL
        a
    }, NULL, msg("Reset to NULL"))
    check.dimnames({
        dimnames(a) <- dn.name
        a
    }, dn.name, msg("Set with names"))
    check.dimnames({
        dimnames(a) <- NULL
        a
    }, NULL, msg("Reset to NULL (2)"))
})

test_that("test.entropy", {
    checkBounds <- function(x) {
        expect_true(x >= 0, "greater than 0")
        expect_true(x <= 1, "lower than 1")
    }
    x <- as.factor(c(rep(1, 5), rep(2, 10), rep(3, 15)))
    e <- entropy(x, x)
    checkBounds(e)
    expect_equal(0, e)
    set.seed(.TestSeed)
    e <- entropy(as.factor(sample(x)), x)
    checkBounds(e)
})

test_that("test.misc", {
    x <- nmfModel()
    m <- slot(x, "misc")
    expect_true(is.list(m) && length(m) == 0L, info = "On empty model misc is an empty list")
    expect_equal(m, misc(x), info = "On empty model misc() returns an empty list")
    x$a <- 3
    expect_equal(list(a = 3), slot(x, "misc"), info = "Setting misc with $ works")
    expect_equal(list(a = 3), misc(x), info = "Getting misc with misc() works")
    expect_equal(3, x$a, info = "Getting misc with $ works")
    expect_equal(NULL, misc(list()), info = "On empty list misc is NULL")
    expect_equal(4, misc(list(misc = 4)), info = "On list with a `misc` element, misc() returns the element")
    expect_equal(NULL, misc(1), info = "On non list object misc() is NULL")
    a <- 1
    attr(a, "misc") <- 2
    expect_equal(2, misc(a), info = "On non list object with `misc` attribute, misc() returns the attribute")
})

test_that("test.NMF.rnmf", {
    n <- 50
    m <- 10
    r <- 3
    a <- nmfModel(r)
    set.seed(.TestSeed)
    expect_error(rnmf(a, as.numeric(NA)), info = "Error thrown if single NA value")
    expect_error(rnmf(a, as.numeric(c(1, NA))), info = "Error thrown if some NA value")
    expect_error(rnmf(a, c(1, 2, 3)), info = "Error thrown if some length greater than 2")
    expect_error(rnmf(a, numeric()), info = "Error thrown if some length is 0")
    set.seed(.TestSeed)
    a <- rnmf(r, n, m)
    checkPlot(basicHM(basis(a)), "Random NMF basis (target dimension)")
    checkPlot(basicHM(coef(a)), "Random NMF coef (target dimension)")
    expect_equal(n, nrow(a), info = "number of rows matches target dimension")
    expect_equal(m, ncol(a), info = "number of columns matches target dimension")
    expect_equal(r, nbasis(a), info = "number of basis matches target dimension")
    set.seed(.TestSeed)
    expect_identical(rnmf(r, c(n, m)), a, info = "calling with numeric target of length 2 is equivalent to separate dimensions")
    set.seed(.TestSeed)
    expect_identical(rnmf(r, n, m), a, info = "calling with numeric target of length 2 is equivalent to separate dimensions")
    a <- nmfModel(r, n, m)
    set.seed(.TestSeed)
    b <- rnmf(a)
    checkPlot(basicHM(basis(b)), "Random NMF basis (target NMF)")
    checkPlot(basicHM(coef(b)), "Random NMF coef (target NMF)")
    expect_equal(nrow(b), nrow(a), info = "number of rows matches NMF target")
    expect_equal(ncol(b), ncol(a), info = "number of columns matches NMF target")
    expect_equal(nbasis(b), nbasis(a), info = "number of basis matches NMF model")
    set.seed(.TestSeed)
    max.entry <- 100
    V <- ExpressionSet(rmatrix(n, m, min = 0, max = max.entry))
    a <- nmfModel(r, n, m)
    check.object(a <- rnmf(a, V), n, m, r, "NMFobject + target ExpressionSet")
    check.dimnames(a, c(dimnames(exprs(V)), list(NULL)), "NMFobject + target ExpressionSet")
    check.object(a <- rnmf(r, V), n, m, r, "rank + target ExpressionSet")
    check.dimnames(a, c(dimnames(exprs(V)), list(NULL)), "rank + target ExpressionSet")
    set.seed(.TestSeed)
    max.entry <- 100
    V <- rmatrix(n, m, min = 0, max = max.entry)
    a <- nmfModel(r, n, m)
    a <- rnmf(a, V)
    checkPlot(basicHM(basis(a)), "Random NMF basis (target matrix)")
    checkPlot(basicHM(coef(a)), "Random NMF coef (target matrix)")
    expect_equal(n, nrow(a), info = "number of rows matches matrix target")
    expect_equal(m, ncol(a), info = "number of columns matches matrix target")
    expect_equal(r, nbasis(a), info = "matrix target: number of basis matches NMF model")
    msg <- function(...) paste("Set max in target matrix -", 
        ...)
    .check_max <- function(x, max.entry) {
        expect_true(max(basis(a)) <= max.entry, msg("Basis maximum entry is OK (<=)"))
        expect_true(max(basis(a)) >= max.entry/2, msg("Basis maximum entry is OK (>=)"))
        expect_true(max(coef(a)) <= max.entry, msg("Coef maximum entry is OK (<=)"))
        expect_true(max(coef(a)) >= max.entry/2, msg("Coef maximum entry is OK (>=)"))
    }
    .check_max(a, max.entry)
    dn <- list(rows = letters[1:nrow(V)], cols = letters[seq(25 - 
        ncol(V) + 1, 25)])
    dimnames(V) <- dn
    check.dimnames(rnmf(a, V), c(dn, list(NULL)), "rnmf on target matrix with dimnames: dimnames are passed")
    check.dimnames(rnmf(a, V, use.dimnames = FALSE), NULL, "rnmf on target matrix with dimnames and use.dimnames=FALSE: dimnames are not passed")
    msg <- function(...) paste("Set max by argument -", ...)
    set.seed(.TestSeed)
    a <- rnmf(r, V, model = "NMFOffset")
    .check_max(a, max.entry)
    msg <- function(...) paste("Set max by argument -", ...)
    set.seed(.TestSeed)
    max.entry <- 5
    a <- rnmf(r, n, m, dist = list(max = max.entry))
    .check_max(a, max.entry)
})

test_that("test.nmfModel", {
    set.seed(.TestSeed)
    n <- as.integer(25)
    m <- as.integer(10)
    r <- as.integer(3)
    check.empty.model <- function(x, n, m, r, msg, class = "NMFstd") {
        check.object(x, n, m, r, msg, class)
        checkOnlyNAs(x, msg)
    }
    expect_error(nmfModel(numeric()), info = "Error if negative rank")
    expect_error(nmfModel(c(1, 2)), info = "Error if rank of length != 1")
    expect_error(nmfModel(r, -1), info = "Error if negative target dimension")
    expect_error(nmfModel(r, 1:3), info = "Error if target dimension of length > 2")
    expect_error(nmfModel(r, 1, -1), info = "Error if target ncol negative")
    expect_error(nmfModel(r, 1, 1:2), info = "Error if target ncol of length > 1")
    expect_error(nmfModel(r, 1, matrix(1, 2, 2)), info = "Error if target ncol not vector")
    check.empty.model(nmfModel(), 0, 0, 0, "Constructor with no arguments returns empty model")
    check.empty.model(nmfModel(ncol = 5), 0, 5, 0, "Constructor with only with ncol specified")
    check.empty.model(nmfModel(r), 0, 0, r, "Constructor with missing target")
    check.empty.model(nmfModel(r, c(n, m)), n, m, r, "Constructor with dimensions as a vector")
    check.empty.model(nmfModel(r, n, m), n, m, r, "Constructor with separate dimensions")
    check.empty.model(nmfModel(r, n), n, n, r, "Constructor with single dimension (nrow)")
    check.empty.model(nmfModel(r, ncol = m), 0, m, r, "Constructor with single dimension (ncol)")
    msg <- function(...) paste("Constructor with target matrix -", 
        ...)
    V <- rmatrix(n, m)
    check.empty.model(nmfModel(r, V), n, m, r, msg("second argument"))
    check.empty.model(nmfModel(V, r), n, m, r, msg("first argument"))
    check.empty.model(nmfModel(V), n, m, 0, msg("single argument"))
    dimnames(V) <- list(rows = letters[1:n], cols = letters[1:m])
    expect_identical(c(dimnames(V), list(NULL)), dimnames(nmfModel(V)), 
        info = msg("dimnames are correctly set"))
    expect_identical(NULL, dimnames(nmfModel(V, use.names = FALSE)), 
        info = msg("dimnames are not used if use.names=FALSE"))
    w <- rmatrix(n, r)
    h <- rmatrix(r, m)
    msg <- function(...) paste("Constructor with target rank and both basis and coef matrices -", 
        ...)
    a <- nmfModel(r, W = w, H = h)
    check.object(a, n, m, r, msg())
    expect_identical(w, basis(a), info = msg("basis matrix is correctly set"))
    expect_identical(h, coef(a), info = msg("coef matrix is correctly set"))
    expect_error(nmfModel(r, c(n + 1, m), W = w, H = h), info = msg("error if bad number of rows in basis matrix"))
    expect_error(nmfModel(r, c(n, m), W = w[, -r], H = h), info = msg("error if bad number of columns in basis matrix"))
    expect_error(nmfModel(r, c(n, m), W = w, H = h[-r, ]), info = msg("error if bad number of rows in coef matrix"))
    expect_error(nmfModel(r, c(n, m + 1), W = w, H = h), info = msg("error if bad number of columns in coef matrix"))
    rmsg <- function(...) msg("reduce rank -", ...)
    expect_warning(a <- nmfModel(r - 1, W = w, H = h), 
                   "(only the first 2 columns of W will be used|only the first 2 rows of H will be used)",
                   all = TRUE)
    check.object(a, n, m, r - 1, rmsg("dimensions are OK"))
    expect_identical(w[, -r], basis(a), info = rmsg("entries for basis are OK"))
    expect_identical(h[-r, ], coef(a), info = rmsg("entries for coef are OK"))
    rmsg <- function(...) msg("reduce nrow -", ...)
    expect_warning(a <- nmfModel(r, n - 1, W = w, H = h), "only the first 24 rows of W will be used")
    check.object(a, n - 1, m, r, rmsg("dimensions are OK"))
    expect_identical(w[-n, ], basis(a), info = rmsg("entries for basis are OK"))
    expect_identical(h, coef(a), info = rmsg("entries for coef are OK"))
    rmsg <- function(...) msg("reduce ncol -", ...)
    expect_warning(a <- nmfModel(r, ncol = m - 1, W = w, H = h), "only the first 9 columns of H will be used")
    check.object(a, n, m - 1, r, rmsg("dimensions are OK"))
    expect_identical(w, basis(a), info = rmsg("entries for basis are OK"))
    expect_identical(h[, -m], coef(a), info = rmsg("entries for coef are OK"))
    msg <- function(...) paste("Constructor with only basis and coef matrices (named arguments) -", 
        ...)
    a <- nmfModel(W = w, H = h)
    check.object(a, n, m, r, msg())
    expect_identical(w, basis(a), info = msg("basis matrix is correctly set"))
    expect_identical(h, coef(a), info = msg("coef matrix is correctly set"))
    expect_error(nmfModel(W = matrix(1, n, r + 1), H = matrix(1, 
        r, m)), info = msg("error if incompatible dimensions (1)"))
    expect_error(nmfModel(W = matrix(1, n, r), H = matrix(1, 
        r + 1, m)), info = msg("error if incompatible dimensions (2)"))
    msg <- function(...) paste("Constructor with only basis and coef matrices (unamed argument) -", 
        ...)
    a <- nmfModel(w, h)
    check.object(a, n, m, r, msg())
    expect_identical(w, basis(a), info = msg("basis matrix is correctly set"))
    expect_identical(h, coef(a), info = msg("coef matrix is correctly set"))
    expect_error(nmfModel(matrix(1, n, r + 1), matrix(1, r, m)), 
        info = msg("error if incompatible dimensions (1)"))
    expect_error(nmfModel(matrix(1, n, r), matrix(1, r + 1, m)), 
        info = msg("error if incompatible dimensions (2)"))
    msg <- function(...) paste("Constructor with target rank and basis matrix only -", 
        ...)
    a <- nmfModel(r, W = w)
    check.object(a, n, 0, r, msg())
    expect_identical(w, basis(a), info = msg("basis matrix is correctly set"))
    expect_true(all(is.na(coef(a))), info = msg("only NA in coef matrix"))
    expect_error(nmfModel(r + 1, W = w), info = msg("error if smaller number of columns in basis matrix"))
    rmsg <- function(...) msg("reduce rank -", ...)
    expect_warning(a <- nmfModel(r - 1, W = w), "only the first 2 columns of W will be used")
    check.object(a, n, 0, r - 1, rmsg("dimensions are OK"))
    expect_identical(w[, -r], basis(a), info = rmsg("entries for basis are OK"))
    expect_true(all(is.na(coef(a))), info = rmsg("only NA in coef matrix"))
    expect_error(nmfModel(r - 1, W = w, force.dim = FALSE), info = msg("error if greater number of columns in basis matrix and force.dim=FALSE"))
    msg <- function(...) paste("Constructor with basis matrix only -", 
        ...)
    a <- nmfModel(W = w)
    check.object(a, n, 0, r, msg())
    expect_identical(w, basis(a), info = msg("basis matrix is correctly set"))
    expect_true(all(is.na(coef(a))), info = msg("only NA in coef matrix"))
    msg <- function(...) paste("Constructor with target rank and coef matrix only -", 
        ...)
    a <- nmfModel(r, H = h)
    check.object(a, 0, m, r, msg())
    expect_identical(h, coef(a), info = msg("coef matrix is correctly set"))
    expect_true(all(is.na(basis(a))), info = msg("only NA in basis matrix"))
    expect_error(nmfModel(r + 1, H = h), info = msg("error if smaller number of rows in coef matrix"))
    rmsg <- function(...) msg("reduce rank -", ...)
    expect_warning(a <- nmfModel(r - 1, H = h), "only the first 2 rows of H will be used")
    check.object(a, 0, m, r - 1, rmsg("dimensions are OK"))
    expect_identical(h[-r, ], coef(a), info = rmsg("coef matrix is correctly set"))
    expect_true(all(is.na(basis(a))), info = rmsg("only NA in basis matrix"))
    expect_error(nmfModel(r - 1, H = h, force.dim = FALSE), info = msg("error if greater number of rows in coef matrix and force.dim=FALSE"))
    msg <- function(...) paste("Constructor with coef matrix only -", 
        ...)
    a <- nmfModel(H = h)
    check.object(a, 0, m, r, msg())
    expect_identical(h, coef(a), info = msg("coef matrix is correctly set"))
    expect_true(all(is.na(basis(a))), info = msg("only NA in basis matrix"))
    msg <- function(...) paste("Constructor with basis matrix and both target dimensions -", 
        ...)
    a <- nmfModel(r, n, m, W = w)
    check.object(a, n, m, r, msg())
    expect_identical(w, basis(a), info = rmsg("entries for basis are OK"))
    expect_true(all(is.na(coef(a))), info = rmsg("only NA in coef matrix"))
    rmsg <- function(...) msg("reduce nrow -", ...)
    expect_warning(a <- nmfModel(r, n - 1, m, W = w), "rows in target is lower than the number of rows in W")
    check.object(a, n - 1, m, r, rmsg("dimensions are OK"))
    expect_identical(w[-n, ], basis(a), info = rmsg("entries for basis are OK"))
    expect_true(all(is.na(coef(a))), info = rmsg("only NA in coef matrix"))
    expect_error(nmfModel(r, n + 1, m, W = w), info = msg("error if smaller number of rows in basis matrix"))
    expect_error(nmfModel(r, n - 1, W = w, force.dim = FALSE), 
        info = msg("error if greater number of rows in basis matrix and force.dim=FALSE"))
    expect_error(nmfModel(r + 1, n, m, W = w), info = msg("error if smaller number of columns in basis matrix"))
    msg <- function(...) paste("Constructor with coef matrix and both target dimensions -", 
        ...)
    a <- nmfModel(r, n, m, H = h)
    check.object(a, n, m, r, msg())
    expect_true(all(is.na(basis(a))), info = msg("only NA in basis matrix"))
    expect_identical(h, coef(a), info = msg("coef matrix is correctly set"))
    rmsg <- function(...) msg("reduce ncol -", ...)
    expect_warning(a <- nmfModel(r, n, m - 1, H = h), "columns in target is lower than the number of columns in H")
    check.object(a, n, m - 1, r, rmsg("dimensions are OK"))
    expect_identical(h[, -m], coef(a), info = rmsg("coef matrix is correctly set"))
    expect_true(all(is.na(basis(a))), info = rmsg("only NA in basis matrix"))
    expect_error(nmfModel(r + 1, n, m, H = h), info = msg("error if smaller number of rows in coef matrix"))
    expect_error(nmfModel(r, n, m - 1, H = h, force.dim = FALSE), 
        info = msg("error if greater number of columns in coef matrix and force.dim=FALSE"))
    expect_error(nmfModel(r, n, m + 1, H = h), info = msg("error if smaller number of columns in coef matrix"))
    check.model.dimnames <- function(x, dn, title) {
        msg <- function(...) paste(title, "-", ...)
        expect_identical(dimnames(x), dn, msg("dimnames are correct"))
        expect_identical(colnames(basis(x)), dn[[3]], msg("colnames of basis matrix are correct"))
        expect_identical(rownames(coef(x)), dn[[3]], msg("rownames of coef matrix are correct"))
    }
    dn <- letters[1:nrow(h)]
    h2 <- h
    rownames(h2) <- dn
    a <- nmfModel(r, H = h2)
    check.model.dimnames(a, list(NULL, NULL, dn), "Basis names are passed from input coef matrix")
    w2 <- w
    colnames(w2) <- dn
    a <- nmfModel(r, W = w2)
    check.model.dimnames(a, list(NULL, NULL, dn), "Basis names are passed from input basis matrix")
    w2 <- w
    colnames(w2) <- dn
    h2 <- h
    rownames(h2) <- dn
    a <- nmfModel(W = w2, H = h2)
    check.model.dimnames(a, list(NULL, NULL, dn), "Basis names are used unchanged if equal in input basis and coef matrices")
    msg <- function(...) paste("Basis names from input basis matrix are used to order the components - ", 
        ...)
    w2 <- w
    colnames(w2) <- dn
    h2 <- h
    rownames(h2) <- rev(dn)
    a <- nmfModel(W = w2, H = h2)
    check.model.dimnames(a, list(NULL, NULL, dn), msg("rownames of input basis are enforced"))
    expect_identical(w2, basis(a), info = msg("basis unchanged"))
    expect_identical(h2[nrow(h):1, ], coef(a), info = msg("coef entries are reordered"))
    msg <- function(...) paste("Basis names from input basis matrix are NOT used to order the components if argument order.basis=FALSE - ", 
        ...)
    w2 <- w
    colnames(w2) <- dn
    h2 <- h
    rownames(h2) <- rev(dn)
    expect_warning(a <- nmfModel(W = w2, H = h2, order.basis = FALSE), "The rownames of the mixture matrix were set")
    check.model.dimnames(a, list(NULL, NULL, dn), msg("rownames of input basis are enforced"))
    expect_identical(w2, basis(a), info = msg("basis unchanged"))
    expect_equivalent(h2, coef(a), info = msg("coef entries are not ordered"))
    msg <- function(...) paste("Basis names from input basis matrix are forced onto to the coef matrix - ", 
        ...)
    w2 <- w
    colnames(w2) <- dn
    h2 <- h
    rownames(h2) <- paste(letters[1:nrow(h)], 2)
    expect_warning(a <- nmfModel(W = w2, H = h2), "The rownames of the mixture matrix were set")
    check.model.dimnames(a, list(NULL, NULL, dn), msg("rownames of input basis are enforced"))
    expect_identical(w2, basis(a), info = msg("basis is unchanged"))
    expect_equivalent(h2, coef(a), info = msg("coef entries are correct"))
})

test_that("test.predict", {
    n <- 100
    m <- 20
    r <- 3
    .msg <- NULL
    mess <- function(...) paste(.msg, ":", ...)
    .msg <- "Artificial matrix"
    V <- matrix(c(c(1, rep(0, n - 1)), c(0, 1, rep(0, n - 2)), 
        rep(c(0, 0, 1, rep(0, n - 3)), m - 2)), n, m)
    res <- .predict.nmf(t(V))
    expect_equal(as.factor(c(1, 2, rep(3, m - 2))), res, info = mess("Return known clusters"))
    .msg <- "Random matrix"
    set.seed(.TestSeed)
    V <- matrix(sapply(sample(n, size = m), function(i) {
        x <- rep(0, n)
        x[i] <- 1
        x
    }), n, m)
    res <- .predict.nmf(V)
    expect_true(is.factor(res), info = mess("Result is a factor"))
    expect_true(length(res) == nrow(V), info = mess("Result is the right size"))
    expect_true(nlevels(res) == ncol(V), info = mess("Result has the right number of levels"))
    .msg <- "NMF model"
    a <- nmfModel(r, V)
    set.seed(.TestSeed)
    a <- rnmf(a)
    res <- predict(a, "samples")
    expect_true(is.factor(res), info = mess("Result is a factor"))
    expect_true(length(res) == ncol(a), info = mess("Result has right size"))
    expect_true(nlevels(res) == nbasis(a), info = mess("Result has right number of levels"))
    .msg <- "Factor"
    res <- predict(a, "features")
    expect_true(is.factor(res), info = mess("Result is a factor"))
    expect_true(length(res) == nrow(a), info = mess("Result has right size"))
    expect_true(nlevels(res) == nbasis(a), info = mess("Result has right number of levels"))
})

test_that("test.purity", {
    checkBounds <- function(x) {
        expect_true(x >= 0, "greater than 0")
        expect_true(x <= 1, "lower than 1")
    }
    x <- as.factor(c(rep(1, 5), rep(2, 10), rep(3, 15)))
    p <- purity(x, x)
    checkBounds(p)
    expect_equal(1, p)
    set.seed(.TestSeed)
    p <- purity(as.factor(sample(x)), x)
    checkBounds(p)
})

test_that("test.silhouette", {
    x <- rmatrix(20, 15)
    cl <- gl(3, 5)
    rownames(x) <- letters[seq(nrow(x))]
    colnames(x) <- LETTERS[seq(ncol(x))]
    si <- NMF:::bigsilhouette(x, cl)
    expect_equal(colnames(x), rownames(si))
    base_si <- silhouette(as.integer(cl), dmatrix = 1 - cor(x))
    expect_identical(rownames(base_si), NULL)
    rownames(base_si) <- colnames(x)
    attr(base_si, "call") <- NULL
    attr(si, "call") <- NULL
    expect_equal(base_si, si)
    
})

test_that("test.sparseness", {
    checkBounds <- function(s) {
        expect_true(s >= 0, "greater than 0")
        expect_true(s <= 1, "lower than 1")
    }
    x <- rep(1, 100)
    s <- sparseness(x)
    checkBounds(s)
    expect_equal(0, s, info = "should be 0")
    x <- c(1, rep(0, 100))
    s <- sparseness(x)
    checkBounds(s)
    expect_equal(1, s)
    n <- 100
    m <- 20
    r <- 3
    set.seed(.TestSeed)
    V <- matrix(rnorm(n * m), n, m)
    s <- sparseness(V)
    checkBounds(s)
    a <- nmfModel(r, V)
    set.seed(.TestSeed)
    a <- rnmf(a, V)
    s <- sparseness(a)
    expect_true(length(s) == 2, info = "Method 'sparseness' returns a 2-length vector")
    checkBounds(s[1])
    checkBounds(s[2])
})

test_that("test.subset", {
    n <- 30
    r <- 5
    p <- 20
    a <- nmfModel(r, n, p)
    a <- rnmf(a)
    expect_true(identical(a[], a), info = "subset [] is OK (identical)")
    expect_true(identical(a[, ], a), info = "subset [,] is OK (identical)")
    expect_true(identical(a[TRUE, ], a), info = "subset [TRUE,] is OK (identical)")
    expect_true(identical(a[, TRUE], a), info = "subset [,TRUE] is OK (identical)")
    expect_true(identical(a[TRUE, TRUE], a), info = "subset [TRUE,TRUE] is OK (identical)")
    expect_true(identical(a[, , TRUE], a), info = "subset [,,TRUE] is OK (identical)")
    expect_true(identical(a[TRUE, TRUE, TRUE], a), info = "subset [TRUE,TRUE,TRUE] is OK (identical)")
    expect_true(identical(a[NULL, ], a[0, ]), info = "subset [NULL,] is OK")
    expect_true(identical(a[, NULL], a[, 0]), info = "subset [,NULL] is OK")
    expect_true(identical(a[NULL, NULL], a[0, 0]), info = "subset [NULL,NULL] is OK")
    expect_true(identical(a[NULL, NULL, NULL], a[0, 0, 0]), info = "subset [NULL,NULL,NULL] is OK")
    expect_error(a[, , ], info = "Error when called with [,,]")
    expect_equal(c(1, p, r), dim(a[1, ]), info = "subset 1 feature is OK (dim)")
    expect_equal(basis(a)[5, , drop = FALSE], basis(a[5, ]), 
        info = "subset 1 feature is OK (basis)")
    expect_equal(coef(a), coef(a[5, ]), info = "subset 1 feature is OK (coef)")
    expect_equal(basis(a)[5, , drop = TRUE], a[5, , drop = TRUE], 
        info = "subset 1 feature dropping is OK (return basis)")
    expect_equal(c(10, p, r), dim(a[1:10, ]), info = "subset more than 1 feature is OK (dim)")
    expect_equal(basis(a)[1:10, , drop = FALSE], basis(a[1:10, 
        ]), info = "subset more than 1 feature is OK (basis)")
    expect_equal(coef(a), coef(a[1:10, ]), info = "subset more than 1 feature is OK (coef)")
    expect_equal(basis(a)[1:10, , drop = TRUE], a[1:10, , drop = TRUE], 
        info = "subset more than 1 feature dropping is OK (return basis)")
    expect_equal(c(n, 1, r), dim(a[, 1]), info = "subset 1 sample is OK (dim)")
    expect_equal(coef(a)[, 5, drop = FALSE], coef(a[, 5]), info = "subset 1 sample is OK (coef)")
    expect_equal(basis(a), basis(a[, 1]), info = "subset 1 sample is OK (basis)")
    expect_equal(coef(a)[, 5, drop = TRUE], a[, 5, drop = TRUE], 
        info = "subset 1 sample dropping is OK (return coef)")
    expect_equal(c(n, 10, r), dim(a[, 1:10]), info = "subset more then 1 sample is OK (dim)")
    expect_equal(coef(a)[, 1:10, drop = FALSE], coef(a[, 1:10]), 
        info = "subset more than 1 sample is OK (coef)")
    expect_equal(basis(a), basis(a[, 1:10]), info = "subset more than 1 sample is OK (basis)")
    expect_equal(coef(a)[, 1:10, drop = TRUE], a[, 1:10, drop = TRUE], 
        info = "subset more than 1 sample dropping is OK (return coef)")
    expect_equal(c(n, p, 1), dim(a[, , 1]), info = "subset 1 basis is OK (dim)")
    expect_equal(coef(a)[3, , drop = FALSE], coef(a[, , 3]), 
        info = "subset 1 basis is OK (coef)")
    expect_equal(basis(a)[, 3, drop = FALSE], basis(a[, , 3]), 
        info = "subset 1 basis is OK (basis)")
    expect_true(identical(a[, , 3, drop = TRUE], a[, , 3, drop = TRUE]), 
        info = "subset 1 basis dropping is OK (do nothing)")
    expect_equal(c(n, p, 3), dim(a[, , 2:4]), info = "subset more than 1 basis is OK (dim)")
    expect_equal(coef(a)[2:4, , drop = FALSE], coef(a[, , 2:4]), 
        info = "subset more than 1 basis is OK (coef)")
    expect_equal(basis(a)[, 2:4, drop = FALSE], basis(a[, , 2:4]), 
        info = "subset more than 1 basis is OK (basis)")
    expect_true(identical(a[, , 2:4, drop = TRUE], a[, , 2:4, 
        drop = FALSE]), info = "subset more than 1 basis dropping is OK (do nothing)")
    expect_equal(c(n, p, 0), dim(a[NULL]), info = "subset basis NULL is OK (dim)")
    expect_identical(basis(a)[, 2], a[2], info = "subset with single index + drop missing returns single basis as vector")
    expect_identical(a[, , 2], a[2, drop = FALSE], info = "subset with single index with drop=FALSE returns the complete NMF object")
    expect_identical(basis(a)[, 2, drop = TRUE], a[2, drop = TRUE], 
        info = "subset with single index + drop=TRUE returns single basis as vector")
    expect_true(is.nmf(a[2:3]), info = "subset with single vector index returns NMF object")
    expect_equal(basis(a)[, 2:3], basis(a[2:3]), info = "subset with single vector index returns subset of NMF object")
    expect_equal(basis(a)[, 2:3, drop = TRUE], a[2:3, drop = TRUE], 
        info = "subset with single vector index + dropping returns correct matrix if length > 1")
    expect_identical(a[, , 2:3], a[2:3, drop = FALSE], info = "subset with single vector index + NOT dropping returns correct matrix if length > 1")
    expect_equal(c(1, 1, r), dim(a[1, 1]), info = "subset 1 feature x 1 sample is OK (dim)")
    expect_equal(coef(a)[, 5, drop = FALSE], coef(a[3, 5]), info = "subset 1 feature x 1 sample is OK (coef)")
    expect_equal(basis(a)[3, , drop = FALSE], basis(a[3, 5]), 
        info = "subset 1 feature x 1 sample is OK (basis)")
    expect_equal(c(10, 7, r), dim(a[10:19, 5:11]), info = "subset more than 1 feature x sample is OK (dim)")
    expect_equal(coef(a)[, 5:11], coef(a[10:19, 5:11]), info = "subset more than 1 feature x sample is OK (coef)")
    expect_equal(basis(a)[10:19, ], basis(a[10:19, 5:11]), info = "subset more than 1 feature x sample is OK (basis)")
})

test_that("test.syntheticNMF", {
    n <- 100
    m <- 20
    r <- 3
    set.seed(.TestSeed)
    d <- syntheticNMF(n, r, m)
    checkPlot(basicHM(d), "Synthetic data plain")
    set.seed(.TestSeed)
    n.offset <- 15
    o <- c(rep(1, n.offset), rep(0, n - n.offset))
    d <- syntheticNMF(n, r, m, offset = o)
    expect_error(syntheticNMF(n, r, m, offset = o[-1]))
    checkPlot(basicHM(d), "Synthetic data offset")
    set.seed(.TestSeed)
    d <- syntheticNMF(n, r, m, noise = TRUE)
    checkPlot(basicHM(d), "Synthetic data with noise")
    set.seed(.TestSeed)
    d <- syntheticNMF(n, r, m, offset = o, noise = TRUE)
    checkPlot(basicHM(d), "Synthetic data with offset and noise")
})

