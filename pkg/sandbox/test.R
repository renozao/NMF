###% Testing script for NMF package
###%
###% @author Renaud Gaujoux
###% @creation 10 April 2009

# load scripts
source('nmf.R')


################################################################
# Test on the Golub dataset
################################################################
# load data
data(esGolub)

# set seed for reproducibility
set.seed(123456)

# set the factorization rank
r <- 3

#####################
# Synthetic data
#####################
syn <- 	syntheticNMF(100, 3, 20)
basicHM(syn)

# with offset
n.offset <- 15
syn <- 	syntheticNMF(100, 3, 20, offset=c(rep(1, n.offset), rep(0, 100-n.offset)))
	

# single NMF run
s <- nmf(esGolub, r) # with default strategy
s
s <- nmf(esGolub, r, 'brunet') # with strategy given by name
s
strat <- nmfStrategy() # build default strategy
strat
# test brunet
strat <- nmfStrategy('brunet') # build given strategy
strat
s <- nmf(esGolub, r, strat) # with strategy given as an object
s
# test lee
strat <- nmfStrategy('lee') # build given strategy
strat
s <- nmf(esGolub, r, strat) # with strategy given as an object
s
# test offset
strat <- nmfStrategy('offset') # build given strategy
strat
s <- nmf(esGolub, r, strat) # with strategy given as an object
s

# multiple NMF runs
n <- 10
run <- nmf(esGolub, r, nrun=n)


# conditional import of package Biobase 
if( "Biobase" %in% rownames(utils::installed.packages()) ){
	#do.call('importClassesFrom', list(Biobase, 'ExpressionSet'))
	importClassesFrom(Biobase, ExpressionSet)	
}


#Package the code into a R-package.
#There will be things to adapt to make it work:
#		
#		* loading of built-in methods (algorithms and seeding methods) should be done in function  @.First.lib@
#* the namespace might create problems with the registry and the package specific options (function @nmf.options@)

###################################
# Tissue heterogeneity
###################################

# load data
load('../package/data/CellType.rda')

esCellType
eSet <- esCellType

# only keep genes with expression that varies some more than a given fold over the mean across the samples
library(genefilter)
fold.from.mean <- function(n=3, fold=3){ function(x){ f <- abs(x/mean(x)); (sum(f  >= fold) >= n)  } }

ff <- filterfun(fold.from.mean(3,3))
genes.ok <- genefilter(eSet, ff) 
sum(genes.ok)
eSet <- subset(eSet, subset= genes.ok)

# create an artificial mixture matrix




plot.NMF.rank <- function(estimates, what=c('all', 'cophenetic', 'rss', 'error') ){
	
	if( what == 'auc' ){
		c.list <- estimates$consensus
		auc <- sapply(c.list, function(C){
					# compute empircal cdf
					cdf <- consensusCDF(C)
					# return area under the curve
					x <- knots(cdf)
					return( diff(x) %*% cdf(x[-length(x)]) )
				})
		
		return( plot(names(c.list), auc) )
	}
	
	measures <- estimates$measures
	what <- match.arg(what, rownames(measures))
	vals <- measures[what,]
	x <- colnames(measures)
	plot(x=x, y=vals, ylab=paste('Quality measure:', what), xlab='Factorization rank', main='NMF quality plot')
	lines(x=x, y=vals)
	
}

###% Compute the value of the cumulative distribution function from a consensus 
###% matrix at a given point or on the ordered values of the matrix
if ( !isGeneric('consensusCDF') ) setGeneric('consensusCDF', function(object, x, ...) standardGeneric('consensusCDF'))
setMethod('consensusCDF', signature(object='matrix', x='missing'),
		function(object, x){
			
			res <- ecdf( as.numeric(object[upper.tri(object)]) )
			return( res )
		}
)

setMethod('consensusCDF', signature(object='matrix', x='numeric'),
		function(object, x){
			
			if( x<0 || x>1 ) stop("Invalid evaluation point: consensus CDF is only defined on [0,1]")
			return( consensusCDF(object)(x) )
			
		}
)

setMethod('consensusCDF', signature(object='NMFSet', x='ANY'),
		function(object, x){
			# return the CDF of the consensus matrix at point x
			return( consensusCDF(consensus(object), x) )
		}
)

if ( !isGeneric('consensusPlot') ) setGeneric('consensusPlot', function(object, ...) standardGeneric('consensusPlot'))
setMethod('consensusPlot', signature(object='matrix'),
		function(object, what=c('all', 'hist', 'ecdf', 'matrix'), ...){
			
		}
)

###############################################
# Re-test C vs. R
###############################################

cmp.c <- function(){
	
	set.seed(123456)
	n <- 3; p <- 3; r <- 2 
	h <- rmatrix(r, p)
	w <- rmatrix(n, r)
	v <- rmatrix(n, p)
	
	eps <- 10^-9	
	.R <- R_std.euclidean.update.h(v, w, h, eps=eps)
	.C <- .Call("euclidean_update_H", v, w, h, eps, FALSE)	
	cat(" - H: ", identical(.R, .C), " : ", all.equal(.R, .C, tol=0), "\n")
	
	.R <- R_std.euclidean.update.w(v, w, h, eps=eps)
	.C <- .Call("euclidean_update_W", v, w, h, eps, NULL, FALSE)
	cat(" - W: ", identical(.R, .C), " : ", all.equal(.R, .C, tol=0), "\n")
}

###################################
# Parallel computations
###################################


library(foreach)

testdoppar <- function(){
registerDoSEQ()
err <- NA 
res <- foreach(i=2:5) %dopar% {
	print(err)
	print(environment())
	print(ls(environment()))
	e <- runif(1)
	if( is.na(err) || e < err)
		err <<- e
	NULL
}
res
}

n <- 50; r <- 3; m <- 20
V <- syntheticNMF(n, r, m, noise=TRUE)

library(foreach)
registerDoSEQ()
res <- foreach(i=2:5) %dopar% {
	f <- nmf(V, i, nrun=4, .opt='vp')
	print(getDoParName())
	f
}

res <- foreach(i=1:5, .combine='as.list') %dopar%{ i }

nw <- 4
ntasks <- rep(5, nw)
set.seed(12345)
res <- foreach(i=1:nw, n=ntasks) %dopar% { nmf(V, 3, 'brunet', nrun=n, .opt='v') }
join(res, .merge=TRUE)


library(doMPI)
cl <- startMPIcluster(verbose=TRUE)
registerDoMPI(cl)
rng <- RNGseq(3, seed=12334)
foreach(i=1:3, s=rng) %dopar%{ message("AA"); rstream.packed(s) <- FALSE; r(s); }
closeCluster(cl)

#####################################

source('../package/R/package.R', chdir=TRUE)
n <- 50; r <- 3; m <- 20
V <- syntheticNMF(n, r, m, noise=TRUE)
res1 <- nmf(V,r, 'R_brunet', nrun=3)
resNew <- new('NMFfitX1', res[[1]],  consensus=res@consensus, totaltime=res@runtime, nrun=res@nrun)
resn <- nmf(V,r, 'R_brunet', nrun=10, .opt='k')
resn <- new('NMFfitXn', res, totaltime=res@runtime)

resM <- nmf(V,r, list('R_lee', 'R_offset', 'R_brunet'), nrun=3)

######################################
# New evar
######################################

library('NMF', lib='../lib')
data(esGolub)
x <- esGolub[1:200,]
res <- nmf(x, 3, seed=123456)
resl <- nmf(x, 3, 'lee', seed=123456)

######################################
# Set seed
######################################

library(doMC)

registerDoMC()
seqrand <- matrix(NA, 5,length(.Random.seed))
set.seed(123); res <- for(n in 1:5){ 	seqrand[n, ] <- .Random.seed; cat(n , ': ', paste(head(.Random.seed), collapse=' '), "\n"); runif(5000*3+38*3); cat(n , ': ', paste(sample(5), collapse=' '), "\n");}


a <- function(seed){
	
	use.seed <- FALSE
	if( !missing(seed) ){
	set.seed(seed)
# - Define a shared memory objects
seed.shared <- shared.big.matrix(1, length(.Random.seed), type='integer', init=NA)
seed.shared[] <- .Random.seed
seed.desc <- bigmemory::describe(seed.shared)				

# - Define a mutex to control the access to the shared memory objects
mut <- rw.mutex()
mut.desc <- bigmemory::describe(mut)
use.seed <- TRUE
}

nrun <- 5
## 2. RUN
foreach(n=1:nrun, .combine=rbind) %dopar% {
		
	res <- NA
	# if only the best fit must be kept then update the shared objects
	if( use.seed ){
		# load shared objects
		seed.shared <- attach.big.matrix(seed.desc)

		##LOCK_MUTEX					
		# retrieve and lock the mutex
		mut <- attach.rw.mutex(mut.desc)	
		rwlock(mut)
		
		# check if the run found a better fit
		seed <- as.integer(seed.shared[])
		#old.seed <- setRNG(kind = "default", seed = seed, normal.kind ="default")
		.Random.seed <<- seed
		res <- .Random.seed		
		cat(n , ': ', paste(head(.Random.seed), collapse=' '), "\n")
		runif(5000*3+38*3);
		cat(n , ': ', paste(sample(5), collapse=' '), "\n")
		
		# update seed
		seed.shared[] <- .Random.seed
						
		# unlock the mutex
		bigmemory::unlock(mut)
		##END_LOCK_MUTEX
	}
	else{
		res <- .Random.seed
		cat(n , ': ', paste(head(.Random.seed), collapse=' '), "\n")
		cat(n , ': ', paste(sample(5), collapse=' '), "\n")		
	}
	
	res
}

}

parrand <- a(123)
all.equal(seqrand, parrand, check.attributes=FALSE)


######################################
# C++ distances
######################################

n <- 5000; p <- 50;

# create random matrices in storage mode double 
x <- matrix(runif(n*p, 1, 100), n, p)
y <- matrix(runif(n*p, 1, 100), n, p)

# create random matrices in storage mode integer 
xi <- x; storage.mode(xi) <- 'integer'
yi <- y; storage.mode(yi) <- 'integer'

# RSS
R_rss <- function(x, y) sum((x-y)^2)
system.time(replicate(1000, R_rss(x,y)))
system.time(replicate(1000, .rss(x,y)))

# KL
R_kl <- function(x, y) sum( ifelse(x==0, y, x * log(x/y) - x + y) );
system.time(replicate(100, R_kl(x,y)))
system.time(replicate(100, .KL(x,y)))

######################################
# New version of bigmemory
######################################

library('bigmemory', lib='../lib')
library('synchronicity', lib='../lib')
source('../package/R/package.R', chd=T)
load('../package/data/esGolub.rda')

res <- nmf(esGolub[1:200,], 3, nrun=4, .opt='v')


######################################
# Parallel RNG
######################################

library(doMC)
library(rstream)
registerDoMC()
#registerDoSEQ()


get.seed <- function(x){	
	state <- .Call("R_RngStreams_GetData", x@stream, PACKAGE = "rstream")
	state
}


setupRNG <- function(i, n, seed){
	
	# check parameters
	if( i <= 0 )
		stop("NMF::setupRNG - invalid value for 'i' [positive value expected]")
	if( n <= 0 )
		stop("NMF::setupRNG - invalid value for 'n' [positive value expected]")
	if( i > n )
		stop("NMF::setupRNG - invalid value for 'i' [must be lower or equal to 'n']")
	
	s <- list(new('rstream.mrg32k3a', seed=seed, force.seed=TRUE))
	if( n > 0 ) 
		s <- c(s, replicate(n-1, new('rstream.mrg32k3a') ))
	
	invisible(rstream.RNG(s[[i]]))
}

# This does NOT ensure reproducibility
set.seed(1)
l <- foreach(i=1:4) %dopar% {
	
	if( i == 1 || i == 2 )
		s <- new('rstream.mrg32k3a', seed=1:6, force.seed=TRUE)
	else 
		s <- new('rstream.mrg32k3a')
	message("After init: ", i, ":", paste(get.seed(s), collpase=', '))	
}

l <- foreach(i=1:4, s=createStream(4, 1:6, TRUE)) %dopar% {
	rstream.packed(s) <- FALSE
	message("Iteration ", i, ": ", r(s))	
}

s <- new('rstream.mrg32k3a', seed=sample(1:999999, 6), force.seed=TRUE)
s <- new('rstream.mrg32k3a', seed=1:6, force.seed=TRUE)
message("Before loop: ", paste(rstream:::.rstream.envir$.rstream.mrg32k3a.DefaultSeed, collapse=", "))
init <- TRUE
l <- foreach(i=1:4) %dopar% {
	
	message("Before ", i, ": ", paste(rstream:::.rstream.envir$.rstream.mrg32k3a.DefaultSeed, collapse=", "))
	if( init ){
		message("Init in ", i)
		s <- new('rstream.mrg32k3a', seed=2:7, force.seed=TRUE)
		init <<- FALSE
	}
	
	if( i > 1 ){
		s <- replicate(i, new('rstream.mrg32k3a'))
		message("In ", i, " : create ", i)
		s <- s[[i]]
	}
	
	message("Iteration ", i, ": ", res <- r(s))
	res
}
message("After loop: ", paste(rstream:::.rstream.envir$.rstream.mrg32k3a.DefaultSeed, collapse=", "))
l
	
#		s <- new('rstream.mrg32k3a')	
#	
#	message("Iteration: ", i)
#	print(s)
#	rstream.sample(s, 2)
	
}

sapply(1:4, function(i){
	
	init <- TRUE
	if( init ){
		s <- new('rstream.mrg32k3a', seed=1:6, force.seed=TRUE)
		init <- FALSE
	}
	else
		s <- new('rstream.mrg32k3a')	
	
	message("Iteration: ", i)
	print(s)
	rstream.sample(s, 2)
	
})


ash <- function (x, ...) 
{
	stopifnot(is.list(x), length(x) == 2)
	n <- length(ord <- unlist(x))
	stopifnot(n == attr(x, "members"))
	n.h <- n - 1L
	labsu <- unlist(labels(x))
	labs <- labsu[sort.list(ord)]
	x <- stats:::.add.dendrInd(x)
	SIMP <- function(d) {
		if (is.leaf(d)) {
			-as.vector(d)
		}
		else {
			j <<- j + 1L
			height[j] <<- attr(d, "height")
			inds[[j]] <<- attr(d, ".indx.")
			attributes(d) <- NULL
			d[] <- lapply(d, SIMP)
			d
		}
	}
	height <- numeric(n.h)
	inds <- vector("list", n.h)
	j <- 0L
	xS <- SIMP(x)
	ii <- sort.list(height)
	merge <- matrix(NA_integer_, 2L, n.h)
	for (k in seq_len(n.h)) {
		if (k < n.h) {
			in.k <- inds[[ii[k]]]
			s <- xS[[in.k]]
		}
		else s <- xS		
		stopifnot(length(s) == 2L, all(vapply(s, is.integer, 
								NA)))
		merge[, k] <- unlist(s)
		if (k < n.h) 
			xS[[in.k]] <- +k
	}
	r <- list(merge = t(merge), height = height[ii], order = ord, 
			labels = labs, call = match.call(), method = NA_character_, 
			dist.method = NA_character_)
	class(r) <- "hclust"
	r
}
ash(m$colDendrogram)

#################################################################
# RCPP STUFFS
#################################################################

static SEXP data = R_NilValue;
SEXP test(SEXP fun){
	
	using namespace Rcpp;
	
	BEGIN_RCPP
	
	if( isNull(fun) )
		return data;
	
	Environment glob = Environment::global_env();
	Function Rfun = glob[as<const char*>(fun)];
	data = Rfun();
	return data;
	END_RCPP
	
}


#################################################################
# NMF UTILS
#################################################################

maxAD <- function(x, y) max(abs(scoef(x) - scoef(y)))

minfit2 <- function(object, y, method=c('adeviance'), what=c('scoef')){
	
	what <- match.arg(what)
	method <- match.arg(method)
	
	m <- which.min(sapply(object, maxAD, y))
	object[[m]]
}

########################################################################################################
# Generic functions to work with BioConductor datasets 
# 
# Author: Renaud Gaujoux
# Date: 02/12/2010
########################################################################################################

library(Biobase)

###% Subsetting ExpressionSet Objects to Only Keep Core Probesets
###% 
###% The function \code{subsetEZ} filters out probes that are \emph{not} mapped 
###% to an Entrez ID.
###% 
###% @param eset an object of class ExpressionSet
###% @param annotation the annotation package, as a character string, that 
###% provides the mapping between probeset IDs and Entrez IDs.
###% @return The filtered ExpressionSet object
###% @author Renaud Gaujoux
###% @rdname bioc
###% @export
subsetEZ <- function(eset, annotation=NULL){
	
	# load annotation package
	if( is.null(annotation) )
		annotation <- annotation(eset)
	
	pkg.db <- annotation
	if( !grepl('\\.db$', pkg.db))
		pkg.db <- paste(pkg.db, '.db', sep='')
	
	library(pkg.db, character.only=TRUE)
	
	# load the ENTREZID
	pkg <- sub('(.*)\\.db$', '\\1', pkg.db)
	ezID <- get(paste(pkg, 'ENTREZID', sep=''))
	
	# return filtered data
	eset[mappedkeys(ezID[featureNames(eset)]), ]
}

###% Filters out Affymetrix control probesets, i.e. the probesets whose ID starts
###% with "AFFX-".
###%  
###% @param eset an ExpressionSet object 
###% @rdname bioc
###% @export
removeAffyControls <- function(eset){
	
	# remove probesets whose id starts with "AFFX-" 
	eset[ !grepl("^AFFX", featureNames(eset)), ]
	
}

###% \code{coreProbes} filters the data to leave only core probesets. For gene 
###% expression data, this means successively applying the functions 
###% \code{\link{removeAffyControls}} and optionally \code{\link{subsetEZ}}.
###%  
###% @param eset an ExpressionSet object 
###% @param ez.only a logical that specify if only probes mapped to an Entrez ID 
###% should be kept. (Default is \code{FALSE})
###% @param ... extra arguments passed to \code{\link{subsetEZ}}.
###% @rdname bioc
###% @export
coreProbes <- function(eset, ez.only=FALSE, ...){
	
	# keep only probes mapped to an Entrez Gene id
	if( ez.only )
		eset <- subsetEZ(eset, ...)	
	
	# remove control probes
	eset <- removeAffyControls(eset)
	
	# return result
	eset
	
}



###########################################
# DEBUGGING PACKAGE NMF
###########################################
pkgname <- "NMF"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library(NMF, lib='../build/lib')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()

rm(list = ls(all.names = TRUE))
RNGkind("default", "default")
kind <- normal.kind <- "default"


pkgname <- "NMF"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library(NMF)

r <- new('rstream.mrg32k3a')
setRNG(r)
library(randtoolbox)
RNGwrap('rstream')
RNGlibs(all=TRUE)
detach('package:randtoolbox', unload=TRUE)
assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()


###########################################
# HEATMAPS
###########################################
source('../pkg/R/package.R', chd=T)
source('../pkg/R/pheatmap.R', chd=T)
set.seed(123)
x <- rmatrix(100,20)
cann <- data.frame(Type=factor(sample(letters[1:3], ncol(x), replace=TRUE))
, Treatment=factor(sample(c('ppp', 'pet', 'sac'), ncol(x), replace=TRUE))
, Findex= runif(ncol(x)))

rann <- data.frame(Type=factor(sample(letters[1:3], nrow(x), replace=TRUE))
		, Treatment=factor(sample(c('ppp', 'pet', 'sac'), nrow(x), replace=TRUE))
		, Findex= runif(nrow(x)))
pdf("heatmaps.pdf")
aheatmap(x)
aheatmap(x, Colv=NA)
aheatmap(x, Rowv=NA)
aheatmap(x, annCol = cann)
aheatmap(x, annRow = rann)
dev.off()


x <- rmatrix(20,10)
res <- nmf(x, 2, seed=123, .opt='v2')
res <- nmf(x, 2, seed=123, nrun=2, .opt='v2-p')
res <- nmf(x, 2, seed=123, nrun=2, .opt='v2p')