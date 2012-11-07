##############################################
# Tests for R package NMF:
#
# - prepare data to be used in Matlab (Octave)
# - read data from tests performed in Matlab (Octave)
#
##############################################



##########################
# Utility functions
##########################
source('octave.R')

compare.results <- function(reffile, object, round=12){
	
	# get results from NMF object
	W <- metagenes(object); H <- metaprofiles(object)
	
	# read reference result from file
	ref <- read.octave(reffile)
	# check validity of reference results
	Wref <- ref$W; Href <- ref$H
	if( !is.matrix(Wref) && dim(Wref) == dim(W) ) 
		stop('Invalid W loaded from results')
	if( !is.matrix(Href) && dim(Href) == dim(H) ) 
		stop('Invalid H loaded from results')
	
	# compute comparison measures
	Wdiff = sum( (W - Wref)^2 )
	Hdiff = sum( (H - Href)^2 )
	res <- c( Wdiff=Wdiff, Wdif.rel = Wdiff/sum(W^2)
		, Hdiff=Hdiff, Hdiff.rel = Hdiff/sum(H^2)
		, Timediff = round(runtime(object)['user.self']/ref$elapsed, 2)
	)
	
	if( is.numeric(round) ) res <- round(res, round)
	return(res)
}

generate.data <- function(n=5000, r=3, p=38, nsample=50){

	# create data directory
	datadir <- paste('data',n,r,p, sep='.')
	dir.create(datadir)
	
	# set random seed
	set.seed(654321)

	# Generate random target matrix
	save.ascii(matrix(runif(n*p), n, p), file.path(datadir,'target.txt'))
	
	# set random seed
	set.seed(123456)
	
	lapply(seq(nsample), function(i){
		# Generate random W.ini matrix
		save.ascii(matrix(runif(n*r), n, r), file.path(datadir, paste('W',i,'txt', sep='.')))
		
		# Generate random H.ini matrix
		save.ascii(matrix(runif(r*p), r, p), file.path(datadir, paste('H', i,'txt', sep='.')))
	})
	
	invisible()
}

read.test.data <- function(file='data.test.oct', n=0, dir='.', source=c('octave', 'R')){

	source <- match.arg(source)
	if( source == 'octave'){
		res <- read.octave(file)
	}else{
		res <- list()
		res$target <- as.matrix(read.table(file.path(dir, 'target.txt')))
		
		wfile <- if( n > 0 ) paste('W', n, 'txt', sep='.') else 'W.txt'
		res$W0 <- as.matrix(read.table(file.path(dir, wfile)))
		hfile <- if( n > 0 ) paste('H', n, 'txt', sep='.') else 'H.txt'
		res$H0 <- as.matrix(read.table(file.path(dir, hfile)))
	}
	
	# return loaded data
	return( res )
}

mean.time <- function(target, r, method, ...){
	mean <- 0
	sapply(seq(10), function(i, ...){
		res <- nmf(target, r, method, ...)
		mean <<- mean + runtime(res)
	}
	, ...
	)
	mean / 10
}


algorithms.methods <- list()
algorithms.methods[['brunet']] <- 'brunet'
algorithms.methods[['snmf/r']] <- 'snmfr'
algorithms.methods[['snmf/l']] <- 'snmfl'


run.test.algorithm <- function(method){
	
	# load the test data
	test.data <- read.test.data()
	
	# create new seed
	start <- nmfModel(W=test.data$W0, H=test.data$H0)
	
	# run NMF algorithms on test data
	res <- nmf(test.data$target, nbasis(start), method, seed=start, verbose=TRUE)
	
	# load reference results
	reffile <- paste('ref', algorithms.methods[[method]], 'oct', sep='.')
	# compare results
	compare.results(reffile, res, round=FALSE)
}

test.cophenetic <- function(){
	
	res <- nmfEstimateRank(esGolub, seq(2,6), 'brunet', nrun=50, conf=TRUE, seed=123456, verbose=T)
	save(res, file='esgolub.estimate2.RData')
	
}

compute.var <- function(res){
	sapply(res$bootstrap.measure, function(boot){
				res <- apply(boot, 1, function(row){
							c(mean=mean(row), sd=sd(row))
						})
				res
			}, simplify=F)
}

if( FALSE ){

	# generate random data
	generate.data()
	
	# load the generated data
	input <- read.data()
	
	# run Brunet algorithm on data
	source('nmf-brunet.r')
	s <- nmfModel(W=input$w, H=input$h)
	Rprof()
	res.brunet <- nmf(input$target, 3, 'brunet', seed=s, verbose=TRUE)
	Rprof(NULL)
	res.brunet <- nmf(input$target, 3, wrap.brunet, seed=s, verbose=TRUE)
	system.time( res.2 <- nmf.brunet(input$target, 3, input$w, input$h, verbose=TRUE) )

	mean.time(input$target, 3, 'brunet', seed=s)
	# compare results
	w.matlab <- read.table('brunet.W.txt')
	
	
	# tests consensus on simple matrix
	cons <- matrix(0, 38,38)
	a <- lapply(seq(50), function(i){
		test <- read.test.data(dir='data.5000.3.38', source='R', n=i)
		start <- nmfModel(W=test$W0, H=test$H0)
		res <- nmf(esGolub, 3, 'brunet', seed=start, verbose=T)
		cons <<- cons + connectivity(res)
		invisible()
	})
}
