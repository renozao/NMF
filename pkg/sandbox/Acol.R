###% Seeding method for Nonnegative Matrix Factorization: Random Acol
###%
###% @author Renaud Gaujoux
###% @created 21 Jul 2009
###% 
###% @include seed.R

###% Seeding method for Nonnegative Matrix Factorization: Random Acol
###%
###% Suppose a NMF: $$V = WH$$
###% Random Acol forms an initialization of each column of the basis matrix W by averaging \code{p} random columns of V. 
###% In some case as of document clustering, it makes more sense to build basis vectors from the given data, 
###% the sparse document vectors themselves, than to form completely dense random basis vectors, as random initialization does. 
###% Random Acol initialization is very inexpensive, and lies between random initialization and centroid initialization
###% in terms of performance. 
###%
###% @references Algorithms, Initializations, and Convergence for the Nonnegative Matrix Factorization
###% Russell Albright, James Cox, David Duling, Amy N Langville, Carl D Meyer 
###%
seed.acol <- function(x, k, p){
	
	# check parameters
	if( p > ncol(x) ) stop("The number of columns to average cannot be greater than the number of columns in the target matrix")
	
	# generate each basis as the average of p random columns
	sapply(seq(k), 
		function(i){
			# randomly select the p columns to average
			sel <- sample(seq(ncol(x)), p)
			# return average of the selected columns
			apply(x[,sel], 1, mean)
		}
	)
}

###% Seeding method for Nonnegative Matrix Factorization: Random C
###%
###% Suppose a NMF: $$V = WH$$
###% The Random C initialization is similar to the Random Acol method, except it chooses p columns at random
###% from the longest (according to the given distance) columns of V, which generally means the densest 
###% columns in the context of text matrices since those matrices are so sparse.
###% The idea is that these might be more likely to be the centroid centers.
###%
###% @references Algorithms, Initializations, and Convergence for the Nonnegative Matrix Factorization
###% Russell Albright, James Cox, David Duling, Amy N Langville, Carl D Meyer 
###%
seed.randomC <- function(x, k, n, p, measure){
	
	# check parameters
	if( p > ncol(x) ) stop("The number of columns to average cannot be greater than the number of columns in the target matrix [p > ncol(x)]")
	if( p > n ) stop("The number of columns to average cannot be greater than the number of vectors in the base set [p > n]")
	if( missing(measure) ) measure <- function(x) sum(x^2)
			
	# compute measure for each columns using the provided distance method
	d <- apply(x, 2, measure)
	# select the n longest vectors
	sel.base <- order(d, decreasing=TRUE)[1:n]
	sapply(seq(k), 
		function(i){
			# select p random vectors amongst the n longest vectors to average
			sel <- sample(sel.base, p)
			# return average of the selected columns
			apply(x[,sel], 1, mean)
		}
	)
}
