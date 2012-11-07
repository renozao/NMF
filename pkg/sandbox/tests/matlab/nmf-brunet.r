#' Computes the Nonegative Matrix Factorization of a matrix.
#'
#' This software and its documentation are copyright 2004 by the
#' Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
#' This software is supplied without any warranty or guaranteed support whatsoever. 
#' Neither the Broad Institute nor MIT can not be responsible for its use, misuse, 
#' or functionality. 
#'
#' @param v N (genes) x M (samples) original matrix to be factorized.
#'           Numerical data only. 
#'           Must be non negative. 
#'           Not all entries in a row can be 0. If so, add a small constant to the 
#'           matrix, eg.v+0.01*min(min(v)),and restart.
#'           
#' @param r the number of desired factors (i.e. the rank of factorization)
#' @param verbose prints iteration count and changes in connectivity matrix elements unless verbose is 0 
#' @return The two factors of the factorization:
#' \item{w }{N x r NMF factor}
#' \item{h }{r x M NMF factor}
#'
#' @note NMF iterations stop when connectivity matrix has not changed 
#'        for 10*stopconv interations. This is experimental and can be
#'        adjusted.
#' @author Jean-Philippe Brunet \email{brunet@@broad.mit.edu}
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @references
#' NMF divergence update equations :
#' Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix 
#' Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.
nmf.brunet <- function(v, r, w, h, theta=0, verbose=FALSE, useMatrix=FALSE, maxIter=2000){
	
	if( FALSE && useMatrix ){ #overload the matrix creating function: use Matrix instead
		stopifnot( require(Matrix) )
		matrix.old <- matrix
		matrix <- Matrix
		v = Matrix(v);
	}	
	
	# test for negative values in v
	if( min(v) < 0 ) stop('matrix entries can not be negative');
	if( min(rowSums(v)) == 0) stop('not all entries in a row can be zero');
	
	n = nrow(v)
	m = ncol(v);
	stopconv = 40;      # stopping criterion (can be adjusted)
	niter = maxIter;     # maximum number of iterations (can be adjusted)

	cons = matrix(0,m,m);
	consold = cons;
	inc = 0;
	j = 0;
	#MATLAB precision
	eps = 2^-52

	#
	# initialize random w and h
	#
	if( missing(w) ) w = matrix(runif(n*r), n, r);
	if( missing(h) ) h = matrix(runif(r*m), r, m); 

	if( theta > 0 ){
		S <- diag(1-theta, r) + theta / r
		message("Using nsNMF")
	}

	for( i in 1:niter ){

		# divergence-reducing NMF iterations
		if( theta > 0 ){ bkp <- w; w <- w %*% S}
		x1 <- colSums(w) #matrix(rep(colSums(w), ncol(h)), nrow=nrow(h), byrow=F);
		h <- h * (t(w) %*% (v / (w%*%h))) / x1;
		
		if( theta > 0 ){ w <- bkp; bkp<-h; h <- S %*% h}
		x2 <- matrix(rep(rowSums(h), nrow(w)), ncol=ncol(w), byrow=T); 
		w <- w * ((v / (w %*% h)) %*% t(h)) / x2;
		if( theta > 0 ){ h <- bkp; w <- apply(w, 2, function(x) x/sum(x))}
		
		# test convergence every 10 iterations
		if( i %% 10 == 0 ){  
			#cat("Test: ");
			j=j+1;

			# adjust small values to avoid undeflow
			if( theta == 0 ){
			h[h<eps] <- eps;
			w[w<eps] = eps;
			}
			# construct connectivity matrix
			index = apply(h, 2, function(x) which( x == max(x) ) )
			#[y,index]=max(h,[],1);   #find largest factor
			#mat1=repmat(index,m,1);  # spread index down
			#mat2=repmat(index',1,m); # spread index right
			cons = outer(index, index, function(x,y) ifelse(x==y, 1,0));

			if( all(cons==consold) ){ # connectivity matrix has not changed
				inc=inc+1; 			  #accumulate count 			                     
			}else inc=0;                         # else restart count
			
			# prints number of changing elements 
			#if( verbose ) cat( sprintf('\t%d\t%d\t%d\n',i,inc,sum(cons != consold)) )
			if( verbose ) cat( sprintf('%d ', sum(cons != consold)) ) 
			
			# assume convergence is connectivity stops changing 
			if(inc>stopconv) break;			

			consold=cons;
		}
	}	
	
	#return result
	res = list(w=w, h=h)
	class(res) <- 'NMF'
	invisible(res)
}

wrap.brunet <- function(V, start, ...){

	res <- nmf.brunet(V, nbasis(start), metagenes(start), metaprofiles(start), verbose=verbose(start), ...)
	metagenes(start) <- res$w
	metaprofiles(start) <- res$h
	start
}

