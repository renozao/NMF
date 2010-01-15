#'
#' SNMF/R  
#'
#' Author: Hyunsoo Kim and Haesun Park, Georgia Insitute of Technology
#'
#' Reference: 
#'
#'   Sparse Non-negative Matrix Factorizations via Alternating 
#'   Non-negativity-constrained Least Squares for Microarray Data Analysis
#'   Hyunsoo Kim and Haesun Park, Bioinformatics, 2007, to appear.
#'
#' This software requires fcnnls.m, which can be obtained from 
#' M. H. Van Benthem and M. R. Keenan, J. Chemometrics 2004; 18: 441-450
#'
#' NMF: min_{W,H} (1/2) || A - WH ||_F^2 s.t. W>=0, H>=0 
#' SNMF/R: NMF with additional sparsity constraints on H
#'
#'   min_{W,H} (1/2) (|| A - WH ||_F^2 + eta ||W||_F^2 
#'                + beta (sum_(j=1)^n ||H(:,j)||_1^2))
#'                s.t. W>=0, H>=0 
#'
#' A: m x n data matrix (m: features, n: data points)
#' W: m x k basis matrix
#' H: k x n coefficient matrix
#'
#' function [W,H,i]=nmfsh_comb(A,k,param,verbose,bi_conv,eps_conv)
#'
#' input parameters:
#'   A: m x n data matrix (m: features, n: data points)
#'   k: desired positive integer k
#'   param=[eta beta]:  
#'      eta (for supressing ||W||_F)
#'         if eta < 0, software uses maxmum value in A as eta. 
#'      beta (for sparsity control)
#'         Larger beta generates higher sparseness on H.
#'         Too large beta is not recommended. 
#'   verbos: verbose = 0 for silence mode, otherwise print output
#'   eps_conv: KKT convergence test (default eps_conv = 1e-4)
#'   bi_conv=[wminchange iconv] biclustering convergence test 
#'        wminchange: the minimal allowance of the change of 
#'        row-clusters  (default wminchange=0)
#'        iconv: decide convergence if row-clusters (within wminchange)
#'        and column-clusters have not changed for iconv convergence 
#'        checks. (default iconv=10)
#'
#' output:
#'   W: m x k basis matrix
#'   H: k x n coefficient matrix
#'   i: the number of iterations
#'
#' sample usage:
#'  [W,H]=nmfsh_comb(amlall,3,[-1 0.01],1);
#'  [W,H]=nmfsh_comb(amlall,3,[-1 0.01],1,[3 10]); 
#'     -- in the convergence check, the change of row-clusters to
#'        at most three rows is allowed.
#'
#' @include fcnnls.R
#' 
#function [W,H,i] 
#nmfsh_comb <- function(A, k, param, verbose=FALSE, bi_conv=c(0, 10), eps_conv=1e-4, version=c('R', 'L')){
.nmfsh_comb <- function(A, k, nmf.fit, eta=-1, beta=0.01, bi_conv=c(0, 10), eps_conv=1e-4, version=c('R', 'L')){

	# retrieve verbosity option
	verbose <- verbose(nmf.fit)
	
	# depending on the version: 
	# in version L: A is transposed while W and H are swapped and transposed
	version <- match.arg(version)
	if( version == 'L' ) A <- t(A) 
	#if( missing(param) ) param <- c(-1, 0.01)
	
	maxiter = 20000; # maximum number of iterations

	m = nrow(A); n = ncol(A); erravg1 = numeric();
	
	#eta=param[1]; beta=param[2]; 
	maxA=max(A); if ( eta<0 ) eta=maxA;
	eta2=eta^2;
	wminchange=bi_conv[1]; iconv=bi_conv[2];
	if ( verbose )
		cat(sprintf("SNMF/%s k=%d eta=%.4e beta (for sparse H)=%.4e wminchange=%d iconv=%d\n",
				version, k,eta,beta,wminchange,iconv));

	idxWold=rep(0, m); idxHold=rep(0, n); inc=0;	
	
	# initialize random W if no starting point is given
	track <- FALSE
	if( missing(nmf.fit) ){
		message('Init W internally (random)')
		W=matrix(runif(m*k), m,k);	
	} else {
		# seed the method (depends on the version to run)
		start <- if( version == 'R' ) basis(nmf.fit) else t(coef(nmf.fit))
		# check compatibility of the starting point with the target matrix
		if( any(dim(start) != c(m,k)) ) 
			stop("Invalid initialization: incompatible dimensions [expected: ", paste(c(m,k), collapse=' x '),", got: ", paste(dim(start), collapse=' x '), " ]")	
		# use the supplied starting point
		W <- start
		
		# use runtime options for tracking
		track <- run.options(nmf.fit, 'error.track')
		track.interval <- run.options(nmf.fit, 'track.interval')
	}
		
	W= apply(W, 2, function(x) x / sqrt(sum(x^2)) );  # normalize columns of W	

	I_k=diag(eta, k); betavec=rep(sqrt(beta), k); nrestart=0;	
	for ( i in 1:maxiter ){

		# min_h ||[[W; 1 ... 1]*H  - [A; 0 ... 0]||, s.t. H>=0, for given A and W.
	  	res = .fcnnls(rbind(W, betavec), rbind(A, rep(0, n)));	  	  	
		H = res[[1]]

		if ( any(rowSums(H)==0) ){ 	  
			cat(sprintf("iter%d: 0 row in H eta=%.4e restart!\n",i,eta));
			nrestart=nrestart+1;
			if ( nrestart >= 10 ){
				print('[*Warning*] too many restarts due to too big beta value...\n');
				break;
			}
						
			# re-initialize random W
			idxWold=rep(0, m); idxHold=rep(0, n); inc=0;
			W=matrix(runif(m*k), m,k);
			W= apply(W, 2, function(x) x / sqrt(sum(x^2)) );  # normalize columns of W	
			next;
		}
		
		# min_w ||[H'; I_k]*W' - [A'; 0]||, s.t. W>=0, for given A and H. 
		res = .fcnnls(rbind(t(H), I_k), rbind(t(A), matrix(0, k,m))); 
		Wt = res[[1]]
		W= t(Wt);		

		# track the error (not computed unless tracking option is enabled in nmf.fit)
		nmf.fit <- trackError(nmf.fit, .snmf.objective(A, W, H, eta, beta), i)
		
		# test convergence every 5 iterations
		if ( (i %% 5==0)  || (i==1) ){
			idxW = apply(W, 1, function(x) which.max(x)) 			
			idxH = apply(H, 2, function(x) which.max(x))			
			changedW=sum(idxW != idxWold); changedH=sum(idxH != idxHold);
			if ( (changedW<=wminchange) && (changedH==0) ) inc=inc+1
			else inc=0

			resmat=pmin(H, (t(W) %*% W) %*% H - t(W) %*% A + matrix(beta, k , k) %*% H); resvec=as.numeric(resmat);
			resmat=pmin(W, W %*% ( H %*% t(H)) - A %*% t(H) + eta2 * W); resvec=c(resvec, as.numeric(resmat));
			conv=sum(abs(resvec)); #L1-norm      
			convnum=sum(abs(resvec)>0);
			erravg=conv/convnum;
			if ( i==1 ){
				erravg1=erravg;
				if( verbose ) cat("\tIter\tInc\tchW\tchH\t---\terravg1\terravg\terravg/erravg1\n")
			}

			if ( verbose || (i %% 1000==0) ) # prints number of changing elements 
			 cat(sprintf("\t%d\t%d\t%d\t%d\t---\terravg1: %.4e\terravg: %.4e\terravg/erravg1: %.4e\n",
				 i,inc,changedW,changedH,erravg1,erravg,erravg/erravg1));
			
			if ( (inc>=iconv) && (erravg<=eps_conv*erravg1) ) break;
			idxWold=idxW; idxHold=idxH; 
		}

	}
	
	# force to compute last error if not already done
	nmf.fit <- trackError(nmf.fit, .snmf.objective(A, W, H, eta, beta), i, force=TRUE)	

	# transpose and reswap the roles
	if( version == 'L' ){
		basis(nmf.fit) <- t(H)
		coef(nmf.fit) <- t(W)
	}
	else{ basis(nmf.fit) <- W; coef(nmf.fit) <- H}
	# set number of iterations
	nmf.fit@extra$iteration=i
	
	return(nmf.fit)	
}


#' Computes the objective value for the SNMF algorithm
.snmf.objective <- function(target, w, h, eta, beta){
	
	1/2 * ( sum( (target - (w %*% h))^2 ) 
				+ eta * sum(w^2) 
				+ beta * sum( colSums( h )^2 )
				)
}

#' Wrapper function to use the SNMF/R algorithm with the NMF package.
#'
.snmf <- function(target, seed, ...){	
	
	# retrieve the version of SNMF algorithm from its name: 
	# it is defined by the last letter in the method's name (in upper case)
	name <- algorithm(seed)
	version <- toupper(substr(name, nchar(name), nchar(name)))
	
	# perform factorization using Kim and Park's algorithm 
	sol <- .nmfsh_comb(target, nbasis(seed), nmf.fit=seed, version=version, ...)
		
	# return solution
	return(sol)
}

# internal function to register the methods
.load.algorithm.snmf <- function(){
	nmfRegisterAlgorithm(.snmf, 'snmf/r', objective='euclidean', overwrite=TRUE)
	nmfRegisterAlgorithm(.snmf, 'snmf/l', objective='euclidean', overwrite=TRUE)
}
