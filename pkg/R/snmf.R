#' M. H. Van Benthem and M. R. Keenan, J. Chemometrics 2004; 18: 441-450
#'
#' Given A and C this algorithm solves for the optimal 
#' K in a least squares sense, using that
#'      A = C*K 
#' in the problem
#'      min ||A-C*K||, s.t. K>=0, for given A and C.
#'
#'
#' @param C the matrix of coefficients
#' @param A the target matrix of observations
#'
#' @return [K, Pset]
#'
.fcnnls <- function(C, A){
# NNLS using normal equations and the fast combinatorial strategy
	#
	# I/O: [K, Pset] = fcnnls(C, A);
	# K = fcnnls(C, A);
	#
	# C is the nObs x lVar coefficient matrix
	# A is the nObs x pRHS matrix of observations
	# K is the lVar x pRHS solution matrix
	# Pset is the lVar x pRHS passive set logical array
	#
	# M. H. Van Benthem and M. R. Keenan
	# Sandia National Laboratories
	#
	# Pset: set of passive sets, one for each column
	# Fset: set of column indices for solutions that have not yet converged
	# Hset: set of column indices for currently infeasible solutions
	# Jset: working set of column indices for currently optimal solutions
	#
	# Check the input arguments for consistency and initializeerror(nargchk(2,2,nargin))
	nObs = nrow(C); lVar = ncol(C);
	if ( nrow(A)!= nObs ) stop('C and A have imcompatible sizes')
	pRHS = ncol(A);
	W = matrix(0, lVar, pRHS);
	iter=0; maxiter=3*lVar;
	# Precompute parts of pseudoinverse
	#CtC = t(C)%*%C; CtA = t(C)%*%A;
	CtC = crossprod(C); CtA = crossprod(C,A);
	
	# Obtain the initial feasible solution and corresponding passive set
	K = .cssls(CtC, CtA);
	Pset = K > 0;
	K[!Pset] = 0;
	D = K;
	# which columns of Pset do not have all entries TRUE?
	Fset = which( colSums(Pset) != lVar );
	#V+# Active set algorithm for NNLS main loop
	oitr=0; # HKim
	while ( length(Fset)>0 ) {
		
		oitr=oitr+1; if ( oitr > 5 ) cat(sprintf("%d ",oitr));# HKim
		
		#Vc# Solve for the passive variables (uses subroutine below)				
		K[,Fset] = .cssls(CtC, CtA[,Fset, drop=FALSE], Pset[,Fset, drop=FALSE]);
		
		# Find any infeasible solutions
		# subset Fset on the columns that have at least one negative entry
		Hset = Fset[ colSums(K[,Fset, drop=FALSE] < 0) > 0 ];
		#V+# Make infeasible solutions feasible (standard NNLS inner loop)
		if ( length(Hset)>0 ){
			nHset = length(Hset);
			alpha = matrix(0, lVar, nHset);
			while ( nHset>0  && (iter < maxiter) ){
				iter = iter + 1; 
				alpha[,1:nHset] = Inf;
				#Vc# Find indices of negative variables in passive set
				ij = which( Pset[,Hset, drop=FALSE] & (K[,Hset, drop=FALSE] < 0) , arr.ind=TRUE);			
				i = ij[,1]; j = ij[,2]
				if ( length(i)==0 ) break;			
				hIdx = (j - 1) * lVar + i; # convert array indices to indexes relative to a lVar x nHset matrix
				negIdx = (Hset[j] - 1) * lVar + i; # convert array indices to index relative to the matrix K (i.e. same row index but col index is stored in Hset)
				
				alpha[hIdx] = D[negIdx] / (D[negIdx] - K[negIdx]);				
				alpha.inf <- alpha[,1:nHset, drop=FALSE]
				minIdx = max.col(-t(alpha.inf)) # get the indce of the min of each row
				alphaMin = alpha.inf[minIdx + (0:(nHset-1) * lVar)]
				alpha[,1:nHset] = matrix(alphaMin, lVar, nHset, byrow=TRUE);
				D[,Hset] = D[,Hset, drop=FALSE] - alpha[,1:nHset, drop=FALSE] * (D[,Hset, drop=FALSE]-K[,Hset, drop=FALSE]);			
				idx2zero = (Hset - 1) * lVar + minIdx; # convert array indices to index relative to the matrix D
				D[idx2zero] = 0;
				Pset[idx2zero] = FALSE;
				K[, Hset] = .cssls(CtC, CtA[,Hset, drop=FALSE], Pset[,Hset, drop=FALSE]);
				# which column of K have at least one negative entry?
				Hset = which( colSums(K < 0) > 0 );
				nHset = length(Hset);
			}
		}
		#V-#
		
		#Vc# Make sure the solution has converged
		#if iter == maxiter, error('Maximum number iterations exceeded'), end
		# Check solutions for optimality
		W[,Fset] = CtA[,Fset, drop=FALSE] - CtC %*% K[,Fset, drop=FALSE];
		# which columns have all entries non-positive
		Jset = which( colSums( (ifelse(!(Pset[,Fset, drop=FALSE]),1,0) * W[,Fset, drop=FALSE]) > 0 ) == 0 );
		Fset = setdiff(Fset, Fset[Jset]);
		
		if ( length(Fset) > 0 ){				
			#Vc# For non-optimal solutions, add the appropriate variable to Pset						
			# get indice of the maximum in each column
			mxidx = max.col( t(ifelse(!Pset[,Fset, drop=FALSE],1,0) * W[,Fset, drop=FALSE]) )
			Pset[ (Fset - 1) * lVar + mxidx ] = TRUE;
			D[,Fset] = K[,Fset, drop=FALSE];
		}		
	}
	#V-#
	
	# return K and Pset
	list(K=K, Pset=Pset)
}
# ****************************** Subroutine****************************
#library(corpcor)
.cssls <- function(CtC, CtA, Pset=NULL){
	# Solve the set of equations CtA = CtC*K for the variables in set Pset
	# using the fast combinatorial approach
	K = matrix(0, nrow(CtA), ncol(CtA));	
	if ( is.null(Pset) || length(Pset)==0 || all(Pset) ){		
		K = solve(CtC) %*% CtA;
		# K = pseudoinverse(CtC) %*% CtA;
		#K=pinv(CtC)*CtA;
	}else{
		lVar = nrow(Pset); pRHS = ncol(Pset);
		codedPset = as.numeric(2.^(seq(lVar-1,0,-1)) %*% Pset);
		sortedPset = sort(codedPset)
		sortedEset = order(codedPset)
		breaks = diff(sortedPset);
		breakIdx = c(0, which(breaks > 0 ), pRHS);
		for( k in seq(1,length(breakIdx)-1) ){
			cols2solve = sortedEset[ seq(breakIdx[k]+1, breakIdx[k+1])];
			vars = Pset[,sortedEset[breakIdx[k]+1]];			
			K[vars,cols2solve] = solve(CtC[vars,vars, drop=FALSE]) %*% CtA[vars,cols2solve, drop=FALSE];			
			#K[vars,cols2solve] = pseudoinverse(CtC[vars,vars]) %*% CtA[vars,cols2solve];
			#TODO: check if this is the right way or needs to be reversed
			#K(vars,cols2solve) = pinv(CtC(vars,vars))*CtA(vars,cols2solve);
		}
	}
	
	# return K
	K
}

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
.nmfsh_comb <- function(A, k, nmf.fit, eta=-1, beta=0.01, bi_conv=c(0, 10), eps_conv=1e-4, version=c('R', 'L'), verbose=FALSE){
		
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
	if( missing(nmf.fit) ){
		message('Init W internally (random)')
		W=matrix(runif(m*k), m,k);	
		nmf.fit <- NULL
	} else {
		
		# seed the method (depends on the version to run)
		start <- if( version == 'R' ) basis(nmf.fit) else t(coef(nmf.fit))
		# check compatibility of the starting point with the target matrix
		if( any(dim(start) != c(m,k)) ) 
			stop("Invalid initialization: incompatible dimensions [expected: ", paste(c(m,k), collapse=' x '),", got: ", paste(dim(start), collapse=' x '), " ]")	
		# use the supplied starting point
		W <- start
		
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
		if( !is.null(nmf.fit) ) 
			nmf.fit <- trackError(nmf.fit, .snmf.objective(A, W, H, eta, beta), i)
		
		# test convergence every 5 iterations
		if ( (i %% 5==0)  || (i==1) ){
			# indice of maximum for each row of W
			idxW = max.col(W)
			# indice of maximum for each column of H
			idxH = max.col(t(H))
			changedW=sum(idxW != idxWold); changedH=sum(idxH != idxHold);
			if ( (changedW<=wminchange) && (changedH==0) ) inc=inc+1
			else inc=0

			resmat=pmin(H, crossprod(W) %*% H - t(W) %*% A + matrix(beta, k , k) %*% H); resvec=as.numeric(resmat);
			resmat=pmin(W, W %*% tcrossprod(H) - A %*% t(H) + eta2 * W); resvec=c(resvec, as.numeric(resmat));
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
	if( !is.null(nmf.fit) ) 
		nmf.fit <- trackError(nmf.fit, .snmf.objective(A, W, H, eta, beta), i, force=TRUE)	

	# transpose and reswap the roles
	if( !is.null(nmf.fit) ){ 
		if( version == 'L' ){
			basis(nmf.fit) <- t(H)
			coef(nmf.fit) <- t(W)
		}
		else{ basis(nmf.fit) <- W; coef(nmf.fit) <- H}
		# set number of iterations performed
		nmf.fit$iteration <- i
		
		return(nmf.fit)	
	}else{
		res <- list(W=W, H=H)
		if( version == 'L' ){
			res$W <- t(H)
			res$H <- t(W)
		}		
		return(invisible(res))
	}
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
	sol <- .nmfsh_comb(target, nbasis(seed), nmf.fit=seed, version=version, verbose=verbose(seed), ...)
		
	# return solution
	return(sol)
}

# internal function to register the methods
.nmf.plugin.snmf <- function(){
	list(newNMFStrategy(.snmf, 'snmf/r', objective='euclidean')
		, newNMFStrategy(.snmf, 'snmf/l', objective='euclidean')
	)
}
