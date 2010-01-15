#V1# threshold=-1

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
	CtC = t(C)%*%C; CtA = t(C)%*%A;
	# Obtain the initial feasible solution and corresponding passive set
	K = .cssls(CtC, CtA);
	Pset = K > 0;
	K[!Pset] = 0;
	D = K;
	Fset = which( apply(Pset, 2, function(x) !all(x)) );
	#V+# Active set algorithm for NNLS main loop
	oitr=0; # HKim
	while ( length(Fset)>0 ) {
		
		oitr=oitr+1; if ( oitr > 5 ) cat(sprintf("%d ",oitr));# HKim

		#Vc# Solve for the passive variables (uses subroutine below)				
		K[,Fset] = .cssls(CtC, CtA[,Fset, drop=FALSE], Pset[,Fset, drop=FALSE]);

		# Find any infeasible solutions
		Hset = Fset[ apply(K[,Fset, drop=FALSE], 2, function(x) any(x < 0)) ];
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
				nK = nrow(K);
				negIdx = (Hset[j] - 1) * nK + i; # convert array indices to index relative to the matrix K (i.e. same row index but col index is stored in Hset)
				
				alpha[hIdx] = D[negIdx] / (D[negIdx] - K[negIdx]);				
				alphaMin = t(apply(alpha[,1:nHset, drop=FALSE], 2, function(x){ idx <- which.min(x); c(x[idx], idx)}))
				minIdx = alphaMin[,2]; alphaMin = alphaMin[,1]; # alphaMin is a vector of length nHset            
				alpha[,1:nHset] = matrix(alphaMin, lVar, nHset, byrow=TRUE);
				D[,Hset] = D[,Hset, drop=FALSE] - alpha[,1:nHset, drop=FALSE] * (D[,Hset, drop=FALSE]-K[,Hset, drop=FALSE]);			
				nD = nrow(D);
				idx2zero = (Hset - 1) * nD + minIdx; # convert array indices to index relative to the matrix D
				D[idx2zero] = 0;
				Pset[idx2zero] = FALSE;
				K[, Hset] = .cssls(CtC, CtA[,Hset, drop=FALSE], Pset[,Hset, drop=FALSE]);
				Hset = which( apply(K, 2, function(x) any(x < 0)) ); nHset = length(Hset);
			}
		}
		#V-#
		   
		#Vc# Make sure the solution has converged
		#if iter == maxiter, error('Maximum number iterations exceeded'), end
		# Check solutions for optimality
		W[,Fset] = CtA[,Fset, drop=FALSE] - CtC %*% K[,Fset, drop=FALSE];
		Jset = which( apply(ifelse(!(Pset[,Fset, drop=FALSE]),1,0) * W[,Fset, drop=FALSE], 2, function(x) all(x <= 0)) );
		Fset = setdiff(Fset, Fset[Jset]);
		
		if ( length(Fset) > 0 ){				
			#Vc# For non-optimal solutions, add the appropriate variable to Pset						
			mxidx = apply(ifelse(!Pset[,Fset, drop=FALSE],1,0) * W[,Fset, drop=FALSE], 2, function(x) which.max(x) )			
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
