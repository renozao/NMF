#' @include registry-algorithms.R
NULL

#' Fast Combinatorial Nonnegative Least-Square
#' 
#' This function solves the following nonnegative least square linear problem
#' using normal equations and the fast combinatorial strategy from \cite{VanBenthem2004}:
#' 
#' \deqn{
#'  \begin{array}{l}
#'  \min \|Y - X K\|_F\\
#'  \mbox{s.t. } K>=0
#'  \end{array}
#' }{min ||Y - X K||_F, s.t. K>=0}
#' 
#' where \eqn{Y} and \eqn{X} are two real matrices of dimension \eqn{n \times p}{n x p} 
#' and \eqn{n \times r}{n x r} respectively, 
#' and \eqn{\|.\|_F}{|.|_F} is the Frobenius norm.
#'
#' The algorithm is very fast compared to other approaches, as it is optimised 
#' for handling multiple right-hand sides.
#' 
#' @details
#' Within the \code{NMF} package, this algorithm is used internally by the
#' SNMF/R(L) algorithm from \cite{KimH2007} to solve general Nonnegative
#' Matrix Factorization (NMF) problems, using alternating nonnegative
#' constrained least-squares. 
#' That is by iteratively and alternatively estimate each matrix factor.
#'
#' The algorithm is an active/passive set method, which rearrange the 
#' right-hand side to reduce the number of pseudo-inverse calculations.  
#' It uses the unconstrained solution \eqn{K_u} obtained from the 
#' unconstrained least squares problem,
#' i.e. \eqn{\min \|Y - X K\|_F^2}{min ||Y - X K||_F^2} , so as to determine
#' the initial passive sets.
#'  
#' The function \code{fcnnls} is provided separately so that it can be 
#' used to solve other types of nonnegative least squares problem. 
#' For faster computation, when multiple nonnegative least square fits 
#' are needed, it is recommended to directly use the function \code{\link{.fcnnls}}.
#' 
#' The code of this function is a port from the original MATLAB code 
#' provided by \cite{KimH2007}.
#' 
#' @inheritParams .fcnnls
#' @param ...  extra arguments passed to the internal function \code{.fcnnls}.
#' Currently not used.
#' @return A list containing the following components:
#' 
#' \item{x}{ the estimated optimal matrix \eqn{K}.} \item{fitted}{ the fitted
#' matrix \eqn{X K}.} \item{residuals}{ the residual matrix \eqn{Y - X K}.}
#' \item{deviance}{ the residual sum of squares between the fitted matrix
#' \eqn{X K} and the target matrix \eqn{Y}. That is the sum of the square
#' residuals.} \item{passive}{ a \eqn{r x p} logical matrix containing the
#' passive set, that is the set of entries in \eqn{K} that are not null (i.e.
#' strictly positive).} \item{pseudo}{ a logical that is \code{TRUE} if the
#' computation was performed using the pseudoinverse. See argument
#' \code{pseudo}.}
#' 
#' @seealso \code{\link{nmf}}
#' @references 
#' 
#' Original MATLAB code from Van Benthem and Keenan, slightly modified by H.
#' Kim:\cr \url{http://www.cc.gatech.edu/~hpark/software/fcnnls.m}
#' 
#' @author 
#' Original MATLAB code : Van Benthem and Keenan
#' 
#' Adaption of MATLAB code for SNMF/R(L): H. Kim
#' 
#' Adaptation to the NMF package framework: Renaud Gaujoux
#' 
#' @keywords optimize multivariate regression
#' @export 
#' @inline
#' @examples
#' 
#' ## Define a random nonnegative matrix matrix
#' n <- 200; p <- 20; r <- 3
#' V <- rmatrix(n, p)
#' 
#' ## Compute the optimal matrix K for a given X matrix
#' X <- rmatrix(n, r)
#' res <- fcnnls(X, V)
#' 
#' ## Compute the same thing using the Moore-Penrose generalized pseudoinverse
#' res <- fcnnls(X, V, pseudo=TRUE)
#' 
#' ## It also works in the case of single vectors
#' y <- runif(n)
#' res <- fcnnls(X, y)
#' # or
#' res <- fcnnls(X[,1], y)
#' 
#' 
setGeneric('fcnnls', function(x, y, ...) standardGeneric('fcnnls') )
#' This method wraps a call to the internal function \code{.fcnnls}, and
#' formats the results in a similar way as other lest-squares methods such 
#' as \code{\link{lm}}. 
#' 
#' @param verbose toggle verbosity (default is \code{FALSE}).
#' 
setMethod('fcnnls', signature(x='matrix', y='matrix'), 
	function(x, y, verbose=FALSE, pseudo=TRUE, ...){
		# load corpcor if necessary
		if( isTRUE(pseudo) ){ 
			library(corpcor)
		}
		
		# call the internal function
		res <- .fcnnls(x, y, verbose=verbose, pseudo=pseudo, ...)
		
		# process the result
		f <- x %*% res$coef
		resid <- y - f
        # set dimnames
        if( is.null(rownames(res$coef)) ) rownames(res$coef) <- colnames(x)

		# wrap up the result
		out <- list(x=res$coef, fitted=f, residuals=resid, deviance=norm(resid, 'F')^2, passive=res$Pset, pseudo=pseudo)
		class(out) <- 'fcnnls'
		out
	}
)
#' Shortcut for \code{fcnnls(as.matrix(x), y, ...)}.
setMethod('fcnnls', signature(x='numeric', y='matrix'), 
	function(x, y, ...){
		fcnnls(as.matrix(x), y, ...)
	}
)
#' Shortcut for \code{fcnnls(x, as.matrix(y), ...)}.
setMethod('fcnnls', signature(x='ANY', y='numeric'), 
	function(x, y, ...){
		fcnnls(x, as.matrix(y), ...)
	}
)

#' @S3method print fcnnls
print.fcnnls <- function(x, ...){
	cat("<object of class 'fcnnls': Fast Combinatorial Nonnegative Least Squares>\n")
	cat("Dimensions:", nrow(x$x)," x ", ncol(x$x), "\n")
	cat("Residual sum of squares:", x$deviance,"\n")
	cat("Active constraints:", length(x$passive)-sum(x$passive),"/", length(x$passive), "\n")
	cat("Inverse method:", 
			if( isTRUE(x$pseudo) ) 'pseudoinverse (corpcor)'
			else if( is.function(x$pseudo) ) str_fun(x$pseudo) 
			else 'QR (solve)', "\n")
	invisible(x)
}


###% M. H. Van Benthem and M. R. Keenan, J. Chemometrics 2004; 18: 441-450
###%
###% Given A and C this algorithm solves for the optimal 
###% K in a least squares sense, using that
###%      A = C*K 
###% in the problem
###%      min ||A-C*K||, s.t. K>=0, for given A and C.
###%
###%
###% @param C the matrix of coefficients
###% @param A the target matrix of observations
###%
###% @return [K, Pset]
###%
#' Internal Routine for Fast Combinatorial Nonnegative Least-Squares
#'
#' @description 
#' This is the workhorse function for the higher-level function 
#' \code{\link{fcnnls}}, which implements the fast nonnegative least-square 
#' algorithm for multiple right-hand-sides from \cite{VanBenthem2004} to solve 
#' the following problem:
#' 
#' \deqn{
#'  \begin{array}{l}
#'  \min \|Y - X K\|_F\\
#'  \mbox{s.t. } K>=0
#'  \end{array}
#' }{min ||Y - X K||_F, s.t. K>=0}
#' 
#' where \eqn{Y} and \eqn{X} are two real matrices of dimension \eqn{n \times p}{n x p} 
#' and \eqn{n \times r}{n x r} respectively, 
#' and \eqn{\|.\|_F}{|.|_F} is the Frobenius norm.
#'
#' The algorithm is very fast compared to other approaches, as it is optimised 
#' for handling multiple right-hand sides.
#' 
#' @param x the coefficient matrix
#' @param y the target matrix to be approximated by \eqn{X K}.
#' @param verbose logical that indicates if log messages should be shown.
#' @param pseudo By default (\code{pseudo=FALSE}) the algorithm uses Gaussian
#' elimination to solve the successive internal linear problems, using the
#' \code{\link{solve}} function.  If \code{pseudo=TRUE} the algorithm uses
#' Moore-Penrose generalized \code{\link[corpcor]{pseudoinverse}} from the
#' \code{corpcor} package instead of \link{solve}.
#' @param eps threshold for considering entries as nonnegative.
#' This is an experimental parameter, and it is recommended to 
#' leave it at 0.
#' 
#' @return A list with the following elements: 
#' 
#' \item{coef}{the fitted coefficient matrix.} 
#' \item{Pset}{the set of passive constraints, as a logical matrix of 
#' the same size as \code{K} that indicates which element is positive.}
#'  
#' @export
.fcnnls <- function(x, y, verbose=FALSE, pseudo=FALSE, eps=0){
	
	# check arguments
	if( any(dim(y) == 0L) ){
		stop("Empty target matrix 'y' [", paste(dim(y), collapse=' x '), "]")
	}
	if( any(dim(x) == 0L) ){
		stop("Empty regression variable matrix 'x' [", paste(dim(x), collapse=' x '), "]")
	}
	
	# map arguments
	C <- x
	A <- y
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
	K = .cssls(CtC, CtA, pseudo=pseudo);
	Pset = K > 0;
	K[!Pset] = 0;
	D = K;
	# which columns of Pset do not have all entries TRUE?
	Fset = which( colSums(Pset) != lVar );
	#V+# Active set algorithm for NNLS main loop
	oitr=0; # HKim
	while ( length(Fset)>0 ) {
		
		oitr=oitr+1; if ( verbose && oitr > 5 ) cat(sprintf("%d ",oitr));# HKim
		
		#Vc# Solve for the passive variables (uses subroutine below)				
		K[,Fset] = .cssls(CtC, CtA[,Fset, drop=FALSE], Pset[,Fset, drop=FALSE], pseudo=pseudo);
		
		# Find any infeasible solutions
		# subset Fset on the columns that have at least one negative entry
		Hset = Fset[ colSums(K[,Fset, drop=FALSE] < eps) > 0 ];
		#V+# Make infeasible solutions feasible (standard NNLS inner loop)
		if ( length(Hset)>0 ){
			nHset = length(Hset);
			alpha = matrix(0, lVar, nHset);
			while ( nHset>0  && (iter < maxiter) ){
				iter = iter + 1; 
				alpha[,1:nHset] = Inf;
				#Vc# Find indices of negative variables in passive set
				ij = which( Pset[,Hset, drop=FALSE] & (K[,Hset, drop=FALSE] < eps) , arr.ind=TRUE);			
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
				K[, Hset] = .cssls(CtC, CtA[,Hset, drop=FALSE], Pset[,Hset, drop=FALSE], pseudo=pseudo);
				# which column of K have at least one negative entry?
				Hset = which( colSums(K < eps) > 0 );
				nHset = length(Hset);
			}
		}
		#V-#
		
		#Vc# Make sure the solution has converged
		#if iter == maxiter, error('Maximum number iterations exceeded'), end
		# Check solutions for optimality
		W[,Fset] = CtA[,Fset, drop=FALSE] - CtC %*% K[,Fset, drop=FALSE];
		# which columns have all entries non-positive
		Jset = which( colSums( (ifelse(!(Pset[,Fset, drop=FALSE]),1,0) * W[,Fset, drop=FALSE]) > eps ) == 0 );
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
	list(coef=K, Pset=Pset)
}
# ****************************** Subroutine****************************
#library(corpcor)
.cssls <- function(CtC, CtA, Pset=NULL, pseudo=FALSE){
	
	# use provided function
	if( is.function(pseudo) ){
		pseudoinverse <- pseudo
		pseudo <- TRUE
	}
	
	# Solve the set of equations CtA = CtC*K for the variables in set Pset
	# using the fast combinatorial approach
	K = matrix(0, nrow(CtA), ncol(CtA));	
	if ( is.null(Pset) || length(Pset)==0 || all(Pset) ){		
		K <- (if( !pseudo ) solve(CtC) else pseudoinverse(CtC)) %*% CtA;
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
			K[vars,cols2solve] <- (if( !pseudo ) solve(CtC[vars,vars, drop=FALSE]) else pseudoinverse(CtC[vars,vars, drop=FALSE])) %*% CtA[vars,cols2solve, drop=FALSE];
			#K[vars,cols2solve] <-  pseudoinverse(CtC[vars,vars, drop=FALSE])) %*% CtA[vars,cols2solve, drop=FALSE];
			#TODO: check if this is the right way or needs to be reversed
			#K(vars,cols2solve) = pinv(CtC(vars,vars))*CtA(vars,cols2solve);
		}
	}
	
	# return K
	K
}

###%
###% SNMF/R  
###%
###% Author: Hyunsoo Kim and Haesun Park, Georgia Insitute of Technology
###%
###% Reference: 
###%
###%   Sparse Non-negative Matrix Factorizations via Alternating 
###%   Non-negativity-constrained Least Squares for Microarray Data Analysis
###%   Hyunsoo Kim and Haesun Park, Bioinformatics, 2007, to appear.
###%
###% This software requires fcnnls.m, which can be obtained from 
###% M. H. Van Benthem and M. R. Keenan, J. Chemometrics 2004; 18: 441-450
###%
###% NMF: min_{W,H} (1/2) || A - WH ||_F^2 s.t. W>=0, H>=0 
###% SNMF/R: NMF with additional sparsity constraints on H
###%
###%   min_{W,H} (1/2) (|| A - WH ||_F^2 + eta ||W||_F^2 
###%                + beta (sum_(j=1)^n ||H(:,j)||_1^2))
###%                s.t. W>=0, H>=0 
###%
###% A: m x n data matrix (m: features, n: data points)
###% W: m x k basis matrix
###% H: k x n coefficient matrix
###%
###% function [W,H,i]=nmfsh_comb(A,k,param,verbose,bi_conv,eps_conv)
###%
###% input parameters:
###%   A: m x n data matrix (m: features, n: data points)
###%   k: desired positive integer k
###%   param=[eta beta]:  
###%      eta (for supressing ||W||_F)
###%         if eta < 0, software uses maxmum value in A as eta. 
###%      beta (for sparsity control)
###%         Larger beta generates higher sparseness on H.
###%         Too large beta is not recommended. 
###%   verbos: verbose = 0 for silence mode, otherwise print output
###%   eps_conv: KKT convergence test (default eps_conv = 1e-4)
###%   bi_conv=[wminchange iconv] biclustering convergence test 
###%        wminchange: the minimal allowance of the change of 
###%        row-clusters  (default wminchange=0)
###%        iconv: decide convergence if row-clusters (within wminchange)
###%        and column-clusters have not changed for iconv convergence 
###%        checks. (default iconv=10)
###%
###% output:
###%   W: m x k basis matrix
###%   H: k x n coefficient matrix
###%   i: the number of iterations
###%
###% sample usage:
###%  [W,H]=nmfsh_comb(amlall,3,[-1 0.01],1);
###%  [W,H]=nmfsh_comb(amlall,3,[-1 0.01],1,[3 10]); 
###%     -- in the convergence check, the change of row-clusters to
###%        at most three rows is allowed.
###%
###% 
#function [W,H,i] 
nmf_snmf <- function(A, x, maxIter= nmf.getOption('maxIter') %||% 20000L, eta=-1, beta=0.01, bi_conv=c(0, 10), eps_conv=1e-4, version=c('R', 'L'), verbose=FALSE){
#nmfsh_comb <- function(A, k, param, verbose=FALSE, bi_conv=c(0, 10), eps_conv=1e-4, version=c('R', 'L')){
	
	# depending on the version: 
	# in version L: A is transposed while W and H are swapped and transposed
	version <- match.arg(version)
	if( version == 'L' ) A <- t(A) 
	#if( missing(param) ) param <- c(-1, 0.01)
	
	m = nrow(A); n = ncol(A); erravg1 = numeric();
	
	#eta=param[1]; beta=param[2]; 
	maxA=max(A); if ( eta<0 ) eta=maxA;
	eta2=eta^2;
	
	# bi_conv
	if( length(bi_conv) != 2 )
		stop("SNMF/", version, "::Invalid argument 'bi_conv' - value should be a 2-length numeric vector")
	wminchange=bi_conv[1]; iconv=bi_conv[2];
	
	## VALIDITY of parameters
	# eps_conv
	if( eps_conv <= 0 )
		stop("SNMF/", version, "::Invalid argument 'eps_conv' - value should be positive")
	# wminchange
	if( wminchange < 0 )
		stop("SNMF/", version, "::Invalid argument 'bi_conv' - bi_conv[1] (i.e 'wminchange') should be non-negative")
	# iconv
	if( iconv < 0 )
		stop("SNMF/", version, "::Invalid argument 'bi_conv' - bi_conv[2] (i.e 'iconv') should be non-negative")
	# beta
	if( beta <=0 )
		stop("SNMF/", version, "::Invalid argument 'beta' - value should be positive")
	##
	
	# initialize random W if no starting point is given
	if( isNumber(x) ){
		# rank is given by x
		k <- x
		message('# NOTE: Initialise W internally (runif)')
		W <- matrix(runif(m*k), m,k);	
		x <- NULL
	} else if( is.nmf(x) ){
		# rank is the number of basis components in x
		k <- nbasis(x)
		# seed the method (depends on the version to run)
		start <- if( version == 'R' ) basis(x) else t(coef(x))
		# check compatibility of the starting point with the target matrix
		if( any(dim(start) != c(m,k)) )
			stop("SNMF/", version, " - Invalid initialization - incompatible dimensions [expected: ", paste(c(m,k), collapse=' x '),", got: ", paste(dim(start), collapse=' x '), " ]")	
		# use the supplied starting point
		W <- start
	}else{
		stop("SNMF/", version, ' - Invalid argument `x`: must be a single numeric or an NMF model [', class(x), ']')
	}
	
	if ( verbose )
		cat(sprintf("--\nAlgorithm: SNMF/%s\nParameters: k=%d eta=%.4e beta (for sparse H)=%.4e wminchange=%d iconv=%d\n",
				version, k,eta,beta,wminchange,iconv));

	idxWold=rep(0, m); idxHold=rep(0, n); inc=0;
		
	# check validity of seed
	if( any(NAs <- is.na(W)) )
		stop("SNMF/", version, "::Invalid initialization - NAs found in the ", if(version=='R') 'basis (W)' else 'coefficient (H)' , " matrix [", sum(NAs), " NAs / ", length(NAs), " entries]")
	
	# normalize columns of W
	W= apply(W, 2, function(x) x / sqrt(sum(x^2)) );	

	I_k=diag(eta, k); betavec=rep(sqrt(beta), k); nrestart=0;
	i <- 0L
	while( i < maxIter){
		i <- i + 1L
		
		# min_h ||[[W; 1 ... 1]*H  - [A; 0 ... 0]||, s.t. H>=0, for given A and W.
	  	res = .fcnnls(rbind(W, betavec), rbind(A, rep(0, n)));	  	  	
		H = res[[1]]

		if ( any(rowSums(H)==0) ){
			if( verbose ) cat(sprintf("iter%d: 0 row in H eta=%.4e restart!\n",i,eta));
			nrestart=nrestart+1;
			if ( nrestart >= 10 ){
				warning("NMF::snmf - Too many restarts due to too big 'beta' value [Computation stopped after the 9th restart]");
				break;
			}
						
			# re-initialize random W
			idxWold=rep(0, m); idxHold=rep(0, n); inc=0; 
			erravg1 <- numeric();# re-initialize base average error
			W=matrix(runif(m*k), m,k);
			W= apply(W, 2, function(x) x / sqrt(sum(x^2)) );  # normalize columns of W	
			next;
		}
		
		# min_w ||[H'; I_k]*W' - [A'; 0]||, s.t. W>=0, for given A and H. 
		res = .fcnnls(rbind(t(H), I_k), rbind(t(A), matrix(0, k,m))); 
		Wt = res[[1]]
		W= t(Wt);		

		# track the error (not computed unless tracking option is enabled in x)
		if( !is.null(x) ) 
			x <- trackError(x, .snmf.objective(A, W, H, eta, beta), niter=i)
		
		# test convergence every 5 iterations OR if the base average error has not been computed yet
		if ( (i %% 5==0)  || (length(erravg1)==0) ){
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
			# compute base average error if necessary
			if ( length(erravg1)==0 )
				erravg1=erravg;
			
			if ( verbose && (i %% 1000==0) ){ # prints number of changing elements
				if( i==1000 ) cat("Track:\tIter\tInc\tchW\tchH\t---\terravg1\terravg\terravg/erravg1\n")
				cat(sprintf("\t%d\t%d\t%d\t%d\t---\terravg1: %.4e\terravg: %.4e\terravg/erravg1: %.4e\n",
					 i,inc,changedW,changedH,erravg1,erravg,erravg/erravg1));
			}
			
			#print(list(inc=inc, iconv=iconv, erravg=erravg, eps_conv=eps_conv, erravg1=erravg1))
			if ( (inc>=iconv) && (erravg<=eps_conv*erravg1) ) break;
			idxWold=idxW; idxHold=idxH; 
		}

	}
	
	if( verbose ) cat("--\n")
	
	# force to compute last error if not already done
	if( !is.null(x) ) 
		x <- trackError(x, .snmf.objective(A, W, H, eta, beta), niter=i, force=TRUE)	

	# transpose and reswap the roles
	if( !is.null(x) ){ 
		if( version == 'L' ){
			.basis(x) <- t(H)
			.coef(x) <- t(W)
		}
		else{ 
			.basis(x) <- W 
			.coef(x) <- H
		}
		# set number of iterations performed
		niter(x) <- i
		
		return(x)	
	}else{
		res <- list(W=W, H=H)
		if( version == 'L' ){
			res$W <- t(H)
			res$H <- t(W)
		}		
		return(invisible(res))
	}
}


###% Computes the objective value for the SNMF algorithm
.snmf.objective <- function(target, w, h, eta, beta){
	
	1/2 * ( sum( (target - (w %*% h))^2 ) 
				+ eta * sum(w^2) 
				+ beta * sum( colSums( h )^2 )
				)
}

snmf.objective <- function(x, y, eta=-1, beta=0.01){
	.snmf.objective(y, .basis(x), .coef(x), eta, beta)
}

###% Wrapper function to use the SNMF/R algorithm with the NMF package.
###%
.snmf <- function(target, seed, maxIter=20000L, eta=-1, beta=0.01, bi_conv=c(0, 10), eps_conv=1e-4, ...){	
	
	# retrieve the version of SNMF algorithm from its name: 
	# it is defined by the last letter in the method's name (in upper case)
	name <- algorithm(seed)
	version <- toupper(substr(name, nchar(name), nchar(name)))
	
	# perform factorization using Kim and Park's algorithm
	ca <- match.call()
	ca[[1L]] <- as.name('nmf_snmf')
	# target
	ca[['A']] <- ca[['target']]
	ca[['target']] <- NULL
	# seed
	ca[['x']] <- ca[['seed']]
	ca[['seed']] <- NULL
	# version 
	ca[['version']] <- version
	# verbose
	ca[['verbose']] <- verbose(seed)
	e <- parent.frame()
	sol <- eval(ca, envir=e)
#	nmf_snmf(target, seed, ..., version = version, verbose = verbose(seed))
		
	# return solution
	return(sol)
}

#' NMF Algorithm - Sparse NMF via Alternating NNLS
#' 
#' NMF algorithms proposed by \cite{KimH2007} that enforces sparsity 
#' constraint on the basis matrix (algorithm \sQuote{SNMF/L}) or the 
#' mixture coefficient matrix (algorithm \sQuote{SNMF/R}).
#' 
#' The algorithm \sQuote{SNMF/R} solves the following NMF optimization problem on 
#' a given target matrix \eqn{A} of dimension \eqn{n \times p}{n x p}:
#' \deqn{
#' \begin{array}{ll}
#' & \min_{W,H} \frac{1}{2} \left(|| A - WH ||_F^2 + \eta ||W||_F^2 
#' 	               + \beta (\sum_{j=1}^p ||H_{.j}||_1^2)\right)\\
#'                 s.t. & W\geq 0, H\geq 0
#' \end{array}
#' }{
#' min_{W,H} 1/2 (|| A - WH ||_F^2 + eta ||W||_F^2 
#' 	               + beta (sum_j ||H[,j]||_1^2))
#' 
#'                 s.t. W>=0, H>=0 
#' }
#' 
#' The algorithm \sQuote{SNMF/L} solves a similar problem on the transposed target matrix \eqn{A},
#' where \eqn{H} and \eqn{W} swap roles, i.e. with sparsity constraints applied to \code{W}. 
#' 
#' @param maxIter maximum number of iterations. 
#' @param eta parameter to suppress/bound the L2-norm of \code{W} and in 
#' \code{H} in \sQuote{SNMF/R} and \sQuote{SNMF/L} respectively.
#' 
#' If \code{eta < 0}, then it is set to the maximum value in the target matrix is used.   
#' @param beta regularisation parameter for sparsity control, which 
#' balances the trade-off between the accuracy of the approximation and the 
#' sparseness of \code{H} and \code{W} in \sQuote{SNMF/R} and \sQuote{SNMF/L} respectively.
#' 
#' Larger beta generates higher sparseness on \code{H} (resp. \code{W}).
#' Too large beta is not recommended.
#' @param bi_conv parameter of the biclustering convergence test.
#' It must be a size 2 numeric vector \code{bi_conv=c(wminchange, iconv)}, 
#' with:
#' \describe{
#' \item{\code{wminchange}:}{the minimal allowance of change in row-clusters.}
#' \item{\code{iconv}:}{ decide convergence if row-clusters 
#' (within the allowance of \code{wminchange})
#' and column-clusters have not changed for \code{iconv} convergence checks.}
#' } 
#' 
#' Convergence checks are performed every 5 iterations.
#' @param eps_conv threshold for the KKT convergence test. 
#' @param ... extra argument not used.
#' 
#' @rdname SNMF-nmf
#' @aliases SNMF/R-nmf
nmfAlgorithm.SNMF_R <- setNMFMethod('snmf/r', .snmf, objective=snmf.objective)
#' @aliases SNMF/L-nmf
#' @rdname SNMF-nmf
nmfAlgorithm.SNMF_L <- setNMFMethod('snmf/l', .snmf, objective=snmf.objective)
