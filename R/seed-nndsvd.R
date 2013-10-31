#' @include registry-seed.R
NULL

###% Seeding method: Nonnegative Double Singular Value Decomposition
###%
###% @author Renaud Gaujoux
###% @creation 17 Jul 2009


###% Auxliary functions
.pos <- function(x){ as.numeric(x>=0) * x }
.neg <- function(x){ - as.numeric(x<0) * x }
.norm <- function(x){ sqrt(drop(crossprod(x))) }


###% This function implements the NNDSVD algorithm described in Boutsidis (2008) for
###% initializattion of Nonnegative Matrix Factorization Algorithms.
###% 
###% @param A the input nonnegative m x n matrix A  
###% @param k the rank of the computed factors W,H
###% @param flag indicates the variant of the NNDSVD Algorithm:
###%        - flag=0 --> NNDSVD
###%        - flag=1 --> NNDSVDa
###%        - flag=2 --> NNDSVDar
###%
###% @note This code is a port from the MATLAB code from C. Boutsidis and E. Gallopoulos kindly provided by the authors for research purposes.
###% Original MATLAB code: http://www.cs.rpi.edu/~boutsc/papers/paper1/nndsvd.m
###% 
###% @references 	C. Boutsidis and E. Gallopoulos, 
###% SVD-based initialization: A head start for nonnegative matrix factorization, 
###% Pattern Recognition, 2007
###% doi:10.1016/j.patcog.2007.09.010
###%
.nndsvd.wrapper <- function(object, x, densify=c('none', 'average', 'random')){
	
	# match parameter 'densify'
	densify <- match.arg(densify)
	flag <- which(densify == c('none', 'average', 'random')) - 1
	res <- .nndsvd.internal(x, nbasis(object), flag)
	
	# update 'NMF' object
	.basis(object) <- res$W; .coef(object) <- res$H	
	
	# return updated object
	object
}

###% Port to R of the MATLAB code from Boutsidis
.nndsvd.internal <- function(A, k, flag=0, LINPACK=FALSE){
	
	#check the input matrix
	if( any(A<0) ) stop('The input matrix contains negative elements !')
	
	#size of input matrix
	size = dim(A);
	m <- size[1]; n<- size[2]
	
	#the matrices of the factorization
	W = matrix(0, m, k);
	H = matrix(0, k, n);
	
	#1st SVD --> partial SVD rank-k to the input matrix A.	
	s = svd(A, k, k, LINPACK=LINPACK);
	U <- s$u; S <- s$d; V <- s$v
	
	#-------------------------------------------------------
	# We also recommend the use of propack for the SVD
	# 1st SVD --> partial SVD rank-k ( propack )
	# OPTIONS.tol  = 0.00001;               % remove comment to this line
	# [U,S,X] = LANSVD(A,k,'L',OPTIONS);    % remove comment to this line 
	#-------------------------------------------------------
			
	#choose the first singular triplet to be nonnegative
	W[,1] = sqrt(S[1]) * abs(U[,1]);         
	H[1,] = sqrt(S[1]) * abs(t(V[,1])); 
					
	# second SVD for the other factors (see table 1 in Boutsidis' paper)
	for( i in seq(2,k) ){
		uu = U[,i]; vv = V[,i];
		uup = .pos(uu); uun = .neg(uu) ;
		vvp = .pos(vv); vvn = .neg(vv);
		n_uup = .norm(uup);
		n_vvp = .norm(vvp) ;
		n_uun = .norm(uun) ;
		n_vvn = .norm(vvn) ;
		termp = n_uup %*% n_vvp; termn = n_uun %*% n_vvn;
		if (termp >= termn){
			W[,i] = sqrt(S[i] * termp) * uup / n_uup; 
			H[i,] = sqrt(S[i] * termp) * vvp / n_vvp;
		}else{		
			W[,i] = sqrt(S[i] * termn) * uun / n_uun; 
			H[i,] = sqrt(S[i] * termn) * vvn / n_vvn;
		}
	}
	
	#------------------------------------------------------------
			
	#actually these numbers are zeros
	W[W<0.0000000001] <- 0;
	H[H<0.0000000001] <- 0;
		
	if( flag==1 ){ #NNDSVDa: fill in the zero elements with the average		
		
		ind1 <- W==0 ;
		ind2 <- H==0 ;
		average <- mean(A); 
		W[ind1] <- average; 
		H[ind2] <- average;
		
	}else if( flag==2  ){#NNDSVDar: fill in the zero elements with random values in the space :[0:average/100]
	
		ind1 <- W==0;
		ind2 <- H==0;
		n1 <- sum(ind1);
		n2 <- sum(ind2);
		
		average = mean(A);
		W[ind1] =  (average * runif(n1, min=0, max=1) / 100); 
		H[ind2] =  (average * runif(n2, min=0, max=1) / 100);
		
	}
	
	# return matrices W and H
	list(W=W, H=H)	
}

###########################################################################
# REGISTRATION
###########################################################################
setNMFSeed('nndsvd', .nndsvd.wrapper, overwrite=TRUE)

