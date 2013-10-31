# Implementation of siNMF from Badea (2008)
# 
# Author: Renaud Gaujoux
# Creation: 09 Jul 2012
###############################################################################

#' @include registry-algorithms.R
NULL

siNMF_R <- function(i, v, data, beta0=1, scale=TRUE, ...){
	
    # retrieve each factor
    w <- basis(data); h <- coef(data);
	# fixed terms
	nb <- nbterms(data); nc <- ncterms(data)
    
	if( i == 1 ){
	    if( !nc )
	  		stop("Method 'siNMF' requires a formula based model")
		
		if( !is.na(beta0) ){
			# compute beta
		    gr <- as.numeric(cterms(data)[[1L]])
		    beta <- beta0 * (norm(v[,gr==1], 'F') / norm(v[,gr==2], 'F'))^2
	    	# make sweeping vector
			vbeta <- rep(1, ncol(v))
			vbeta[gr==2] <- sqrt(beta)
		    staticVar('vbeta', vbeta, init=TRUE)
			# sweep data
			staticVar('v', sweep(v, 2L, vbeta, '*', check.margin=FALSE), init=TRUE)
		}
		# store non-fixed coef indexes
		staticVar('icoef', icoef(data), init=TRUE)
    }
    
    #precision threshold for numerical stability
    eps <- 10^-9
    
	sh <- h 
	if( !is.na(beta0) ){
		# retrieved swept matrix
		sv <- staticVar('v')
		vbeta <- staticVar('vbeta')
		# sweep h with beta
		sh <- sweep(h, 2L, vbeta, '*', check.margin=FALSE)
	}
	
    # compute standard euclidean updates
	w <- nmf_update.euclidean.w(sv, w, sh, eps=eps, nbterms=nb, ncterms=nc, copy=TRUE)
    h <- nmf_update.euclidean.h(v, w, h, eps=eps, nbterms=nb, ncterms=nc, copy=TRUE)
	
	# normalize columns of w
	if( scale ){
		icoef <- staticVar('icoef')
		wb <- w[, icoef]
		d <- sqrt(colSums(wb^2))
	    w[, icoef] <- sweep(wb, 2L, d, '/') 
		h[icoef, ] <- sweep(h[icoef, ], 1L, d, '*')
	}
	
	.coef(data) <- h
	.basis(data) <- w
    data
}

setNMFMethod('.siNMF', 'lee', Update=siNMF_R)

siNMF <- function(i, v, data, beta0=1, scale=TRUE, eps=10^-9, ...){
	
    # retrieve each factor
    w <- basis(data); h <- coef(data);
	# fixed terms
	nb <- nbterms(data); nc <- ncterms(data)
    
	if( i == 1 ){
	    if( !nc )
	  		stop("Method 'siNMF' requires a formula based model")
		
		vbeta <- NULL
		if( !is.na(beta0) ){
			# compute beta
		    gr <- cterms(data)[[1L]]
			gr <- droplevels(gr)
			# make sweeping vector
			vbeta <- rep(1, ncol(v))
			idx <- split(1:ncol(v), gr)
			# compute base value from first level
			beta <- beta0 * norm(v[,idx[[1]]], 'F')^2
			vbeta <- lapply(idx[-1], function(j){
				rep(beta / norm(v[,j], 'F')^2, length(j))
			})
			vbeta <- c(rep(1, length(idx[[1]])), unlist(vbeta, use.names=FALSE))
			vbeta <- vbeta[order(unlist(idx))]
		}
		# store weights
		staticVar('beta', vbeta, init=TRUE)
		# store non-fixed coef indexes
		staticVar('icoef', icoef(data), init=TRUE)
    }
    
    # retrieve weights
	beta <- staticVar('beta')
	
    # compute standard euclidean updates
	w <- nmf_update.euclidean.w(v, w, h, eps=eps, weight=beta, nbterms=nb, ncterms=nc, copy=FALSE)
    h <- nmf_update.euclidean.h(v, w, h, eps=eps, nbterms=nb, ncterms=nc, copy=FALSE)
	
	# normalize columns of w
	if( scale ){
		icoef <- staticVar('icoef')
		wb <- w[, icoef]
		d <- sqrt(colSums(wb^2))
	    w[, icoef] <- sweep(wb, 2L, d, '/', check.margin=FALSE) 
		h[icoef, ] <- sweep(h[icoef, ], 1L, d, '*', check.margin=FALSE)
	}
	
	.coef(data) <- h
	.basis(data) <- w
    data
}

nmfAlgorithm.siNMF <- setNMFMethod('siNMF', 'lee', Update=siNMF)
