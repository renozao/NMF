# TODO: Add comment
# 
# Author: Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###############################################################################


############################################
# OPTIMIZATION WITH C/C++
############################################

library(inline)

file_get_content <- function(f){ paste(readLines(f), collapse="\n") }

lee.optim.H <- cfunction(signature(v="matrix", w="matrix", h="matrix", eps="numeric"), file_get_content('lee.cpp'))
lee.optim.W <- cfunction(signature(v="matrix", w="matrix", h="matrix", eps="numeric"), file_get_content('lee.W.cpp'))

off.optim.H <- cfunction(signature(v="matrix", w="matrix", h="matrix", eps="numeric", offset='numeric'), file_get_content('lee.cpp'), otherdefs='#define WITH_OFFSET')
off.optim.W <- cfunction(signature(v="matrix", w="matrix", h="matrix", eps="numeric", offset='numeric'), file_get_content('lee.W.cpp'), otherdefs='#define WITH_OFFSET')

div.optim <- function(i, v, data, ...)
{
	if(i==1) print("use optimized")
	# retrieve each factor
	w <- basis(data); h <- coef(data);
	
	# standard divergence-reducing NMF update for H
	h <- funxH(v, w, h)
	
	# standard divergence-reducing NMF update for W
	w <- funxW(v, w, h)
	
	#every 10 iterations: adjust small values to avoid underflow 
	if( i %% 10 == 0 ){
		#precision threshold for numerical stability
		eps <- .Machine$double.eps
		h[h<eps] <- eps;
		w[w<eps] <- eps;
	}
		
	#return the modified data
	basis(data) <- w; coef(data) <- h;	
	return(data)
	
}

speed.l <- function(){
	
	cat("H\n")
	lee.optim.H(target, W0, H0, eps)
	body(lee.optim.H)[2:3] <- NULL
	print(system.time({sapply(seq(1030), function(x){ a <- std.euclidean.update.h(target, W0, H0, eps=eps); NULL}) }))
	print(system.time({ sapply(seq(1030), function(x){ a <- lee.optim.H(target, W0, H0, eps); NULL}) }))
	
	cat("W\n")
	lee.optim.W(target, W0, H0, eps)
	body(lee.optim.W)[2:3] <- NULL
	print(system.time({sapply(seq(1030), function(x){ a <- std.euclidean.update.w(target, W0, H0, eps=eps); NULL}) }))
	print(system.time({ sapply(seq(1030), function(x){ a <- lee.optim.W(target, W0, H0, eps); NULL}) }))
	
	cat("NMF\n")
	lee.optim <- function(i, v, data, rescale=TRUE, ...)
	{
		if(i==1) print("use optimized")
		# retrieve each factor
		w <- basis(data); h <- coef(data);	
		
		#precision threshold for numerical stability
		eps <- 10^-9
		
		# compute the estimate WH
		#wh <- estimate(data)
		
		# euclidean-reducing NMF iterations	
		# H_au = H_au (W^T V)_au / (W^T W H)_au
		#h <- pmax(h * (t(w) %*% v),eps) / ((t(w) %*% w) %*% h + eps);
		h <- lee.optim.H(v, w, h, eps=eps)
		
		# update H and recompute the estimate WH
		#metaprofiles(data) <- h
		#wh <- estimate(data)
		
		# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
		#w <- pmax(w * (v %*% t(h)), eps) / (w %*% (h %*% t(h)) + eps);
		w <- lee.optim.W(v, w, h, eps=eps)
		#rescale columns TODO: effect of rescaling? the rescaling makes the update with offset fail
		if( rescale ) w <- sweep(w, 2L, colSums(w), "/", check.margin=FALSE)
		
		#return the modified data
		basis(data) <- w; coef(data) <- h;	
		return(data)
	}
	
	lo <- nmfAlgorithm('lee');
	lo@Update <- lee.optim;
	name(lo) <- 'lee.optim'

	nmf(target, 3, seed=123456, method='lee')
	print(runtime(res))
	res <- nmf(target, 3, seed=123456, method=lo)
	print(runtime(res))
	
}

check.l <- function(){
	eps <- 10^-9
	
	# test H
	resRefH <- std.euclidean.update.h(target, W0, H0, eps=eps)
	resH <- lee.optim.H(target, W0, H0, eps)
	cat('H ', all.equal(resH, resRefH), "\n")
	
	# test W
	resRefW <- std.euclidean.update.w(target, W0, H0, eps=eps)
	resW <- lee.optim.W(target, W0, H0, eps)
	cat('W ', all.equal(resW, resRefW), "\n")	
}


source('package.R')

if( FALSE ){

source('../package/R/optim-C.R', chd=T)
target <- as.matrix(read.table('target.txt'))
storage.mode(target) <- 'double'
W0 <- as.matrix(read.table('W.txt'))
H0 <- as.matrix(read.table('H.txt'))
s <- nmfModelW=W0, H=H0)
eps <- 10^-9


# test H
resRefH <- std.divergence.update.h(target, W0, H0)
resH <- funxH(target, W0, H0)
all.equal(resH, resRefH)

# test W
resRefW <- std.divergence.update.w(target, W0, H0)
resW <- funxW(target, W0, H0)
all.equal(resW, resRefW)


body(funxH)[2:3] <- NULL
system.time({sapply(seq(670), function(x){ h <- std.divergence.update.h(target, W0, H0); NULL}) })
system.time({ sapply(seq(670), function(x){h <- funxH(target, W0, H0); NULL}) })

body(funxW)[2:3] <- NULL
system.time({sapply(seq(670), function(x){ w <- std.divergence.update.w(target, W0, H0); NULL}) })
system.time({ sapply(seq(670), function(x){w <- funxW(target, W0, H0); NULL}) })


nmf(target, 3, seed=123456, method=bo)
nmf(target, 3, seed=123456, method='brunet')

#http://subclipse.tigris.org/update_1.6.x

#################
# LEE
#################

# test H
resRefH <- std.euclidean.update.h(target, W0, H0, eps=eps)
resH <- lee.optim.H(target, W0, H0, eps)
all.equal(resH, resRefH)
body(lee.optim.H)[2:3] <- NULL
system.time({sapply(seq(1030), function(x){ h <- std.euclidean.update.h(target, W0, H0, eps=eps); NULL}) })
system.time({ sapply(seq(1030), function(x){h <- lee.optim.H(target, W0, H0, eps); NULL}) })

# test W
resRefW <- std.euclidean.update.w(target, W0, H0, eps=eps)
resW <- lee.optim.W(target, W0, H0, eps)
all.equal(resW, resRefW)
body(lee.optim.W)[2:3] <- NULL
system.time({sapply(seq(1030), function(x){ h <- std.euclidean.update.w(target, W0, H0, eps=eps); NULL}) })
system.time({ sapply(seq(1030), function(x){h <- lee.optim.W(target, W0, H0, eps); NULL}) })

# Lee Optimized	
lo <- new('NMFStrategyIterative', name='olee', objective='euclidean'
		, Update='CO_nmf.update.lee'
		, Stop='connectivity'
)


reso <- nmf(target, 3, seed=123456, method=lo)
res <- nmf(target, 3, seed=123456, method='lee')
all.equal(basis(res), basis(reso))
all.equal(coef(res), coef(reso))
runtime(reso)/runtime(res())

##################
# OFFSET
##################
off <- rep(1, nrow(W0))
so <- nmfModelmodel='NMFOffset', W=W0, H=H0, offset=off)
resRefH <- std.euclidean.update.h(target, W0, H0, wh=W0%*%H0 + off, eps=eps)
resH <- off.optim.H(target, W0, H0, eps, off)
all.equal(resH, resRefH)
body(off.optim.H)[2:3] <- NULL
n <- 760
t <- system.time({sapply(seq(n), function(x){ h <- std.euclidean.update.h(target, W0, H0, wh=W0%*%H0 + off, eps=eps); NULL}) });t;t/n
t <- system.time({ sapply(seq(n), function(x){h <- off.optim.H(target, W0, H0, eps, off); NULL}) }); t;t/n


resRefW <- std.euclidean.update.w(target, W0, H0, wh=W0%*%H0 + off, eps=eps)
resW <- off.optim.W(target, W0, H0, eps, off)
all.equal(resW, resRefW)
body(off.optim.W)[2:3] <- NULL

# NMF with offset Optimized
oo <- new('NMFStrategyIterative', name='ooffset', objective='euclidean'
		, model = 'NMFOffset'
		, Update='CO_nmf.update.offset'
		, Stop='connectivity'
)


resO <- nmf(target, seed=so, method=oo, maxIter=1)
res <- nmf(target,  seed=so, method='off', maxIter=1)
res <- nmf(target, seed=so, method='off')


}
