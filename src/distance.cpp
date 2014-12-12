#ifndef NMF_DISTANCE_H // include header only once
#define NMF_DISTANCE_H
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>

extern "C" {

	SEXP Euclidean_rss( SEXP x, SEXP y);
	SEXP KL_divergence( SEXP x, SEXP y);

	// compute silhouette on large dataset
	SEXP big_silhouette(SEXP x, SEXP N, SEXP idx, SEXP method);

}

//define helper macro
#define both_non_NA(a,b) (!ISNAN(a) && !ISNAN(b))

// include versions for double-double storage.mode
#include "distance.cpp"
// include versions for double-integer storage.mode
#define NMF_ARG2_INT
#include "distance.cpp"
#undef NMF_ARG2_INT

// include versions for integer-* storage.mode
#define NMF_ARG1_INT
	#include "distance.cpp"

	// include versions for integer-integer storage.mode
	#define NMF_ARG2_INT
	#include "distance.cpp"
	#undef NMF_ARG2_INT
#undef NMF_ARG1_INT

// define the exported version of RSS (for SEXP)
SEXP Euclidean_rss ( SEXP x, SEXP y){

	// retrieve dimensions
	int n = INTEGER(GET_DIM(x))[0];
	int p = INTEGER(GET_DIM(x))[1];

    if( INTEGER(GET_DIM(y))[0] != n )
    	error("non-conformable arrays (rows)");
    if( INTEGER(GET_DIM(y))[1] != p )
    	error("non-conformable arrays (columns)");

	if( TYPEOF(x) == REALSXP ){// x is double
		if( TYPEOF(y) == REALSXP )// x and y are double
			return rss( NUMERIC_POINTER(x), NUMERIC_POINTER(y), n, p);
		else// x is double, y is integer
			return rss( NUMERIC_POINTER(x), INTEGER_POINTER(y), n, p);
	}else{
		if( TYPEOF(y) == REALSXP ) // x is integer, y is double
			return rss( INTEGER_POINTER(x), NUMERIC_POINTER(y), n, p);
		else // x is integer, y is integer
			return rss( INTEGER_POINTER(x), INTEGER_POINTER(y), n, p);
	}

}

// define the exported version of KL (for SEXP)
SEXP KL_divergence ( SEXP x, SEXP y){

	// retrieve dimensions
	int n = INTEGER(GET_DIM(x))[0];
	int p = INTEGER(GET_DIM(x))[1];

    if( INTEGER(GET_DIM(y))[0] != n )
    	error("non-conformable arrays (rows)");
    if( INTEGER(GET_DIM(y))[1] != p )
    	error("non-conformable arrays (columns)");

	if( TYPEOF(x) == REALSXP ){// x is double
		if( TYPEOF(y) == REALSXP )// x and y are double
			return KL( NUMERIC_POINTER(x), NUMERIC_POINTER(y), n, p);
		else// x is double, y is integer
			return KL( NUMERIC_POINTER(x), INTEGER_POINTER(y), n, p);
	}else{
		if( TYPEOF(y) == REALSXP ) // x is integer, y is double
			return KL( INTEGER_POINTER(x), NUMERIC_POINTER(y), n, p);
		else // x is integer, y is integer
			return KL( INTEGER_POINTER(x), INTEGER_POINTER(y), n, p);
	}

}

/**
 * Computes silhouette on Large Data
 *
 * @param x input centered and normalised matrix with samples in columns
 * @param i index of each cluster
 * @param method distance method
 */
SEXP big_silhouette(SEXP x, SEXP N, SEXP idx, SEXP method = R_NilValue){

	// process arguments
	int n = INTEGER(GET_DIM(x))[0];
	int p = INTEGER(GET_DIM(x))[1];
	int n_k = *INTEGER_POINTER(N);
	int n_item = p; int l_item = n;
	const double* m_x = NUMERIC_POINTER(x);
	const int* p_idx = INTEGER_POINTER(idx);

	// allocate result object: matrix of intra and inter cluster distance
	int nprotect = 0;
	SEXP res = PROTECT(allocMatrix(REALSXP, n_item, n_k)); nprotect++;
	double* p_res = NUMERIC_POINTER(res);
	double* ptr = p_res;
	for(int i=n_item * n_k; i>0; --i) *ptr++ = 0;

	// loop over item
	for(int j=0; j<n_item; ++j){
		const double* p_item = m_x + j * l_item;

		// compute distance to all other items
		for(int l=(j+1); l<n_item; ++l){

			// compute distance to item
			long double d = 0;
			const double* p_a = p_item;
			const double* p_b = m_x + l * l_item;
			for(int i=0; i<l_item; i++){
				d += *p_a++ * *p_b++;
			}
			d = 1 - d / (l_item - 1);

			// add to distance matrix
			int a_k = p_idx[j] - 1;
			int b_k = p_idx[l] - 1;
			p_res[j + b_k * n_item] += d;
			p_res[l + a_k * n_item] += d;
		}
	}

	UNPROTECT(nprotect);
	return res;


}

#define NMF_DISTANCE_DONE
#else // END OF NMF_DISTANCE_H
#ifndef NMF_DISTANCE_DONE // START DEFINITION OF FUNCTIONS

// Internal function that computes the RSS
SEXP rss(
	#ifdef NMF_ARG1_INT
	int*
	#else
	double*
	#endif
	px,
	#ifdef NMF_ARG2_INT
	int*
	#else
	double*
	#endif
	py
	, int n, int p){

    double dev=0, dist=0;
    double xval, yval;
    //int count = 0;
    for(int i=n-1; i>=0; --i) {
    	for(int j=p-1; j>=0; --j) {
    		xval = px[i + j*n]; yval = py[i + j*n];
			if (both_non_NA(xval, yval)) {
				dev = xval - yval;
				if (!ISNAN(dev)) {
					dist += dev * dev;
					//count++;
				}
				else return ScalarReal(NA_REAL);
			}
			else return ScalarReal(NA_REAL);
		}
    }
    //if (count == 0) return ScalarReal(NA_REAL);
    return ScalarReal(dist);
}

// Internal function that computes the KL divergence
SEXP KL(
	#ifdef NMF_ARG1_INT
	int*
	#else
	double*
	#endif
	px,
	#ifdef NMF_ARG2_INT
	int*
	#else
	double*
	#endif
	py
	, int n, int p){

    double dev=0, dist=0;
    double xval, yval;
    for(int i=n-1; i>=0; --i) {
    	for(int j=p-1; j>=0; --j) {
    		xval = px[i + j*n]; yval = py[i + j*n];
    		if( xval == 0 )
    			dev = yval;
    		else if (both_non_NA(xval, yval))
				dev = xval * log((double) xval / yval) - xval + yval;
			else return ScalarReal(NA_REAL);

    		// only add and continue if the term is not NA
    		if ( R_FINITE(dev) )
    			dist += dev;
    		else
    			return ScalarReal(dev);
		}
    }
    return ScalarReal(dist);
}

#endif //END ifndef NMF_DISTANCE_DONE
#endif //END ifdef NMF_DISTANCE_H

