#ifndef NMF_EUCLIDEAN_H // include header only once
#define NMF_EUCLIDEAN_H

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>

extern "C" {

	// NMF from Lee and Seung (based on Euclidean norm)
	SEXP euclidean_update_H ( SEXP v, SEXP w, SEXP h, SEXP eps, SEXP nbterms, SEXP ncterms, SEXP dup);
	SEXP euclidean_update_W ( SEXP v, SEXP w, SEXP h, SEXP eps, SEXP weight, SEXP nbterms, SEXP ncterms, SEXP dup);

	// NMF with offset
	SEXP offset_euclidean_update_H ( SEXP v, SEXP w, SEXP h, SEXP offset, SEXP eps, SEXP dup);
	SEXP offset_euclidean_update_W ( SEXP v, SEXP w, SEXP h, SEXP offset, SEXP eps, SEXP dup);

}

//////////////////////////////////
// STANDARD NMF: LEE
//////////////////////////////////

// include version for both double/integer storage.mode (defined as templates)
#include "euclidean.cpp"

// define the exported versions (for SEXP)
SEXP euclidean_update_H ( SEXP v, SEXP w, SEXP h, SEXP eps
		, SEXP nbterms=ScalarInteger(0), SEXP ncterms=ScalarInteger(0)
		, SEXP dup=ScalarLogical(1)){

	if( TYPEOF(v) == REALSXP ){
		return euclidean_update_H(NUMERIC_POINTER(v), w, h, eps
				, *INTEGER(nbterms), *INTEGER(ncterms)
				, *LOGICAL(dup));
	}else{
		return euclidean_update_H(INTEGER_POINTER(v), w, h, eps
				, *INTEGER(nbterms), *INTEGER(ncterms)
				, *LOGICAL(dup));
	}
}

//////////////////////////////////
// NMF WITH WEIGHT
//////////////////////////////////
#define NMF_WITH_WEIGHT
#include "euclidean.cpp"
#undef NMF_WITH_WEIGHT

SEXP euclidean_update_W ( SEXP v, SEXP w, SEXP h, SEXP eps
		, SEXP weight = R_NilValue
		, SEXP nbterms=ScalarInteger(0), SEXP ncterms=ScalarInteger(0)
		, SEXP dup=ScalarLogical(1)){

	int nb = *INTEGER(nbterms), nc = *INTEGER(ncterms);
	bool copy = *LOGICAL(dup);
	if( TYPEOF(v) == REALSXP ){
		if( isNull(weight) ){
			return euclidean_update_W(NUMERIC_POINTER(v), w, h, eps, nb, nc, copy);
		}else{
			return weuclidean_update_W(NUMERIC_POINTER(v), w, h, eps, weight, nb, nc, copy);
		}
	}else{
		if( isNull(weight) ){
			return euclidean_update_W(INTEGER_POINTER(v), w, h, eps, nb, nc, copy);
		}else{
			return weuclidean_update_W(INTEGER_POINTER(v), w, h, eps, weight, nb, nc, copy);
		}
	}
}

//////////////////////////////////
// NMF WITH OFFSET
//////////////////////////////////
#define NMF_WITH_OFFSET

// include version for both double/integer storage.mode (defined as templates)
#include "euclidean.cpp"

#undef NMF_WITH_OFFSET

// define the exported versions (for SEXP)
SEXP offset_euclidean_update_H ( SEXP v, SEXP w, SEXP h, SEXP offset, SEXP eps, SEXP dup=ScalarLogical(1)){

	if( TYPEOF(v) == REALSXP )
		return offset_euclidean_update_H(NUMERIC_POINTER(v), w, h, offset, eps, *LOGICAL(dup));
	else
		return offset_euclidean_update_H(INTEGER_POINTER(v), w, h, offset, eps, *LOGICAL(dup));
}

SEXP offset_euclidean_update_W ( SEXP v, SEXP w, SEXP h, SEXP offset, SEXP eps, SEXP dup=ScalarLogical(1)){

	if( TYPEOF(v) == REALSXP )
		return offset_euclidean_update_W(NUMERIC_POINTER(v), w, h, offset, eps, *LOGICAL(dup));
	else
		return offset_euclidean_update_W(INTEGER_POINTER(v), w, h, offset, eps, *LOGICAL(dup));
}

#define NMF_EUCLIDEAN_DONE
#else // END OF NMF_EUCLIDEAN_H
#ifndef NMF_EUCLIDEAN_DONE // START DEFINITION OF FUNCTIONS

/**
 * Euclidean norm based multiplicative update for the mixture coefficients matrix H
 * from Lee and Seung.
 * Also used in the NMF with Offset algorithm
 *
 * Note: for performance reason the dimension names are NOT conserved.
 */
#ifndef NMF_WITH_WEIGHT

template <typename T_Rnumeric>
static
#ifdef NMF_WITH_OFFSET
SEXP offset_euclidean_update_H (
#else
SEXP euclidean_update_H (
#endif
		T_Rnumeric* pV, SEXP w, SEXP h
#ifdef NMF_WITH_OFFSET
		, SEXP offset
#endif
		, SEXP eps
#ifndef NMF_WITH_OFFSET
		, int nbterms=0, int ncterms=0
#endif
		, int dup=1)
{
	SEXP res;
	int nprotect = 0;

	double eps_val = *NUMERIC_POINTER(eps);

	// retrieve dimensions
	int n = INTEGER(GET_DIM(w))[0];
	int r = INTEGER(GET_DIM(w))[1];
	int p = INTEGER(GET_DIM(h))[1];
	// get number of non-fixed terms
	int vr =
#ifdef NMF_WITH_OFFSET
	r;
#else
	r - ncterms;
#endif

	// duplicate H (keeping attributes)
	PROTECT( res = (dup != 0 ? duplicate(h) : h) ); nprotect++;

	// define internal pointers
	double* pW = NUMERIC_POINTER(w);
	double* pH = NUMERIC_POINTER(h);
	double* p_res = NUMERIC_POINTER(res);
	double* pH_buffer = (double*) R_alloc(r, sizeof(double));

	// extra variables in the case of an optional offset
	#ifdef NMF_WITH_OFFSET
	double *pOffset = NULL, *den_addon = NULL;

	if( offset != R_NilValue ){
		pOffset = NUMERIC_POINTER(offset);
		den_addon = (double*) R_alloc(r, sizeof(double));
		//memset(den_addon, 0, r);
	}
	#endif

	// auxiliary temporary variable
	double temp = 0.0;

	// Pre-compute symmetric matrix t(W)W
	// -> allocate internal memory as a upper triangular in column major
	double* p_tWW = (double*) R_alloc((int) (r*(r+1))/2, sizeof(double));
	double* p_row = NULL;
	for( int i=r-1; i>=0; --i){
		p_row = pW + i*n;

		#ifdef NMF_WITH_OFFSET
		den_addon[i] = 0.0;
		#endif

		for( int j=r-1; j>=0; --j){
			temp = 0.0;
			for( int u=n-1; u>=0; --u){
				temp += p_row[u] * pW[u + j*n];

				#ifdef NMF_WITH_OFFSET
				if( pOffset != NULL && j==0 )
					den_addon[i] += p_row[u] * pOffset[u];
				#endif
			}

			p_tWW[((j+1)*j)/2 + i] = temp;
		}
	}

	// H_au = H_au (W^T V)_au / (W^T W H)_au
	// Compute update of H column by column
	for (int j=p-1; j>=0; --j){

		for(int i=vr-1; i>=0; --i){ // compute value for H_ij (only non-fixed entries)

			// numerator
			double numerator = 0.0;
			for( int u=n-1; u>=0; --u)
				numerator += pW[u + i*n] * pV[u + j*n];

			double den = 0.0;
			for( int l=r-1; l>=0; --l){ // use all entries (fixed and non-fixed)
				// bufferize jth-column of H, as it can be changed at the end of the current i-loop
				if( i==vr-1 )
					pH_buffer[l] = pH[l + j*r];
				den += p_tWW[i > l ? ((i+1)*i)/2 + l : ((l+1)*l)/2 + i] * pH_buffer[l];
			}

			// add offset addon if necessary
			#ifdef NMF_WITH_OFFSET
			if( pOffset != NULL )
				den += den_addon[i];
			#endif

			// multiplicative update
			p_res[i + j*r] = ((temp = pH_buffer[i] * numerator) > eps_val ? temp : eps_val) / (den + eps_val);
		}
	}

	// return result
	UNPROTECT(nprotect);
	return res;
}
#endif

/**
 * Euclidean norm based multiplicative update for the basis matrix W
 * from Lee and Seung.
 * Also used in the NMF with Offset algorithm
 *
 * Note: for performance reason the dimension names are NOT conserved.
 */
template <typename T_Rnumeric>
static SEXP
#ifdef NMF_WITH_OFFSET
	offset_euclidean_update_W
#else
#ifdef NMF_WITH_WEIGHT
	weuclidean_update_W
#else
	euclidean_update_W
#endif
#endif
	(T_Rnumeric* pV, SEXP w, SEXP h
#ifdef NMF_WITH_OFFSET
		, SEXP offset
#endif
		, SEXP eps
#ifdef NMF_WITH_WEIGHT
		, SEXP weight
#endif
#ifndef NMF_WITH_OFFSET
		, int nbterms=0, int ncterms=0
#endif
		, int dup=1)
{
	SEXP res;
	int nprotect = 0;

	// setup variables for enforcing a limit Inf on the entries
	double limInf = *NUMERIC_POINTER(eps);

	// retrieve dimensions
	int n = INTEGER(GET_DIM(w))[0];
	int r = INTEGER(GET_DIM(w))[1];
	int p = INTEGER(GET_DIM(h))[1];

	// duplicate H (keeping attributes)
	//PROTECT(res = duplicate(w)); nprotect++;
	PROTECT(res = (dup != 0 ? duplicate(w) : w) ); nprotect++;

	// define internal pointers to data
	double* pW = NUMERIC_POINTER(w);
	double* pH = NUMERIC_POINTER(h);
	double* p_res = NUMERIC_POINTER(res);
	double* pW_buffer = (double*) R_alloc(r, sizeof(double));

	// extra variables in the case of an optional offset
	#ifdef NMF_WITH_OFFSET
	double *pOffset = NULL, *rowSumsH = NULL;

	if( offset != R_NilValue ){
		pOffset = NUMERIC_POINTER(offset);

		// pre-compute the row sums of H
		rowSumsH = (double*) R_alloc(r, sizeof(double));
		for( int i=r-1; i>=0; --i){
			rowSumsH[i] = 0.0;
			for( int j=p-1; j>=0; --j){
				rowSumsH[i] += pH[i + j*r];
			}
		}
	}
	#endif

#ifdef NMF_WITH_WEIGHT
	// take sample weights into account
	double* p_weight = !isNull(weight) ? NUMERIC_POINTER(weight) : NULL;
	double beta = -1.0;
	if( p_weight == NULL ){// <=> no weights
		beta = 1.0;
	}
	else if( length(weight) == 1 ){// all weighted are the same
		// NB: theoretically this is equivalent to weight=1, but may be used
		// to test it in practice (with the numerical adjustments via eps)
		beta = *p_weight;
	}
	// fill weight vector with single value
	if( beta > 0 ){
		double* pw = p_weight = (double*) R_alloc(p, sizeof(double));
		for(int i=0; i<p; ++i, ++pw){
			*pw = beta;
		}
	}
#endif

	// auxiliary temporary variable
	double temp = 0.0;

	// Pre-compute symetric matrix Ht(H)
	// -> allocate internal memory as a lower triangular in column major
	double* p_HtH = (double*) R_alloc((int) (r*(r+1))/2, sizeof(double));
	for( int j=r-1; j>=0; --j){
		for( int i=j; i<r; ++i){
			temp = 0.0;
			for( int u=p-1; u>=0; --u){
				temp += pH[j + u*r] * pH[i + u*r]
#ifdef NMF_WITH_WEIGHT
				                    * p_weight[u]
#endif
				;
			}

			p_HtH[((i+1)*i)/2 + j] = temp;
		}
	}

	// W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration
	// Compute update of W row by row
	double numerator = 0.0;
	double den = 0.0;
	for(int i=n-1; i>=0; --i){

		for (int j=r-1; j>=0; --j){// compute value for W_ij

			// numerator
			numerator = 0.0;
			for( int u=p-1; u>=0; --u){
				numerator += pV[i + u*n] * pH[j + u*r]
#ifdef NMF_WITH_WEIGHT
				                         * p_weight[u]
#endif
				;
			}

			den = 0.0;
			for( int l=r-1; l>=0; --l){
				// bufferize ith-row of W, as it can be changed at the end of the current j-loop
				if( j==r-1 )
					pW_buffer[l] = pW[i + l*n];
				// compute index of the stored value of t(w)w_[iH,l] in column major
				den += pW_buffer[l] * p_HtH[l < j ? ((j+1)*j)/2 + l : ((l+1)*l)/2 + j];
			}

			// add offset addon if necessary
			#ifdef NMF_WITH_OFFSET
			if( pOffset != NULL )
				den += pOffset[i] * rowSumsH[j];
			#endif

			// multiplicative update
			temp = pW_buffer[j] * numerator;
			p_res[i + j*n] = ( temp < limInf ? limInf : temp ) / (den + limInf);
		}
	}

	// return result
	UNPROTECT(nprotect);
	return res;
}

#endif //END ifndef NMF_EUCLIDEAN_DONE
#endif //END ifdef NMF_EUCLIDEAN_H
