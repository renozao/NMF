#ifndef NMF_EUCLIDEAN_H // include header only once
#define NMF_EUCLIDEAN_H

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>

extern "C" {

	// NMF from Lee and Seung (based on Euclidean norm)
	SEXP euclidean_update_H ( SEXP v, SEXP w, SEXP h, SEXP eps);
	SEXP euclidean_update_W ( SEXP v, SEXP w, SEXP h, SEXP eps);

	// NMF with offset
	SEXP offset_euclidean_update_H ( SEXP v, SEXP w, SEXP h, SEXP offset, SEXP eps);
	SEXP offset_euclidean_update_W ( SEXP v, SEXP w, SEXP h, SEXP offset, SEXP eps);

}

//////////////////////////////////
// STANDARD NMF: LEE
//////////////////////////////////

// include version for double storage.mode
#include "euclidean.cpp"

// include version for integer storage.mode
#define NMF_INT
#include "euclidean.cpp"
#undef NMF_INT

// define the exported versions (for SEXP)
SEXP euclidean_update_H ( SEXP v, SEXP w, SEXP h, SEXP eps){

	if( TYPEOF(v) == REALSXP )
		return euclidean_update_H(NUMERIC_POINTER(v), w, h, eps);
	else
		return euclidean_update_H(INTEGER_POINTER(v), w, h, eps);
}

SEXP euclidean_update_W ( SEXP v, SEXP w, SEXP h, SEXP eps){

	if( TYPEOF(v) == REALSXP )
		return euclidean_update_W(NUMERIC_POINTER(v), w, h, eps);
	else
		return euclidean_update_W(INTEGER_POINTER(v), w, h, eps);
}

//////////////////////////////////
// NMF WITH OFFSET
//////////////////////////////////
#define NMF_WITH_OFFSET
// include version for double storage.mode
#include "euclidean.cpp"

// include version for integer storage.mode
#define NMF_INT
#include "euclidean.cpp"
#undef NMF_INT
#undef NMF_WITH_OFFSET

// define the exported versions (for SEXP)
SEXP offset_euclidean_update_H ( SEXP v, SEXP w, SEXP h, SEXP offset, SEXP eps){

	if( TYPEOF(v) == REALSXP )
		return offset_euclidean_update_H(NUMERIC_POINTER(v), w, h, offset, eps);
	else
		return offset_euclidean_update_H(INTEGER_POINTER(v), w, h, offset, eps);
}

SEXP offset_euclidean_update_W ( SEXP v, SEXP w, SEXP h, SEXP offset, SEXP eps){

	if( TYPEOF(v) == REALSXP )
		return offset_euclidean_update_W(NUMERIC_POINTER(v), w, h, offset, eps);
	else
		return offset_euclidean_update_W(INTEGER_POINTER(v), w, h, offset, eps);
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
#ifdef NMF_WITH_OFFSET
SEXP offset_euclidean_update_H (
		#ifdef NMF_INT
		int*
		#else
		double*
		#endif
		pV
		, SEXP w, SEXP h, SEXP offset, SEXP eps)
#else
SEXP euclidean_update_H (
		#ifdef NMF_INT
		int*
		#else
		double*
		#endif
		pV
		, SEXP w, SEXP h, SEXP eps)
#endif
{
	SEXP res;
	int nprotect = 0;

	double eps_val = *NUMERIC_POINTER(eps);

	// retrieve dimensions
	int n = INTEGER(GET_DIM(w))[0];
	int r = INTEGER(GET_DIM(w))[1];
	int p = INTEGER(GET_DIM(h))[1];

	// duplicate H (keeping attributes)
	PROTECT(res = duplicate(h)); nprotect++;

	// define internal pointers
	//double* pV = NUMERIC_POINTER(v);
	double* pW = NUMERIC_POINTER(w);
	double* pH = NUMERIC_POINTER(h);
	double* p_res = NUMERIC_POINTER(res);

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
	double temp = 0;

	// Pre-compute symetric matrix t(W)W
	// -> allocate internal memory as a upper triangular in column major
	double* p_tWW = (double*) R_alloc((int) (r*(r+1))/2, sizeof(double));
	double* p_row = NULL;
	for( int i=r-1; i>=0; --i){
		p_row = pW + i*n;

		#ifdef NMF_WITH_OFFSET
		den_addon[i] = 0;
		#endif

		for( int j=r-1; j>=0; --j){
			temp = 0;
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
	// Compute update of H row by row
	for(int i=r-1; i>=0; --i){

		for (int j=p-1; j>=0; --j){ // compute value for H_ij

			// numerator
			double numerator = 0;
			for( int u=n-1; u>=0; --u)
				numerator += pW[u + i*n] * pV[u + j*n];

			double den = 0;
			for( int l=r-1; l>=0; --l){
				den += p_tWW[i > l ? ((i+1)*i)/2 + l : ((l+1)*l)/2 + i] * pH[l + j*r];
			}

			// add offset addon if necessary
			#ifdef NMF_WITH_OFFSET
			if( pOffset != NULL )
				den += den_addon[i];
			#endif

			// multiplicative update
			p_res[i + j*r] = ((temp = pH[i + j*r] * numerator) > eps_val ? temp : eps_val) / (den + eps_val);
		}
	}

	// return result
	UNPROTECT(nprotect);
	return res;
}

/**
 * Euclidean norm based multiplicative update for the basis matrix W
 * from Lee and Seung.
 * Also used in the NMF with Offset algorithm
 *
 * Note: for performance reason the dimension names are NOT conserved.
 */
#ifdef NMF_WITH_OFFSET
SEXP offset_euclidean_update_W (
		#ifdef NMF_INT
		int*
		#else
		double*
		#endif
		pV
		, SEXP w, SEXP h, SEXP offset, SEXP eps)
#else
SEXP euclidean_update_W (
		#ifdef NMF_INT
		int*
		#else
		double*
		#endif
		pV
		, SEXP w, SEXP h, SEXP eps)
#endif
{
	SEXP res;
	int nprotect = 0;

	double eps_val = *NUMERIC_POINTER(eps);

	// retrieve dimensions
	int n = INTEGER(GET_DIM(w))[0];
	int r = INTEGER(GET_DIM(w))[1];
	int p = INTEGER(GET_DIM(h))[1];

	// duplicate H (keeping attributes)
	PROTECT(res = duplicate(w)); nprotect++;

	// define internal pointers to data
	//double* pV = NUMERIC_POINTER(v);
	double* pW = NUMERIC_POINTER(w);
	double* pH = NUMERIC_POINTER(h);
	double* p_res = NUMERIC_POINTER(res);

	// extra variables in the case of an optional offset
	#ifdef NMF_WITH_OFFSET
	double *pOffset = NULL, *rowSumsH = NULL;

	if( offset != R_NilValue ){
		pOffset = NUMERIC_POINTER(offset);

		// pre-compute the row sums of H
		rowSumsH = (double*) R_alloc(r, sizeof(double));
		for( int i=r-1; i>=0; --i){
			rowSumsH[i] = 0;
			for( int j=p-1; j>=0; --j){
				rowSumsH[i] += pH[i + j*r];
			}
		}
	}
	#endif

	// auxiliary temporary variable
	double temp = 0;

	// Pre-compute symetric matrix Ht(H)
	// -> allocate internal memory as a lower triangular in column major
	double* p_HtH = (double*) R_alloc((int) (r*(r+1))/2, sizeof(double));
	for( int j=r-1; j>=0; --j){
		for( int i=j; i<r; ++i){
			temp = 0;
			for( int u=p-1; u>=0; --u){
				temp += pH[j + u*r] * pH[i + u*r];
			}

			p_HtH[((i+1)*i)/2 + j] = temp;
		}
	}

	// W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration
	// Compute update of H column by column
	double numerator = 0;
	double den = 0;
	for (int j=r-1; j>=0; --j){
		for(int i=n-1; i>=0; --i){// compute value for H_ij

			// numerator
			numerator = 0;
			for( int u=p-1; u>=0; --u)
				numerator += pV[i + u*n] * pH[j + u*r];

			den = 0;
			for( int l=r-1; l>=0; --l){
				// compute index of the stored value of t(w)w_[iH,l] in column major
				den += pW[i + l*n] * p_HtH[l < j ? ((j+1)*j)/2 + l : ((l+1)*l)/2 + j];
			}

			// add offset addon if necessary
			#ifdef NMF_WITH_OFFSET
			if( pOffset != NULL )
				den += pOffset[i] * rowSumsH[j];
			#endif

			// multiplicative update
			p_res[i + j*n] = ((temp = pW[i + j*n] * numerator) > eps_val ? temp : eps_val) / (den + eps_val);
		}
	}

	// return result
	UNPROTECT(nprotect);
	return res;
}

#endif //END ifndef NMF_EUCLIDEAN_DONE
#endif //END ifdef NMF_EUCLIDEAN_H
