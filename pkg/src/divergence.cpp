#ifndef NMF_DIVERGENCE_H // include header only once
#define NMF_DIVERGENCE_H

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>

extern "C" {

	SEXP divergence_update_H ( SEXP v, SEXP w, SEXP h, SEXP dup);
	SEXP divergence_update_W ( SEXP v, SEXP w, SEXP h, SEXP dup);

}

// include version for both double/integer storage.mode (defined as templates)
#include "divergence.cpp"

// define the exported versions (for SEXP)
SEXP divergence_update_H ( SEXP v, SEXP w, SEXP h, SEXP dup=ScalarLogical(1)) {

	if( TYPEOF(v) == REALSXP )
		return divergence_update_H(NUMERIC_POINTER(v), w, h
				, *LOGICAL(dup));
	else
		return divergence_update_H(INTEGER_POINTER(v), w, h
				, *LOGICAL(dup));
}

SEXP divergence_update_W ( SEXP v, SEXP w, SEXP h, SEXP dup=ScalarLogical(1)) {

	if( TYPEOF(v) == REALSXP )
		return divergence_update_W(NUMERIC_POINTER(v), w, h
				, *LOGICAL(dup));
	else
		return divergence_update_W(INTEGER_POINTER(v), w, h
				, *LOGICAL(dup));
}

#define NMF_DIVERGENCE_DONE
#else // END OF NMF_DIVERGENCE_H
#ifndef NMF_DIVERGENCE_DONE // START DEFINITION OF FUNCTIONS

/**
 * Divergence based multiplicative update for the mixture coefficients matrix H
 * from Brunet et al. algorithm.
 *
 * @param pV target matrix
 * @param w basis vector matrix
 * @param h mixture coefficient matrix to be updated
 * @param limInf limit inf for the updated entries. Applied only if limInf > 0
 * @param dup boolean (flag) that specifies if the update must be perform directly on w or
 * on a duplicated version of w
 *
 * @return the updated mixture coefficient matrix.
 *

  */
template <typename T_Rnumeric>
static SEXP divergence_update_H ( T_Rnumeric* pV, SEXP w, SEXP h, int dup=1)
{
	SEXP res;
	int nprotect = 0;

	// retrieve dimensions from W and H
	int n = INTEGER(GET_DIM(w))[0];
	int r = INTEGER(GET_DIM(w))[1];
	int p = INTEGER(GET_DIM(h))[1];

	// duplicate H (keeping attributes) or modify in place
	PROTECT(res = (dup != 0 ? duplicate(h) : h) ); nprotect++;

	// define internal pointers
	double* pW = NUMERIC_POINTER(w);
	double* pH = NUMERIC_POINTER(h);
	double* p_res = NUMERIC_POINTER(res);

	// allocate internal memory
	double* sumW = (double*) R_alloc(r, sizeof(double)); // will store column sums of W
	double* pWH = (double*) R_alloc(n, sizeof(double)); // will store the currently used column of WH

	// Compute update of H column by column
	for(int jH=0; jH < p; jH++){

		for (int iH=0; iH < r; iH++){ // compute value for H_ij

			// initialise values
			double tmp_res = 0.0;
			double &w_sum = sumW[iH];
			if( jH == 0 ) w_sum = 0.0;

			// compute cross-product w_.i by (v/wh)_.j
			for( int u=0; u<n; u++){

				// The jH-th column of WH is used to compute all elements of H_.j
				// => compute once and store the result for using for the next rows
				double wh_term = pWH[u];
				if( iH == 0 ){
					wh_term = 0.0;
					for (int k=0; k<r; k++){
						wh_term += pW[u + k*n] * pH[k + jH*r];
					}
					wh_term = pV[u + jH*n] / wh_term;
					pWH[u] = wh_term;
				}

				tmp_res +=  pW[u + iH*n] * wh_term;

				// compute sum of iH-th column of W (done only once)
				if( jH == 0 ) w_sum += pW[u + iH*n];
			}

			// multiplicative update
			p_res[iH + jH*r] = pH[iH + jH*r] * tmp_res / w_sum;
		}
	}

	// return result
	UNPROTECT(nprotect);
	return res;

}


/**
 * Divergence based multiplicative update for the basis matrix W
 * from Brunet et al. algorithm.
 *
 * @param pV target matrix
 * @param w basis vector matrix to be updated
 * @param h mixture coefficient matrix
 * @param limInf limit inf for the updated entries. Applied only if limInf > 0
 * @param dup boolean (flag) that specifies if the update must be perform directly on h or
 * on a duplicated version of h
 *
 * @return the updated basis vector matrix.
 *
 */

template <typename T_Rnumeric>
static SEXP divergence_update_W ( T_Rnumeric* pV, SEXP w, SEXP h, int dup=1)
{

	SEXP res;
	int nprotect = 0;

	// retrieve dimensions
	int n = INTEGER(GET_DIM(w))[0];
	int r = INTEGER(GET_DIM(w))[1];
	int p = INTEGER(GET_DIM(h))[1];

	// duplicate W (keeping attributes)
	PROTECT(res = (dup != 0 ? duplicate(w) : w) ); nprotect++;

	// define internal pointers
	double* pW = NUMERIC_POINTER(w);
	double* pH = NUMERIC_POINTER(h);
	double* p_res = NUMERIC_POINTER(res);

	// allocate internal memory
	double* sumH = (double*) R_alloc(r, sizeof(double)); // will store the row sums of H
	double* pWH = (double*) R_alloc(p, sizeof(double)); // will store currently used row of WH

	// Compute update of W row by row
	for(int iW=0; iW < n; iW++){

		for (int jW=0; jW < r; jW++){ // compute value for W_ij

			// initialise values
			double tmp_res = 0.0;
			double &h_sum = sumH[jW];
			if( iW == 0 ) h_sum = 0.0;

			// compute cross-product (v/wh)_i. by h_j.
			for( int u=0; u<p; u++){

				// The iW-th row of WH is used to compute all elements of W_i.
				// => compute once and store the result for using for the next columns
				if( jW == 0 ){
					double wh_term = 0.0;
					for (int k=0; k<r; k++){
						wh_term += pW[iW + k*n] * pH[k + u*r];
					}
					wh_term = pV[iW + u*n] / wh_term;
					pWH[u] = wh_term;
				}

				tmp_res +=  pH[jW + u*r] * pWH[u];

				// compute sum of jW-th row of H (done only once)
				if( iW == 0 ) h_sum += pH[jW + u*r];
			}

			// multiplicative update
			p_res[iW + jW*n] = pW[iW + jW*n] * tmp_res / h_sum;
		}

	}

	// return result
	UNPROTECT(nprotect);
	return res;
}

#endif //END ifndef NMF_DIVERGENCE_DONE
#endif //END ifdef NMF_DIVERGENCE_H
