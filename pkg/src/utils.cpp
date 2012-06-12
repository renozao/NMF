#ifndef NMF_UTILS_H // include header only once
#define NMF_UTILS_H

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>

extern "C" {

	/** Returns the pointer address of 'x' as a character string*/
	SEXP ptr_address (SEXP x);

	/** Clone an object 'x'*/
	SEXP clone_object (SEXP x);

	/** pmin in place with 'y' being a single numeric value*/
	SEXP ptr_pmin (SEXP x, SEXP y, SEXP skip);

	/** Apply inequality constraints in place. */
	SEXP ptr_neq_constraints(SEXP x, SEXP constraints, SEXP ratio=R_NilValue, SEXP value=R_NilValue);

	/** Minimum per column*/
	SEXP colMin(SEXP x);

	/** Maximum per row*/
	SEXP colMax(SEXP x);

	/** Test if an external pointer is null.
	 *
	 * Function taken from the package bigmemory (v4.2.11).
	 */
	SEXP ptr_isnil(SEXP address)
	{
	  void *ptr = R_ExternalPtrAddr(address);
	  SEXP ret = PROTECT(NEW_LOGICAL(1));
	  LOGICAL_DATA(ret)[0] = (ptr==NULL) ? (Rboolean)TRUE : Rboolean(FALSE);
	  UNPROTECT(1);
	  return(ret);
	}

}

// define the exported versions (for SEXP)
SEXP ptr_address (SEXP x){

	SEXP ans = R_NilValue;
	char tmp[15];
	PROTECT(ans = allocVector(STRSXP, 1));
	sprintf(tmp, "%p", (void *) x);
	SET_STRING_ELT(ans, 0, mkChar(tmp));
	UNPROTECT(1);
	return ans;
}

SEXP clone_object (SEXP x){

	return Rf_duplicate(x);

}

SEXP ptr_pmin(SEXP x, SEXP y, SEXP skip=R_NilValue){

	int n = length(x);
	double* p_x = ( isNull(x) ? NULL : NUMERIC_POINTER(x) );
	double lim = isNull(y) ? -1.0 : *NUMERIC_POINTER(y);

	// backup skipped values
	int n_skip = length(skip);
	int ncol = isNull(GET_DIM(x)) ? 1 : INTEGER(GET_DIM(x))[1];
	int nrow = n / ncol;
	double* old_value = NULL;
	int* p_skip = NULL;

	if( !isNull(skip) ){
		old_value = (double*) R_alloc(n_skip*ncol, sizeof(double));
		p_skip = INTEGER_POINTER(skip);
		for(int k=ncol-1; k>=0; --k){
			for(int i=n_skip-1; i>=0; --i){
				//Rprintf("skip %i x %i\n", i, k);
				int is = p_skip[i]-1;
				double val = p_x[k*nrow + is];
				old_value[k*n_skip + i] = val;
			}
		}
	}

	// apply limit inf to all values
	double* p_x2 = p_x + n-1;
	for(int i=n-1; i>=0; --i){
		if( *p_x2 < lim )
			*p_x2 = lim;
		--p_x2;
	}
	p_x2 = NULL;

	// restore skipped values
	if( !isNull(skip) ){
		for(int k=ncol-1; k>=0; --k){
			for(int i=n_skip-1; i>=0; --i){
				//Rprintf("restore %i x %i\n", i, k);
				int is = p_skip[i]-1;
				p_x[k*nrow + is] = old_value[k*n_skip + i];
			}
		}
	}


	// return modified x
	return x;
}

/** Apply inequality constraints in place. */
SEXP ptr_neq_constraints(SEXP x, SEXP constraints, SEXP ratio, SEXP value){

	double* p_x = ( isNull(x) ? NULL : NUMERIC_POINTER(x) );
	double d_ratio = isNull(ratio) ? 0 : *NUMERIC_POINTER(ratio);
	double* p_value = ( isNull(value) ? NULL : NUMERIC_POINTER(value) );
	double eps = sqrt(DOUBLE_EPS);

	// get dimensions
	int ncol = isNull(GET_DIM(x)) ? 1 : INTEGER(GET_DIM(x))[1];
	int nrow = isNull(GET_DIM(x)) ? length(x) : INTEGER(GET_DIM(x))[0];
	int nc = length(constraints);
	if( nc != ncol )
		error("There must be as many elements in list `constraints` as columns in `x`.");

	// apply each set of constraints (from first to last)
	double* _xj = p_x; // pointer to marked column
	double* _x_last = p_x + (ncol - 1) * nrow; // pointer to last column
	for(int j=0; j<nc; ++j){
		SEXP c_j = VECTOR_ELT(constraints, j);
		int n = length(c_j);
		int* p_i = INTEGER_POINTER(c_j);

		// apply the constraint on each row in the set
		for(int k=n-1; k>=0; --k){
			double lim = d_ratio != 0.0 ? _xj[p_i[k]-1] / d_ratio - eps : 0.0;
			if( lim < 0 )
				lim = 0;

			// apply constraints on each column
			// pointer to current row in last column
			double* _xi = _x_last + p_i[k]-1;
			for(int l=ncol-1; l>=0; --l){
				//Rprintf("Before: xi=%f > lim=%f ? => ", lim, *_xi);
				if( l != j && *_xi > lim ){ // constrain column to 'lim'
					*_xi = lim;
				}else if( l == j && p_value != NULL ){ // constrain column to 'value'
					*_xi = *p_value;
				}
				//Rprintf("xi=%f\n", *_xi);
				// move to previous column
				_xi -= nrow;
			}
			_xi = NULL;
		}
		// move to next marked column
		_xj += nrow;
	}

	// return modified x
	return x;
}


template<class T> inline void colMin(T* x, int n, int p, T* res, const T& NA_value){

	// do nothing if there is no data or fill with NAs
	if( n <= 0 ){
		if( p <= 0 )
			return;
		for(int j=p-1; j>=0; --j, ++res)
			*res = NA_value;
	}

	for(int j=p-1; j>=0; --j, ++res){
		*res = *(x++);
		for(int i=n-2; i>=0; --i, ++x){
			if( *res > *x )
				*res = *x;
		}
	}

}

template<class T> inline void colMax(T* x, int n, int p, T* res, const T& NA_value){

	// do nothing if there is no data or fill with NAs
	if( n <= 0 ){
		if( p <= 0 )
			return;
		for(int j=p-1; j>=0; --j, ++res)
			*res = NA_value;
	}

	for(int j=p-1; j>=0; --j, ++res){
		*res = *(x++);
		for(int i=n-2; i>=0; --i, ++x){
			if( *res < *x )
				*res = *x;
		}
	}

}

/**
 * Minimum per column
 */
SEXP colMin(SEXP x){

	SEXP ans, dims;

	// check that the argument is a matrix
	dims = GET_DIM(x);
	if (dims == R_NilValue)
		error("a matrix-like object is required as argument to 'colMin'");
	// check that it is a numeric data
	if (!isNumeric(x))
		error("a numeric object is required as argument to 'colMin'");

	// get the dimension of the input matrix
	int n = INTEGER(dims)[0];
	int p = INTEGER(dims)[1];

	if( TYPEOF(x) == REALSXP ){
		// allocate memory for the result (a vector of length the number of columns of x)
		PROTECT(ans = allocVector(REALSXP, p));
		colMin(NUMERIC_POINTER(x), n, p, NUMERIC_POINTER(ans), NA_REAL);
		UNPROTECT(1);
	}
	else{
		// allocate memory for the result (a vector of length the number of columns of x)
		PROTECT(ans = allocVector(INTSXP, p));
		colMin(INTEGER_POINTER(x), n, p, INTEGER_POINTER(ans), NA_INTEGER);
		UNPROTECT(1);
	}

	return ans;
}

/**
 * Maximum per column
 */
SEXP colMax(SEXP x){

	SEXP ans, dims;

	// check that the argument is a matrix
	dims = GET_DIM(x);
	if (dims == R_NilValue)
		error("a matrix-like object is required as argument to 'colMax'");
	// check that it is a numeric data
	if (!isNumeric(x))
		error("a numeric object is required as argument to 'colMax'");

	// get the dimension of the input matrix
	int n = INTEGER(dims)[0];
	int p = INTEGER(dims)[1];

	if( TYPEOF(x) == REALSXP ){
		// allocate memory for the result (a vector of length the number of columns of x)
		PROTECT(ans = allocVector(REALSXP, p));
		colMax(NUMERIC_POINTER(x), n, p, NUMERIC_POINTER(ans), NA_REAL);
		UNPROTECT(1);
	}
	else{
		// allocate memory for the result (a vector of length the number of columns of x)
		PROTECT(ans = allocVector(INTSXP, p));
		colMax(INTEGER_POINTER(x), n, p, INTEGER_POINTER(ans), NA_INTEGER);
		UNPROTECT(1);
	}

	return ans;
}

#endif //END ifdef NMF_UTILS_H
