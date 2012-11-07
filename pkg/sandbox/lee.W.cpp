	SEXP res;
	int nprotect = 0;

	double eps_val = *NUMERIC_POINTER(eps);

	// retrieve dimensions
	int n = INTEGER(GET_DIM(v))[0];
	int r = INTEGER(GET_DIM(w))[1];
	int p = INTEGER(GET_DIM(v))[1];

	// duplicate H (keeping attributes)
	PROTECT(res = duplicate(w)); nprotect++;
	DUPLICATE_ATTRIB(res, w);

	// define internal pointers to data
	double* pV = NUMERIC_POINTER(v);
	double* pW = NUMERIC_POINTER(w);
	double* pH = NUMERIC_POINTER(h);
	double* p_res = NUMERIC_POINTER(res);

	// extra variables in the case of an optional offset
	#ifdef WITH_OFFSET
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
			#ifdef WITH_OFFSET
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
