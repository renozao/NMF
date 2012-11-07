	SEXP res;
	int nprotect = 0;

	double eps_val = *NUMERIC_POINTER(eps);

	// retrieve dimensions
	int n = INTEGER(GET_DIM(v))[0];
	int r = INTEGER(GET_DIM(w))[1];
	int p = INTEGER(GET_DIM(v))[1];

	// duplicate H (keeping attributes)
	PROTECT(res = duplicate(h)); nprotect++;
	DUPLICATE_ATTRIB(res, h);

	// define internal pointers
	double* pV = NUMERIC_POINTER(v);
	double* pW = NUMERIC_POINTER(w);
	double* pH = NUMERIC_POINTER(h);
	double* p_res = NUMERIC_POINTER(res);

	// extra variables in the case of an optional offset
	#ifdef WITH_OFFSET
	double *pOffset = NULL, *den_addon = NULL;

	if( offset != R_NilValue ){
		pOffset = NUMERIC_POINTER(offset);
		den_addon = (double*) R_alloc(r, sizeof(double));
		memset(den_addon, 0, r);
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
		for( int j=r-1; j>=0; --j){
			temp = 0;
			for( int u=n-1; u>=0; --u){
				temp += p_row[u] * pW[u + j*n];

				#ifdef WITH_OFFSET
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
			#ifdef WITH_OFFSET
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
