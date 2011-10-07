#ifndef RNGTOOLS_H // include header only once
#define RNGTOOLS_H

// check GCC version for compatibility with Rcpp
#if (__GNUC__ < 4 ) || ((__GNUC__ == 4) && (__GNUC_MINOR__ < 2) )
#define NO_RNGTOOLS
#endif

#ifdef NO_RNGTOOLS
	#warning "SKIP COMPILATION of `rngtools`"
#else

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>
#include <R_ext/Rdynload.h>
#include <RcppCommon.h>
#include <Rcpp.h>

extern "C" {

	/** Set/Get the current provider of user-supplied RNGs */
	SEXP rngtools_setProvider(SEXP name, SEXP reload_hooks=R_NilValue);
	SEXP rngtools_getProvider();

	/** Set/Get the next current provider of user-supplied RNGs */
	SEXP rngtools_setNextProvider(SEXP name);
	SEXP rngtools_getNextProvider();

	/** Tell if the RNGwrapper library is in use */
	SEXP rngtools_detectProvider();

}

//SEXP rngtools_getSymbol(SEXP pkg){
//
//	DllInfo* info = (DllInfo *) R_ExternalPtrAddr(pkg);
//
//	int i;
//	int n = info->numCallSymbols;
//	for(i = 0; i < n;) {
//		for(int k=0; k<6 && i<n; k++, i++)
//			Rprintf(" %s", info->CallSymbols[i].name);
//		Rprintf("\n");
//	}
//	return( R_NilValue );
//}

#define NDEBUG
#ifndef NDEBUG
#define DEBUG_LEVEL 2
#define DEBUG_VERB(x) x
#if DEBUG_LEVEL > 1
#define DEBUG_VERB2(x) x
#else
#define DEBUG_VERB2(x)
#endif
#else
#define DEBUG_VERB(x)
#define DEBUG_VERB2(x)
#endif


#define RNG_unif_rand "user_unif_rand"
#define RNG_unif_init "user_unif_init"
#define RNG_unif_nseed "user_unif_nseed"
#define RNG_unif_seedloc "user_unif_seedloc"

#define RNGWRAP_ENV Environment::namespace_env("NMF")
//#define RNGWRAP_ENV Environment("package:rngtools")
//#define RNGWRAP_ENV Environment::global_env()


/**
 * Helper class to locally backup and restore R-core RNG settings
 */
class RNGRestorationPoint{

	Rcpp::IntegerVector m_seed;

public:
	RNGRestorationPoint(){
		using namespace Rcpp;

		// Update Random seed (on R side)
		PutRNGstate();

		// backup value of .Random.seed
		DEBUG_VERB( Rprintf("### Backup .Random.seed\n"); )
		Environment glob = Environment::global_env();
		IntegerVector Random_seed = glob[".Random.seed"];
		m_seed = clone(Random_seed);
	}

	~RNGRestorationPoint(){
		using namespace Rcpp;
		DEBUG_VERB( Rprintf("### Restore .Random.seed\n"); )

		// Restore value of .Random.seed
		Environment glob = Environment::global_env();
		glob.assign(".Random.seed", m_seed);

		// Update internal state (on C side)
		GetRNGstate();
	}

};

/**
 * Load user-supplied hook
 */
template<typename T> void load_RNG_hook(T &fun, const char* provider, const char* hook, bool required = false){

	// check input
	if( *provider == 0  ){
		error("load_user_RNG - Invalid empty user-supplied RNG provider");
	}
	if( *hook == 0  ){
		error("load_user_RNG - Invalid empty user-supplied RNG hook");
	}

	DEBUG_VERB( Rprintf("rngtools::load_user_RNG - Lookup hook '%s' in provider '%s' ... ", hook, provider); )
	fun = (T) R_FindSymbol(hook, provider, NULL);
	if( required && !fun ){
		Rprintf("ERROR\n");
		std::string err("Could not find hook '");
		err += hook;
		err += "' in provider '";
		err += provider;
		err += "'";
		error(err.c_str());
	}
	DEBUG_VERB( else if( !required ){
		Rprintf("SKIP [not found]\n");
	}
	else
		Rprintf("OK [%p]\n", fun);
	)
}

/**
 * Static object that contains the data (i.e. pointers to RNG hooks)
 * from the of the package that provides the current definition for
 * user-supplied RNG hooks: user_unif_rand, etc...
 */
typedef void (*UnifInitFun)(Int32);
class RNGprovider{

public:
	std::string name;
	DL_FUNC User_unif_fun;
	DL_FUNC User_unif_nseed;
	DL_FUNC	User_unif_seedloc;
	UnifInitFun User_unif_init;

	/**
	 * Default constructor.
	 */
	RNGprovider()
	: name(""), User_unif_fun(NULL)
	, User_unif_nseed(NULL), User_unif_seedloc(NULL)
	, User_unif_init(NULL){
	}

	void reset(){
		name = "";
		User_unif_fun = User_unif_nseed = User_unif_seedloc = NULL;
		User_unif_init = NULL;
	}

	SEXP set(SEXP provider, const char* desc = "TODO"){

		std::string new_provider(Rcpp::as<const char*>(provider));
		DEBUG_VERB( Rprintf("rngtools::setProvider - Try changing RNG %s provider from '%s' to '%s' ... "
							, desc, name.c_str(), new_provider.c_str()); )

		SEXP old = Rcpp::wrap(name);
		if( name == new_provider ){
			DEBUG_VERB( Rprintf("SKIP\n"); )
			return( old );
		}
		DEBUG_VERB( else Rprintf("\n"); )

		// load data from the new provider
		RNGprovider tmp(new_provider.c_str());
		// simple copy the new provider data into the static variable
		*this = tmp;

		DEBUG_VERB(	Rprintf("rngtools::setProvider - DONE [provider: '%s']\n", name.c_str()); )

		return( old );

	}

	/**
	 * Load data from a package that provides RNG hooks
	 */
	void load(const char* provider){

		if( !provider || *provider == 0 )
			reset();
		else{
			name = provider;
			// load each hook from the package
			load_RNG_hook(User_unif_fun, provider, RNG_unif_rand, true);
			load_RNG_hook(User_unif_init, provider, RNG_unif_init);
			load_RNG_hook(User_unif_nseed, provider, RNG_unif_nseed);
			load_RNG_hook(User_unif_seedloc, provider, RNG_unif_seedloc);
		}

	}

	/**
	 * Constructor that loads data from a package that provides RNG hooks
	 */
	RNGprovider(const char* provider){
		load(provider);
	}
};

/**
 * Static local variable that stores data on the currently "active"
 * user-supplied RNG.
 *
 * The currently "active" library is the library that would be called by the
 * wrapper hooks if the RNGkind were "user-supplied" AND the RNGwrapper library
 * is the last loaded RNG library.
 *
 */
static RNGprovider _current_provider;

/*
 * Static local variable used to temporary store data on the next current RNG
 * library.
 *
 * It is used in .setRNG to ensure that the following call to RNGkind uses the
 * correct RNG libraries:
 * - use the current provider when calling runif(1) to generate the seed
 * - use the next current when possibly calling user_unif_init, etc...
 */
static RNGprovider _next_provider;

/** Fixes issues that appear when the RNGwrapper library is reloaded.
 *
 * This function is only called when reloading the RNGwrapper library AND when
 * some user-supplied RNG is in use.
 * Its purpose is to update the cached pointers to RNG hooks
 *
 * This is performed as follows:
 * - set a static flag to TRUE that makes user_unif_rand reset it to FALSE.
 * - call unif_rand
 * - if the RNGwrapper is in use, then user_unif_rand will be called and the
 * static flag will change value. In any case on stores the
 */
static bool _do_runif_calltest = false;

#define RUNIF_CACHE_MAX 3
class cache{

	double _cache[RUNIF_CACHE_MAX];
	int _ncached;
	int _ncalled;

public:

	/** Default constructor */
	cache()
	: _ncached(0), _ncalled(0){
		reset();
	}

	void reset(){
		// reset cache as it will not be used
		_ncached = _ncalled = 0;
		for(int k=0; k<RUNIF_CACHE_MAX; k++)
			_cache[k] = -1.0;
	}

	inline int size() const{
		return( _ncached );
	}

	inline int ncall() const{
		return( _ncalled );
	}

	inline double add(double value){
		_cache[_ncached++] = value;
		return value;
	}

	inline bool is_empty() const{
		return( _ncached == 0 );
	}

	inline double* get(){
		if( is_empty() )
			return( NULL );

		int idx = _ncalled++;
		if( _ncalled == _ncached )
			_ncalled = _ncached = 0;

		return( &_cache[idx] );
	}

	double operator[](int i)const{
		return( _cache[i] );
	}

	double& operator[](int i){
		return( _cache[i] );
	}
};

static cache _runif_cache;
//static double _runif_cache[RUNIF_CACHE_MAX] = {-1,-1,-1};
//static int _runif_ncache = 0; // number of results from runif that were cached
//static int _cache_ncall = 0; // number of values from runif that were retrieved from the cache

//static void reset_runif_cache(){
//	// reset cache as it will not be used
//	_runif_ncache = _cache_ncall = 0;
//	for(int k=0; k<RUNIF_CACHE_MAX; k++)
//		_runif_cache[k] = -1;
//}

// static flag to allow running a fake set.seed
static bool _runif_fake_init = false;
static bool _runif_init_cache = false;

/**
 * Identify the currently active user-supplied RNG provider by the length and
 * value of .Random.seed.
 *
 * This lookup is not guaranteed to succeed as RNG libraries often do not provide
 * the hook 'user_unif_seedloc', but if it does identify a unique provider then
 * it makes the RNGwrapper library compatible with the given RNG library.
 *
 */
static int searchRNG_bySeed(const Rcpp::CharacterVector& libs, Rcpp::LogicalVector& flags){

	using namespace Rcpp;

	DEBUG_VERB( Rprintf("# Search matching seed ... \n"); )

	int unique_match = -1;
	bool all_have_seedloc = true;
	for(int i=libs.length()-1; i>=0; --i){
		const char* provider = libs[i];

		// discard indexes already flagged by previous searches
		if( !flags[i] ){
			DEBUG_VERB( Rprintf("\t* Provider '%s' ... SKIP [discarded]\n", provider); )
			continue;
		}
		DEBUG_VERB( Rprintf("\t* Provider '%s' ... ", provider); )

		// load the user_unif_seedloc hook from provider
		DL_FUNC seedloc = R_FindSymbol(RNG_unif_seedloc, provider, NULL);
		if( !seedloc ){
			all_have_seedloc = false;
			DEBUG_VERB( Rprintf("SKIP [no seedloc hook]\n"); )
			continue;
		}
		// load the user_unif_nseed hook from provider
		DL_FUNC nseed = R_FindSymbol(RNG_unif_nseed, provider, NULL);
		if( !nseed ){
			all_have_seedloc = false;
			DEBUG_VERB( Rprintf("SKIP [no nseed hook]\n"); )
			continue;
		}

		// load value of .Random.seed
		IntegerVector Random_seed( findVarInFrame(R_GlobalEnv, R_SeedsSymbol) );
		int ns = *(int*) nseed();

		// look for differences in the seed's length
		if( ns != Random_seed.length()-1 ){
			DEBUG_VERB( Rprintf("FAILED [length]\n"); )
			// flag the index to be discarded in the following searches
			flags[i] = false;
			continue;
		}

		// look for differences in the seed's value
		Int32* seed = static_cast<Int32*>(seedloc());
		bool diff = false;
		for(int k=ns-1; k>=0; --k){
			if( seed[k] != Random_seed[k+1] ){
				diff = true;
				break;
			}
		}
		if( diff ){
			DEBUG_VERB( Rprintf("FAILED [value]\n"); )
			// flag the index to be discarded in the following searches
			flags[i] = false;
			continue;
		}

		DEBUG_VERB( Rprintf("MATCH\n"); )

		// flag the index that uniquely passed
		if( unique_match < 0 ) unique_match = i; // it's the first index that passed: store it
		else unique_match = libs.length(); // it's the second index that passed: stop storing

	}
	DEBUG_VERB( Rprintf("# DONE\n"); )

	// return something meaningful if a unique match was found and all the libraries
	// provides a hook for user_unif_seedloc.
	if( all_have_seedloc && unique_match >= 0 && unique_match < libs.length() )
		return unique_match;
	else
		return -1;
}

/**
 * Load the external or custom state vectors for a vector of RNG providers.
 *
 */
void load_RNGstates(const Rcpp::CharacterVector& libs, const Rcpp::LogicalVector& flags, std::vector<SEXP>& states){

	using namespace Rcpp;

	BEGIN_RCPP

		// reset the state vector
		states.clear();

		// load some environments
		Environment base = Environment("package:base");
		Function r_Call = base[".Call"];
		Environment ns_digest = Environment("package:digest");
		Function r_digest = ns_digest["digest"];
		Environment RNGwrap_env = RNGWRAP_ENV;
		Environment glob = Environment::global_env();

		for( int i=0; i<libs.length(); ++i){
			const char* provider = libs[i];

			states.push_back(R_NilValue);
			if( !flags[i] ){
				DEBUG_VERB( Rprintf("\t* Provider '%s' ... SKIP [discarded]\n", provider); )
				continue;
			}

			DEBUG_VERB( Rprintf("\t* Provider '%s' ... ", provider); )
			if( !strcmp(provider, "rlecuyer") ){

				SEXP tmp = r_Call("r_get_current_stream", Named("PACKAGE", "rlecuyer"));
				NumericVector seed(clone(VECTOR_ELT(tmp,0)));
				states[i] = wrap(seed);
				DEBUG_VERB(	Rprintf("OK [state: %s]\n", as<const char*>(r_digest(seed)) ); )

			}
			else if( !strcmp(provider, "rstream") ){

				Environment rstream_ns = Environment::namespace_env("rstream");
				if( !rstream_ns.exists(".rstream.envir") ){
					DEBUG_VERB( Rprintf("SKIP [not initialised]\n"); )
					continue;
				}
				Environment rstream_env = rstream_ns[".rstream.envir"];
				if( !rstream_env.exists(".rstream.current") ){
					DEBUG_VERB( Rprintf("SKIP [not in use]\n"); )
					continue;
				}
				SEXP current = rstream_env[".rstream.current"];
					RObject stream(current);
					try{

						SEXP tmp = r_Call("R_RngStreams_GetData", stream.slot("stream"), Named("PACKAGE","rstream"));
						NumericVector seed(clone(tmp));
						states[i] = wrap(seed);
						DEBUG_VERB(	Rprintf("OK [state: %s]\n", as<const char*>(r_digest(seed)) ); )

					}catch(...){
						DEBUG_VERB( Rprintf("SKIP [bad pointer]\n"); )
					continue;
				}

			}else{

				IntegerVector Random_seed = glob[".Random.seed"];
				// only load the state if Random_seed contains specific data (first element gives the RNG kinds)
				if( Random_seed.length() > 1 ){
					states[i] = wrap(clone(Random_seed));
					DEBUG_VERB(
						DEBUG_VERB(	Rprintf("OK [.Random.seed: %s]\n", as<const char*>(r_digest(Random_seed)) ); )
					)
				}else{
					DEBUG_VERB( Rprintf("SKIP [unknown]\n"); )
				}

			}

		}

	VOID_END_RCPP

}


/**
 * Runs a fake set.seed to force re-loading of the wrapper hook for user_unif_rand
 */
static void fake_set_seed(const char* kind = NULL){

	using namespace Rcpp;
	DEBUG_VERB( Rprintf("# Run fake initialization to refresh RNG.c::User_unif_fun ... \n"); )
	Environment base = Environment("package:base");
	Function setseed = base["set.seed"];
	// turn on the fake run flag
	_runif_fake_init = true;
	if( kind != NULL ){
		// this will make 2 fake calls to the init hooks
		setseed(0, Named("kind", kind));
	}else{
		// this will make a single fake call to the init hooks
		setseed(0);
	}
	// turn off the fake run flag
	_runif_fake_init = false;
	DEBUG_VERB( Rprintf("DONE\n"); )

}

/**
 * Auto-detect the currently active user-supplied RNG provider.
 */
SEXP rngtools_detectProvider(){

	using namespace Rcpp;

	BEGIN_RCPP

	DEBUG_VERB( Rprintf("# RNGwrapper: Trying to detect current user-supplied RNG ...\n"); )

	// identify the RNG library that was in use before reloading
	//1. get all the RNG libs (they are in loading order: the last one is the
	// more likely to be the one that was in use
	Environment RNGwrap_env = RNGWRAP_ENV;
	Function RNGlibs = RNGwrap_env["RNGlibs"];
	DEBUG_VERB( Rprintf("# Loading RNG libraries\n"); )
	CharacterVector libs = as<CharacterVector>(RNGlibs());

	// if there is no loaded RNG provider then reset the hooks
	if( libs.length() == 0 ){
		DEBUG_VERB( Rprintf("# Detected no RNG provider\n"); )
		rngtools_setProvider(wrap(""), wrap(true));
		return( wrap("") );
	}
	// if there is only one possible library then this must be the one
	if( libs.length() == 1 ){
		const char* provider = libs[0];
		DEBUG_VERB( Rprintf("# Detected unique RNG provider '%s'\n", provider); )

		// load the detected RNG provider, reloading the hooks
		rngtools_setProvider( wrap(provider), wrap(true));

		// return the identified provider
		return( wrap(provider) );
	}

	DEBUG_VERB( Rprintf("# Allocate memory for hooks\n"); )
	std::vector<DL_FUNC> lib_hooks;
	lib_hooks.reserve(libs.length());

	// create a vector of flags: 1=lookup, 0=discard
	LogicalVector flags(Dimension(libs.length()), true);
	int match = -1;

	// Load R-level function: runif
	Environment stats = Environment("package:stats");
	Function stats_runif = stats["runif"];

	//2. Search a matching RNG by seed
	match = searchRNG_bySeed(libs, flags);

	match = -1;
	// no unique match was found: carry on the search
	if( match < 0 ){

		// find index for currently wrapped RNG
		DEBUG_VERB( Rprintf("# Lookup for currently wrapped RNG '%s' ... ", _current_provider.name.c_str()); )
		int current_index = -1;
		for(int i=libs.length()-1; i>=0; --i){
			const char* provider = libs[i];
			if( _current_provider.name ==  provider ){
				current_index = i;
				break;
			}
		}
		if( current_index < 0 ){
			DEBUG_VERB( Rprintf("ERROR\n"); )
			error("Could not find the currently wrapped RNG provider in the provider list");
		}
		DEBUG_VERB( Rprintf("OK [index:%i]\n", current_index); )

		DEBUG_VERB( Rprintf("# Store external state for each RNG provider ... \n"); )
		// 2. Load states for known RNG libraries, call `runif` and see if it changes
		// the state.
		std::vector<SEXP> states_0;
		states_0.reserve(libs.length());
		load_RNGstates(libs, flags, states_0);
		DEBUG_VERB( Rprintf("DONE\n"); )

		// 3. check if the active RNG was wrapped already, in
		// this case, the provider is known as it was restored in the R function
		// RNGwrap.
		DEBUG_VERB( Rprintf("# Drawing from runif ... "); )
		_do_runif_calltest = true;
		_runif_cache.reset();
		_runif_cache.add( as<double>( stats_runif(1) ) );
		DEBUG_VERB( Rprintf("OK [%f]\n", _runif_cache[0]); )

		// if there is now a cached value it means that the RNG was already wrapped:
		// run a fake set.seed to enforce User_unif_fun to be refreshed in RNG.c [line 240]
		// NB: because the RNGwrapper library is now the last loaded RNG library
		// the other hooks are actually already active.
		DEBUG_VERB( Rprintf("# Check if active RNG is the one currently wrapped ... "); )
		if( !_do_runif_calltest ){
			DEBUG_VERB( Rprintf("YES\n"); )
			match = current_index;
		}else{
			_do_runif_calltest = false;
			DEBUG_VERB( Rprintf("NO\n"); )
		}

		// 4. Reload the states for the RNG libraries and look for any change
		if( match < 0 ){
			DEBUG_VERB( Rprintf("# Load internal states for each RNG provider ... \n"); )
			std::vector<SEXP> states_1;
			states_1.reserve(libs.length());
			load_RNGstates(libs, flags, states_1);
			DEBUG_VERB( Rprintf("DONE\n"); )

			DEBUG_VERB( Rprintf("# Check changes in internal states ... \n"); )
			int unique_match = -1;
			for( int i=libs.length()-1; i>=0; --i){
				const char* provider = libs[i];

				DEBUG_VERB( Rprintf("\t* Provider '%s' ... ", provider); )
				//Environment base = Environment("package:base");
				//Function _print = base["print"];
				//_print(states_0[i]);
				//_print(states_1[i]);

				if( !flags[i] ){
					DEBUG_VERB( Rprintf("SKIP [discarded]\n"); )
					continue;
				}

				// look for differences only if the states have been loaded
				if( !isNull(states_0[i]) && !isNull(states_1[i]) ){
					bool diff = false;
					if( TYPEOF(states_0[i]) == REALSXP ){
						NumericVector s0(states_0[i]);
						NumericVector s1(states_1[i]);
						diff = any( s1 != s0 ).is_true();
					}else{
						IntegerVector s0(states_0[i]);
						IntegerVector s1(states_1[i]);
						diff = any( s1 != s0 ).is_true();
					}
					if( diff ){// found differences
						DEBUG_VERB( Rprintf("YES\n"); )
						if( unique_match < 0 ) unique_match = i;
						else unique_match = libs.length();

					}else{
						flags[i] = false;
						DEBUG_VERB( Rprintf("NO [identical]\n"); )
					}
				}
				DEBUG_VERB( else Rprintf("SKIP [not comparable]\n"); )

				// set match to a meaningful value if a unique match was found
				match = unique_match >= 0 && unique_match < libs.length() ? unique_match : -1;
				if( unique_match >= libs.length() )
					warning("Multiple RNG providers changed their state after runif draw");
			}
		}
	}

	// 5. Draw once from each RNG provider and once from runif:
	// look for possibility of prediction of the result by pointer obtain
	// from the first draw
	if( match < 0 ){

		// allocate memory for storing the drawn values and pointers
		NumericVector runif_val(libs.length());
		std::vector<double*> runif_ptr;
		runif_ptr.reserve(libs.length());

		DEBUG_VERB( Rprintf("# Drawing once from each RNG provider ... \n"); )
		for( int i=libs.length()-1; i>=0; --i){

			const char* provider = libs[i];
			runif_val[i] = -1;
			runif_ptr[i] = NULL;

			if( !flags[i] ){
				DEBUG_VERB( Rprintf("\t* Drawing from provider '%s' ... SKIP [discarded]\n", provider); )
				continue;
			}

			DEBUG_VERB( Rprintf("\t* Drawing from provider '%s' ... ", provider); )
			// load the user_unif_rand hook from RNG library i
			lib_hooks[i] = R_FindSymbol(RNG_unif_rand, provider, NULL);

			// try to get a value from it if present
			if( !lib_hooks[i] ){
				DEBUG_VERB( Rprintf("FAILED [invalid hook]\n"); )
				continue;
			}

			try{
				// call the provider's unif_rand hook
				GetRNGstate();
				double* val = runif_ptr[i] = static_cast<double*>( lib_hooks[i]() );
				PutRNGstate();
				if( !val ){
					DEBUG_VERB( Rprintf("FAILED [null pointer]\n"); )
					continue;
				}
				runif_val[i] = *val;
				DEBUG_VERB( Rprintf("OK [%f]\n", *val); )
			}catch(...){
				// do nothing
			}
		}

		// put a dummy value for the second cached RNG value. It will be filled
		// with its correct value once the RNG is identified
		_runif_cache.add(-1);

		// Draw once again from the active RNG
		//NB: there is no issue with the cache being called since we already
		// tested that the RNG wrapper hook is not called by runif.
		DEBUG_VERB( Rprintf("# Drawing from current runif ... "); )
		double runif_val2 = _runif_cache.add(as<double>( stats_runif(1) ));
		DEBUG_VERB( Rprintf("OK [%f]\n", runif_val2); )

		DEBUG_VERB( Rprintf("# Check each RNG provider's ability to predict next drawn ... \n"); )
		int unique_match = -1;
		for( int i=libs.length()-1; i>=0; --i){

			const char* provider = libs[i];

			if( !flags[i] ){
				DEBUG_VERB( Rprintf("\t* Comparing values for provider '%s' ... SKIP [discarded]\n", provider); )
				continue;
			}

			DEBUG_VERB( Rprintf("\t* Comparing values for provider '%s' ... ", provider); )

			// skip bad values
			if( runif_ptr[i] == NULL || runif_val[i] == -1 ){
				DEBUG_VERB( Rprintf("FAILED [invalid values]\n"); )
				continue;
			}

			try{
				// if the previous pointer predicts the value, we've got a match
				if( *runif_ptr[i] == runif_val2 ){

					DEBUG_VERB( Rprintf("YES [value]\n"); )
					if( unique_match < 0 ) unique_match = i;
					else if(unique_match < libs.length() )
						unique_match = libs.length() + unique_match;

				}else if( runif_ptr[i]+1 != NULL && *(runif_ptr[i]+1) == runif_val2 ){// maybe the internal state only shifts

					DEBUG_VERB( Rprintf("YES [shifted value]\n"); )
					if( unique_match < 0 ) unique_match = i;
					else if(unique_match < libs.length() )
						unique_match = libs.length() + unique_match;

				}
			}catch(...){
				// do nothing
			}
		}

		// set match to a meaningful value: use first match in case of multiple match
		if( unique_match >= libs.length() ){
			warning("Multiple RNG providers changed their pointer's value after runif draw: using first match");
			match = unique_match - libs.length();
		}
		else if( unique_match >= 0 )
			match = unique_match;

		// set the second cached RNG value with its correct value
		if( match >= 0 )
			_runif_cache[1] = runif_val[match];
	}


	const char* provider = NULL;

	// if the active RNG provider was not found, then set the last loaded RNG
	// as the currently wrapped provider and throw a warning about it.
	//NB: this means that the active RNG is still called by the R function `runif`.
	// The change of provider will only occur at the next call to `RNGkind` or
	// `set.seed`.
	if( match < 0 ){
		const char* last = as<const char*>(RNGlibs(1));
		warning("Could not identify the active user-supplied RNG provider: "
				"using the last loaded RNG library '%s' as the next provider."
				" [NB: restoration will not be possible]"
				, last);

		// reset the next RNG provider
		rngtools_setNextProvider( wrap(last) );
		// return an empty string
		return( R_NilValue );

	}else{
		// the current user-supplied RNG provider was identified
		provider = libs[match];
		DEBUG_VERB( Rprintf("# Detected RNG provider '%s'\n", provider); )

		// load the detected RNG provider, reloading the hooks
		rngtools_setProvider( wrap(provider), wrap(true));

		// return the identified provider
		return( wrap(provider) );
	}


	END_RCPP

}

/**
 * Wrapper call to the current user-supplied hook 'unif_user_rand'
 */
double* user_unif_rand(){

	DEBUG_VERB2( Rprintf("Calling hook '%s' from provider '%s'\n", RNG_unif_rand, _current_provider.name.c_str()); )

	// check that the hook was correctly initialized
	if( !_current_provider.User_unif_fun ){
		std::string err("Current hook 'user_unif_rand' is invalid [provider:'");
		err += _current_provider.name;
		err += "']";
		error(err.c_str());
		return( NULL );
	}else if( _do_runif_calltest ){ // switch test flag to show that this hook was call
		DEBUG_VERB( Rprintf("# This is a testing identification call\n"); )
		_do_runif_calltest = false;
	}else if( !_runif_cache.is_empty() ){

		int n = _runif_cache.ncall()+1;
		int ntotal = _runif_cache.size();
		double* res = _runif_cache.get();
		Rprintf("# Using CACHED runif value %i/%i: %f\n", n, ntotal, *res);
		return res;
	}

	// call the RNG provider's hook
	return (double*) _current_provider.User_unif_fun();
}

/**
 * Wrapper call to the current user-supplied hook 'user_unif_init'
 */
void user_unif_init(Int32 seed){

	if( _runif_fake_init ){
		DEBUG_VERB( Rprintf("Calling hook %s(%i) from provider '%s' ... SKIP [fake call]\n"
				, RNG_unif_init, seed, _current_provider.name.c_str()); )
		return;
	}

	// reset the cache for runif
	_runif_cache.reset();

	// use and reset next provider if set
	if( _next_provider.name != "" ){
		_current_provider = _next_provider;
		_next_provider.reset();
	}

	DEBUG_VERB( Rprintf("Calling hook %s(%i) from provider '%s' ... "
			, RNG_unif_init, seed, _current_provider.name.c_str()); )

	// Use a "bug" in R-core RNG management to enforce the use of the RNGwrapper
	// library:
	// in RNG.c, 'user_unif_init' is looked up in the loaded DLL table when calling
	// either set.seed or RNGkind. Some RNG libraries do not provide this hook
	// so that the RNGwrapper hook might be called in their place.
	// This allows us to test if this call to 'user_unif_init' is consistent
	// with the active user-supplied RNG provider which necessarily is the last
	// loaded library that provides 'user_unif_rand':
	// if not, one forces the RNGwrapper library to reload the correct RNG provider.
	// The following calls to 'user_unif_nseed' and 'user_unif_seedloc' should
	// correctly resolve to the RNGwrapper hooks.
	//
	DL_FUNC urand = R_FindSymbol(RNG_unif_rand, "", NULL);
	// this should return a non null function pointer as either RNGwrapper or
	// or the library responsible for this call to be made (a test for non-nullity
	// is already made in RNG.c [line 240 in R-2.13.0] before calling user_unif_init)
	if( !urand )
		error("rngtools::unif_init - Could not find symbol 'user_unif_rand'");
	if( urand != (DL_FUNC) user_unif_rand ){
		using namespace Rcpp;

		BEGIN_RCPP

			DEBUG_VERB( Rprintf("\n") );

			//Use RNGlib and not getNativeSymbolInfo due to bug in Windows which returns $ package=NULL
			Environment env = RNGWRAP_ENV;
			Function RNGlib = env["RNGlib"];
			const char* provider = as<const char*>(RNGlib());

			Rprintf("# rngtools::unif_init - Detected unwrapped user-supplied RNG provider '%s'\n", provider);

			// Wrap the detected provider			
			Function RNGwrap = env["RNGwrap"];
			RNGwrap(provider);
			// NOTES:
			// 1. We do not pass the value of seed, since we know this RNG provider
			// has no hook user_unif_init.
			// 2. the hook for user_unif_rand will still be provided by the
			// original library, until the next call to set.seed. However the
			// initialisation procedure in RNG.c [line 239] will carry on
			// and call the relevant wrapped hooks.
			return;

		VOID_END_RCPP
	}

	// call the provider's hook if present
	if( _current_provider.User_unif_init ){
		_current_provider.User_unif_init(seed);
		DEBUG_VERB( Rprintf("DONE\n") );
	}
	DEBUG_VERB( else Rprintf("SKIP [no hook]\n") );
}

/**
 * Wrapper call to the current user-supplied hook 'user_unif_nseed'
 */
static int nseed;
int* user_unif_nseed(void){


	DEBUG_VERB( Rprintf("Calling hook '%s' from provider '%s' ... ", RNG_unif_nseed, _current_provider.name.c_str()); )

	if( _current_provider.User_unif_nseed ){ // call the provider's hook if present
		DEBUG_VERB( Rprintf("NOW") );
		int* res = static_cast<int*>(_current_provider.User_unif_nseed());
		DEBUG_VERB( Rprintf(" [%i]\n", *res) );
		return res;
	}
	DEBUG_VERB( Rprintf("SKIP [no hook]\n") );

	// one needs to return a valid value in any case: 0 should not have any effect.
	nseed = 0;
	return( &nseed );
}

/**
 * Wrapper call to the current user-supplied hook 'user_unif_seedloc'
 */
static int seedata;
int* user_unif_seedloc(void){

	DEBUG_VERB( Rprintf("Calling hook '%s' from provider '%s' ... ", RNG_unif_seedloc, _current_provider.name.c_str()); )

	// call the provider's hook if present
	if( _current_provider.User_unif_seedloc ){
		DEBUG_VERB( Rprintf("NOW\n") );
		return static_cast<int*>(_current_provider.User_unif_seedloc());
	}

	// one needs to return a valid value that will not have any effect.
	DEBUG_VERB( Rprintf("SKIP [no hook]\n") );
	seedata = 0;
	return( &seedata );
}

/**
 * Set the current provider of user-supplied RNGs
 */
SEXP rngtools_setProvider(SEXP name, SEXP reload_hooks){

	BEGIN_RCPP

	SEXP ans = _current_provider.set(name, "current");

	// run a fake initialisation to reload the hooks if requested
	using namespace Rcpp;
	if( !isNull(reload_hooks) && as<bool>(reload_hooks) ){
		Environment base = Environment("package:base");
		Function rngKind = base["RNGkind"];
		CharacterVector kind = rngKind();

		DEBUG_VERB(
			Environment pkg_env = RNGWRAP_ENV;
			Function rnginfo = pkg_env["RNGinfo"];
			Rprintf("### Start RNG is:\n"); rnginfo();
		)

		// if the active RNG is user-supplied, then simply run a fake set.seed:
		// this will reload the RNG hooks and hence change the active user RNG
		bool is_base = strcmp(kind[0], "user-supplied");
		if( !is_base )
			fake_set_seed();
		else{
			// if the active RNG is a base RNG then on must ensure that the
			// user-supplied hooks are called, by temporarily changing the
			// active RNG to user-supplied
			DEBUG_VERB( Rprintf("### RNG is base\n"); )
			// setup restoration point for RNG settings
			RNGRestorationPoint rp;
			// do a fake set.seed changing to user-supplied
			fake_set_seed("user-supplied");
		}
		DEBUG_VERB( Rprintf("### end RNG is:\n"); rnginfo(); )
	}

	return ans;

	END_RCPP

} /* end of rngtools_setProvider() */

/**
 * Get the current provider of user-supplied RNGs
 */
SEXP rngtools_getProvider(){

	DEBUG_VERB( Rprintf("rngtools::getProvider - Return provider '%s'\n", _current_provider.name.c_str()); )
	return( Rcpp::wrap(_current_provider.name) );

} /* end of rngtools_getProvider() */


/**
 * Set the current provider of user-supplied RNGs
 */
SEXP rngtools_setNextProvider(SEXP name){

	return( _next_provider.set(name, "next") );

} /* end of rngtools_setNextProvider() */

/**
 * Get the current provider of user-supplied RNGs
 */
SEXP rngtools_getNextProvider(){

	DEBUG_VERB( Rprintf("rngtools::getProvider - Return provider '%s'\n", _next_provider.name.c_str()); )
	return( Rcpp::wrap(_next_provider.name) );

} /* end of rngtools_getProvider() */

#endif //END_:NO_RNGTOOLS
#endif //END_:RNGTOOLS_H
