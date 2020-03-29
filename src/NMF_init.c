#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
/* .Call calls */
extern SEXP big_silhouette(SEXP, SEXP, SEXP, SEXP);
extern SEXP clone_object(SEXP);
extern SEXP c_colMax(SEXP);
extern SEXP c_colMin(SEXP);
extern SEXP divergence_update_H(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP divergence_update_W(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Euclidean_rss(SEXP, SEXP);
extern SEXP euclidean_update_H(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP euclidean_update_W(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP KL_divergence(SEXP, SEXP);
extern SEXP offset_euclidean_update_H(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP offset_euclidean_update_W(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ptr_address(SEXP);
extern SEXP c_ptr_isnil(SEXP);
extern SEXP ptr_neq_constraints(SEXP, SEXP, SEXP, SEXP);
extern SEXP ptr_pmax(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"big_silhouette",            (DL_FUNC) &big_silhouette,            4},
  {"clone_object",              (DL_FUNC) &clone_object,              1},
  {"c_colMax",                    (DL_FUNC) &c_colMax,                    1},
  {"c_colMin",                    (DL_FUNC) &c_colMin,                    1},
  {"divergence_update_H",       (DL_FUNC) &divergence_update_H,       6},
  {"divergence_update_W",       (DL_FUNC) &divergence_update_W,       6},
  {"Euclidean_rss",             (DL_FUNC) &Euclidean_rss,             2},
  {"euclidean_update_H",        (DL_FUNC) &euclidean_update_H,        7},
  {"euclidean_update_W",        (DL_FUNC) &euclidean_update_W,        8},
  {"KL_divergence",             (DL_FUNC) &KL_divergence,             2},
  {"offset_euclidean_update_H", (DL_FUNC) &offset_euclidean_update_H, 6},
  {"offset_euclidean_update_W", (DL_FUNC) &offset_euclidean_update_W, 6},
  {"ptr_address",               (DL_FUNC) &ptr_address,               1},
  {"c_ptr_isnil",                 (DL_FUNC) &c_ptr_isnil,                 1},
  {"ptr_neq_constraints",       (DL_FUNC) &ptr_neq_constraints,       4},
  {"ptr_pmax",                  (DL_FUNC) &ptr_pmax,                  3},
  {NULL, NULL, 0}
};

void R_init_NMF(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
