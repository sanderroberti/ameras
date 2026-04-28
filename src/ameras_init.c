#include <stdlib.h> // for NULL
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void C_dldd_clogit(void *, void *, void *, void *, void *);
extern void C_dldd_prophaz(double*, double*, int*, double*, int*, double*);
extern void C_loglik_prophaz_rcpp(double*, double*, double*, double*, int*, int*, int*, double*, double*);


static const R_CMethodDef CEntries[] = {
  {"C_dldd_clogit",    (DL_FUNC) &C_dldd_clogit,   5},
  {"C_dldd_prophaz",   (DL_FUNC) &C_dldd_prophaz,  6},
  {"C_loglik_prophaz_rcpp",     (DL_FUNC) &C_loglik_prophaz_rcpp,    9},
  {NULL, NULL, 0}
};

#ifdef __cplusplus
extern "C" {
#endif

extern SEXP _ameras_compute_ERCsum_clogit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
extern "C" {
#endif

extern SEXP _ameras_compute_ERCsum_prophaz(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#ifdef __cplusplus
}
#endif

static const R_CallMethodDef CallEntries[] = {
  {"_ameras_compute_ERCsum_clogit", (DL_FUNC) &_ameras_compute_ERCsum_clogit, 6},
  {"_ameras_compute_ERCsum_prophaz", (DL_FUNC) &_ameras_compute_ERCsum_prophaz, 8},
  {NULL, NULL, 0}
};

void R_init_ameras(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}


