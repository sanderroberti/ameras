#include <stdlib.h> // for NULL
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#ifdef __cplusplus
extern "C" {
#endif

extern SEXP _ameras_dldd_clogit(SEXP, SEXP);
extern SEXP _ameras_dldd_prophaz(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ameras_compute_ERCsum_clogit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ameras_compute_ERCsum_prophaz(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ameras_loglik_prophaz_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#ifdef __cplusplus
}
#endif

static const R_CallMethodDef CallEntries[] = {
  {"_ameras_dldd_clogit", (DL_FUNC) &_ameras_dldd_clogit, 2},
  {"_ameras_dldd_prophaz", (DL_FUNC) &_ameras_dldd_prophaz, 4},
  {"_ameras_compute_ERCsum_clogit", (DL_FUNC) &_ameras_compute_ERCsum_clogit, 6},
  {"_ameras_compute_ERCsum_prophaz", (DL_FUNC) &_ameras_compute_ERCsum_prophaz, 9},
  {"_ameras_loglik_prophaz_rcpp", (DL_FUNC) &_ameras_loglik_prophaz_rcpp, 6},
  {NULL, NULL, 0}
};

void R_init_ameras(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

