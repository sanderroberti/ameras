#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void C_dldd_clogit(void *, void *, void *, void *, void *);
extern void C_dldd_prophaz(double*, double*, int*, double*, int*, double*);
extern void C_compute_ERCmatrix_clogit(void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_compute_ERCmatrix_prophaz(void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_loglik_prophaz_rcpp(double*, double*, double*, double*, int*, int*, int*, double*, double*);


static const R_CMethodDef CEntries[] = {
  {"C_dldd_clogit",    (DL_FUNC) &C_dldd_clogit,   5},
  {"C_dldd_prophaz",   (DL_FUNC) &C_dldd_prophaz,  6},
  {"C_compute_ERCmatrix_clogit",   (DL_FUNC) &C_compute_ERCmatrix_clogit,  8},
  {"C_compute_ERCmatrix_prophaz",  (DL_FUNC) &C_compute_ERCmatrix_prophaz, 8},
  {"C_loglik_prophaz_rcpp",     (DL_FUNC) &C_loglik_prophaz_rcpp,    9},
  {NULL, NULL, 0}
};

void R_init_ameras(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}


