#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP R_corHuberAdj(SEXP, SEXP, SEXP);
extern SEXP R_corHuberBi(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_corHuberUni(SEXP, SEXP, SEXP);
extern SEXP R_corMatHuber(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_fastGrplars(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_fastLars(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_fastLasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_fastSparseLTS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_findSmallest(SEXP, SEXP);
extern SEXP R_partialOrder(SEXP, SEXP);
extern SEXP R_sparseSubsets(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_testCStep(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_testKeepBest(SEXP, SEXP, SEXP);
extern SEXP R_testLasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"R_corHuberAdj",   (DL_FUNC) &R_corHuberAdj,    3},
  {"R_corHuberBi",    (DL_FUNC) &R_corHuberBi,     5},
  {"R_corHuberUni",   (DL_FUNC) &R_corHuberUni,    3},
  {"R_corMatHuber",   (DL_FUNC) &R_corMatHuber,    4},
  {"R_fastGrplars",   (DL_FUNC) &R_fastGrplars,    5},
  {"R_fastLars",      (DL_FUNC) &R_fastLars,       9},
  {"R_fastLasso",     (DL_FUNC) &R_fastLasso,      9},
  {"R_fastSparseLTS", (DL_FUNC) &R_fastSparseLTS, 12},
  {"R_findSmallest",  (DL_FUNC) &R_findSmallest,   2},
  {"R_partialOrder",  (DL_FUNC) &R_partialOrder,   2},
  {"R_sparseSubsets", (DL_FUNC) &R_sparseSubsets,  9},
  {"R_testCStep",     (DL_FUNC) &R_testCStep,      9},
  {"R_testKeepBest",  (DL_FUNC) &R_testKeepBest,   3},
  {"R_testLasso",     (DL_FUNC) &R_testLasso,      8},
  {NULL, NULL, 0}
};

void R_init_robustHD(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
