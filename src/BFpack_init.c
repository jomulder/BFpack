#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(compute_rcet)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(compute_rcet2)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(draw_ju)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(estimate_postmeancov_fisherz)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"compute_rcet",                 (DL_FUNC) &F77_NAME(compute_rcet),                  6},
    {"compute_rcet2",                (DL_FUNC) &F77_NAME(compute_rcet2),                 9},
    {"draw_ju",                      (DL_FUNC) &F77_NAME(draw_ju),                       6},
    {"estimate_postmeancov_fisherz", (DL_FUNC) &F77_NAME(estimate_postmeancov_fisherz), 22},
    {NULL, NULL, 0}
};

void R_init_BFpack(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

