#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(bivrecur)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mprovar)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(onesamp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xmproee)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ymproee)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"bivrecur",  (DL_FUNC) &F77_NAME(bivrecur),  27},
    {"mprovar",   (DL_FUNC) &F77_NAME(mprovar),   16},
    {"onesamp",   (DL_FUNC) &F77_NAME(onesamp),   17},
    {"xmproee",   (DL_FUNC) &F77_NAME(xmproee),    9},
    {"ymproee",   (DL_FUNC) &F77_NAME(ymproee),   10},
    {NULL, NULL, 0}
};

void R_init_BivRec(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}