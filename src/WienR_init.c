
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP dWiener(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dtdWiener(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dadWiener(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dvdWiener(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dwdWiener(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pWiener(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dapWiener(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dvpWiener(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dwpWiener(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP dxdWiener(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dxpWiener(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP dDiffusion7(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pDiffusion7(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP dxdDiffusion7(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dxpDiffusion7(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP randWiener(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    {"dWiener", (DL_FUNC) &dWiener, 10},
    {"dtdWiener", (DL_FUNC) &dtdWiener, 10},
    {"dadWiener", (DL_FUNC) &dadWiener, 10},
    {"dvdWiener", (DL_FUNC) &dvdWiener, 10},
    {"dwdWiener", (DL_FUNC) &dwdWiener, 10},
    {"pWiener", (DL_FUNC) &pWiener, 10},
    {"dapWiener", (DL_FUNC) &dapWiener, 10},
    {"dvpWiener", (DL_FUNC) &dvpWiener, 10},
    {"dwpWiener", (DL_FUNC) &dwpWiener, 10},

    {"dxdWiener", (DL_FUNC) &dxdWiener, 10},
    {"dxpWiener", (DL_FUNC) &dxpWiener, 10},

    {"dDiffusion7", (DL_FUNC) &dDiffusion7, 15},
    {"pDiffusion7", (DL_FUNC) &pDiffusion7, 15},

    {"dxdDiffusion7", (DL_FUNC) &dxdDiffusion7, 14},
    {"dxpDiffusion7", (DL_FUNC) &dxpDiffusion7, 14},

    {"randWiener", (DL_FUNC) &randWiener, 16},

    {NULL, NULL, 0}
};

void R_init_WienR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
