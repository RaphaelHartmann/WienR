
// Chair of Social Psychology, University of Freiburg
// Authors: Raphael Hartmann and Christoph Klauer

#include "tools.h"
#include "derivs.h"
#include "methods.h"
#include <Rinternals.h>

extern "C" {

	SEXP dWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *w = REAL(re4);
		double *sv = REAL(re5);
		double eps = REAL(re6)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP pdf = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP logpdf = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rpdf = REAL(pdf);
		double *Rlogpdf = REAL(logpdf);

		/* calculate the derivatives */
		PDF(t, a, v, w, sv, eps, resp, K, N, epsFLAG, Rpdf, Rlogpdf, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,pdf);
		SET_VECTOR_ELT(out,1,logpdf);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("pdf"));
		SET_STRING_ELT(names,1,Rf_mkChar("logpdf"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}




extern "C" {

	SEXP dtdWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *w = REAL(re4);
		double *sv = REAL(re5);
		double eps = REAL(re6)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP deriv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rderiv = REAL(deriv);


		/* calculate the derivatives */
		dtPDF(t, a, v, w, sv, eps, resp, K, N, epsFLAG, Rderiv, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,deriv);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("deriv"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}




extern "C" {

	SEXP dadWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *w = REAL(re4);
		double *sv = REAL(re5);
		double eps = REAL(re6)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP deriv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rderiv = REAL(deriv);


		/* calculate the derivatives */
		daPDF(t, a, v, w, sv, eps, resp, K, N, epsFLAG, Rderiv, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,deriv);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("deriv"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}




extern "C" {

	SEXP dvdWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *w = REAL(re4);
		double *sv = REAL(re5);
		double eps = REAL(re6)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP deriv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rderiv = REAL(deriv);


		/* calculate the derivatives */
		dvPDF(t, a, v, w, sv, eps, resp, K, N, epsFLAG, Rderiv, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,deriv);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("deriv"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}




extern "C" {

	SEXP dwdWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *w = REAL(re4);
		double *sv = REAL(re5);
		double eps = REAL(re6)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP deriv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rderiv = REAL(deriv);


		/* calculate the derivatives */
		dwPDF(t, a, v, w, sv, eps, resp, K, N, epsFLAG, Rderiv, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,deriv);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("deriv"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}




extern "C" {
  
  SEXP dsvdWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {
    
    /* define input variables */
    double *t = REAL(re1);
    double *a = REAL(re2);
    double *v = REAL(re3);
    double *w = REAL(re4);
    double *sv = REAL(re5);
    double eps = REAL(re6)[0];
    
    int *resp = INTEGER(in1);
    int K = INTEGER(in2)[0];
    int N = INTEGER(in3)[0];
    int NThreads = INTEGER(in4)[0];
    
    int epsFLAG = INTEGER(bo1)[0];
    
    
    /* declare R objects for output */
    int outCnt = 0, prtCnt = 0;
    SEXP deriv = PROTECT(Rf_allocVector(REALSXP, N));
    outCnt++;
    SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
    prtCnt = outCnt + 1;
    
    
    /* declare C++ pointers for R objects */
    double *Rderiv = REAL(deriv);
    
    
    /* calculate the derivatives */
    dsvPDF(t, a, v, w, sv, eps, resp, K, N, epsFLAG, Rderiv, NThreads);
    
    
    /* set elements of list out */
    SET_VECTOR_ELT(out,0,deriv);
    
    
    /* make name vector and set element names */
    SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
    prtCnt++;
    SET_STRING_ELT(names,0,Rf_mkChar("deriv"));
    
    Rf_setAttrib(out,R_NamesSymbol,names);
    
    
    /* Unprotect the out and names objects */
    UNPROTECT(prtCnt);
    
    return(out);
  }
  
}



extern "C" {

	SEXP pWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *w = REAL(re4);
		double eps = REAL(re5)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP cdf = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP logcdf = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rcdf = REAL(cdf);
		double *Rlogcdf = REAL(logcdf);


		/* calculate the derivatives */
		CDF(t, a, v, w, eps, resp, K, N, epsFLAG, Rcdf, Rlogcdf, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,cdf);
		SET_VECTOR_ELT(out,1,logcdf);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("cdf"));
		SET_STRING_ELT(names,1,Rf_mkChar("logcdf"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}




extern "C" {

	SEXP dapWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *w = REAL(re4);
		double eps = REAL(re5)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP deriv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rderiv = REAL(deriv);

		/* calculate the derivatives */
		daCDF(t, a, v, w, eps, resp, K, N, epsFLAG, Rderiv, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,deriv);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("deriv"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}




extern "C" {

	SEXP dvpWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *w = REAL(re4);
		double eps = REAL(re5)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP deriv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rderiv = REAL(deriv);


		/* calculate the derivatives */
		dvCDF(t, a, v, w, eps, resp, K, N, epsFLAG, Rderiv, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,deriv);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("deriv"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}




extern "C" {

	SEXP dwpWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *w = REAL(re4);
		double eps = REAL(re5)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP deriv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rderiv = REAL(deriv);


		/* calculate the derivatives */
		dwCDF(t, a, v, w, eps, resp, K, N, epsFLAG, Rderiv, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,deriv);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("deriv"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}




extern "C" {

	SEXP dxdWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *w = REAL(re4);
		double *sv = REAL(re5);
		double eps = REAL(re6)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP da = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dw = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rda = REAL(da);
		double *Rdv = REAL(dv);
		double *Rdw = REAL(dw);


		/* calculate the derivatives */
		dxPDF(t, a, v, w, sv, eps, resp, K, N, epsFLAG, Rda, Rdv, Rdw, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,da);
		SET_VECTOR_ELT(out,1,dv);
		SET_VECTOR_ELT(out,2,dw);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("da"));
		SET_STRING_ELT(names,1,Rf_mkChar("dv"));
		SET_STRING_ELT(names,2,Rf_mkChar("dw"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}




extern "C" {

	SEXP dxpWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *w = REAL(re4);
		double eps = REAL(re5)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP da = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dw = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rda = REAL(da);
		double *Rdv = REAL(dv);
		double *Rdw = REAL(dw);


		/* calculate the derivatives */
		dxCDF(t, a, v, w, eps, resp, K, N, epsFLAG, Rda, Rdv, Rdw, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,da);
		SET_VECTOR_ELT(out,1,dv);
		SET_VECTOR_ELT(out,2,dw);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("da"));
		SET_STRING_ELT(names,1,Rf_mkChar("dv"));
		SET_STRING_ELT(names,2,Rf_mkChar("dw"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------7-PARAM DIFFUSION-------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

extern "C" {

	SEXP dDiffusion7(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP re7, SEXP re8, SEXP re9, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP in5, SEXP in6, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *to = REAL(re4);
		double *w = REAL(re5);
		double *sw = REAL(re6);
		double *sv = REAL(re7);
		double *st = REAL(re8);
		double eps = REAL(re9)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];
		int choice = INTEGER(in5)[0];
		int Neval = INTEGER(in6)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP value = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP value_ln = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP err = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rvalue = REAL(value);
		double *Rvalue_ln = REAL(value_ln);
		double *Rerr = REAL(err);


		/* calculate the PDF or derivatives */
		PDF7(choice, t, resp, a, v, to, w, sw, sv, st, eps, K, N, epsFLAG, Rvalue, Rvalue_ln, Rerr, NThreads, Neval);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,value);
		SET_VECTOR_ELT(out,1,value_ln);
		SET_VECTOR_ELT(out,2,err);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		if (choice == 0) {
			SET_STRING_ELT(names,0,Rf_mkChar("pdf"));
			SET_STRING_ELT(names,1,Rf_mkChar("logpdf"));
			SET_STRING_ELT(names,2,Rf_mkChar("err"));
		} else {
			Rvalue_ln = nullptr;
			SET_STRING_ELT(names,0,Rf_mkChar("deriv"));
			SET_STRING_ELT(names,1,Rf_mkChar("deriv_ln"));
			SET_STRING_ELT(names,2,Rf_mkChar("err"));
		}


		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}


extern "C" {

	SEXP pDiffusion7(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP re7, SEXP re8, SEXP re9, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP in5, SEXP in6, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *to = REAL(re4);
		double *w = REAL(re5);
		double *sw = REAL(re6);
		double *sv = REAL(re7);
		double *st = REAL(re8);
		double eps = REAL(re9)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];
		int choice = INTEGER(in5)[0];
		int Neval = INTEGER(in6)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP value = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP value_ln = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP err = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rvalue = REAL(value);
		double *Rvalue_ln = REAL(value_ln);
		double *Rerr = REAL(err);


		/* calculate the PDF or derivatives */
		CDF7(choice, t, resp, a, v, to, w, sw, sv, st, eps, K, N, epsFLAG, Rvalue, Rvalue_ln, Rerr, NThreads, Neval);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,value);
		SET_VECTOR_ELT(out,1,value_ln);
		SET_VECTOR_ELT(out,2,err);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		if (choice == 0) {
			SET_STRING_ELT(names,0,Rf_mkChar("cdf"));
			SET_STRING_ELT(names,1,Rf_mkChar("logcdf"));
			SET_STRING_ELT(names,2,Rf_mkChar("err"));
		} else {
			Rvalue_ln = nullptr;
			SET_STRING_ELT(names,0,Rf_mkChar("deriv"));
			SET_STRING_ELT(names,1,Rf_mkChar("deriv_ln"));
			SET_STRING_ELT(names,2,Rf_mkChar("err"));
		}

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}





extern "C" {

	SEXP dxdDiffusion7(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP re7, SEXP re8, SEXP re9, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP in5, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *t0 = REAL(re4);
		double *w = REAL(re5);
		double *sw = REAL(re6);
		double *sv = REAL(re7);
		double *st = REAL(re8);
		double eps = REAL(re9)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];
		int Neval = INTEGER(in5)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP da = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dt0 = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dw = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dsw = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dsv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dst = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP err = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rda = REAL(da);
		double *Rdv = REAL(dv);
		double *Rdt0 = REAL(dt0);
		double *Rdw = REAL(dw);
		double *Rdsw = REAL(dsw);
		double *Rdsv = REAL(dsv);
		double *Rdst = REAL(dst);
		double *Rerr = REAL(err);


		/* calculate the derivatives */
		dxPDF7(t, resp, a, v, t0, w, sw, sv, st, eps, K, N, epsFLAG, Rda, Rdv, Rdt0, Rdw, Rdsw, Rdsv, Rdst, Rerr, NThreads, Neval);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,da);
		SET_VECTOR_ELT(out,1,dv);
		SET_VECTOR_ELT(out,2,dt0);
		SET_VECTOR_ELT(out,3,dw);
		SET_VECTOR_ELT(out,4,dsw);
		SET_VECTOR_ELT(out,5,dsv);
		SET_VECTOR_ELT(out,6,dst);
		SET_VECTOR_ELT(out,7,err);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("da"));
		SET_STRING_ELT(names,1,Rf_mkChar("dv"));
		SET_STRING_ELT(names,2,Rf_mkChar("dt0"));
		SET_STRING_ELT(names,3,Rf_mkChar("dw"));
		SET_STRING_ELT(names,4,Rf_mkChar("dsw"));
		SET_STRING_ELT(names,5,Rf_mkChar("dsv"));
		SET_STRING_ELT(names,6,Rf_mkChar("dst"));
		SET_STRING_ELT(names,7,Rf_mkChar("err"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}




extern "C" {

	SEXP dxpDiffusion7(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP re7, SEXP re8, SEXP re9, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP in5, SEXP bo1) {

		/* define input variables */
		double *t = REAL(re1);
		double *a = REAL(re2);
		double *v = REAL(re3);
		double *t0 = REAL(re4);
		double *w = REAL(re5);
		double *sw = REAL(re6);
		double *sv = REAL(re7);
		double *st = REAL(re8);
		double eps = REAL(re9)[0];

		int *resp = INTEGER(in1);
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];
		int Neval = INTEGER(in5)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP da = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dt0 = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dw = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dsw = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dsv = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP dst = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP err = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rda = REAL(da);
		double *Rdv = REAL(dv);
		double *Rdt0 = REAL(dt0);
		double *Rdw = REAL(dw);
		double *Rdsw = REAL(dsw);
		double *Rdsv = REAL(dsv);
		double *Rdst = REAL(dst);
		double *Rerr = REAL(err);


		/* calculate the derivatives */
		dxCDF7(t, resp, a, v, t0, w, sw, sv, st, eps, K, N, epsFLAG, Rda, Rdv, Rdt0, Rdw, Rdsw, Rdsv, Rdst, Rerr, NThreads, Neval);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,da);
		SET_VECTOR_ELT(out,1,dv);
		SET_VECTOR_ELT(out,2,dt0);
		SET_VECTOR_ELT(out,3,dw);
		SET_VECTOR_ELT(out,4,dsw);
		SET_VECTOR_ELT(out,5,dsv);
		SET_VECTOR_ELT(out,6,dst);
		SET_VECTOR_ELT(out,7,err);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("da"));
		SET_STRING_ELT(names,1,Rf_mkChar("dv"));
		SET_STRING_ELT(names,2,Rf_mkChar("dt0"));
		SET_STRING_ELT(names,3,Rf_mkChar("dw"));
		SET_STRING_ELT(names,4,Rf_mkChar("dsw"));
		SET_STRING_ELT(names,5,Rf_mkChar("dsv"));
		SET_STRING_ELT(names,6,Rf_mkChar("dst"));
		SET_STRING_ELT(names,7,Rf_mkChar("err"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}



/* SAMPLING */

extern "C" {

	SEXP randWiener(SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP re7, SEXP re8, SEXP re9, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP in5, SEXP in6, SEXP bo1, SEXP bo2) {

		/* define input variables */
		double a = REAL(re2)[0];
		double v = REAL(re3)[0];
		double w = REAL(re4)[0];
		double sv = REAL(re5)[0];
		double sw = REAL(re6)[0];
		double eps = REAL(re7)[0];
		double bound = REAL(re8)[0];
		double *list = REAL(re9);

		int R = INTEGER(in1)[0];
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];
		int choice = INTEGER(in5)[0];
		int *list_nrs = INTEGER(in6);

		int epsFLAG = INTEGER(bo1)[0];
		int storeFLAG = INTEGER(bo2)[0];

		if (choice != 1) storeFLAG = 0;


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP q = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP resp = PROTECT(Rf_allocVector(INTSXP, N));
		outCnt++;
		SEXP store1 = PROTECT(Rf_allocVector(VECSXP, 15));
		prtCnt++;
		SEXP store2 = PROTECT(Rf_allocVector(VECSXP, 15));
		prtCnt++;
		int add = storeFLAG ? (R==0 ? 2 : 1) : 0;

		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt+add));
		prtCnt += outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rq = REAL(q);
		int *Rresp = INTEGER(resp);
		int use_store = 0;

		/* convert ars_archiv if needed */
		ars_archiv ars_arch1, ars_arch2;
		if (list_nrs[0] != 0) {

			ars_archiv temp_arch1;
			std::vector<point> hstr(list_nrs[1]);
			std::vector<piece> lowstr(list_nrs[2]);
			std::vector<piece> uppstr(list_nrs[3]);
			std::vector<double> sstr(list_nrs[4]);
			int i;
			for (i = 0; i < list_nrs[1]; ++i) {
				hstr[i].x = list[i];
				hstr[i].h = list[i + list_nrs[1]];
				hstr[i].dh = list[i + list_nrs[1] + list_nrs[1]];
			}
			temp_arch1.hstore.assign(hstr.begin(), hstr.end());
			int offset = 3 * list_nrs[1];
			for (i = 0; i < list_nrs[2]; ++i) {
				lowstr[i].z = list[i + offset];
				lowstr[i].slope = list[i + offset + list_nrs[2]];
				lowstr[i].absc = list[i + offset + 2*list_nrs[2]];
				lowstr[i].center = list[i + offset + 3*list_nrs[2]];
			}
			temp_arch1.lowerstore.assign(lowstr.begin(), lowstr.end());
			offset += 4*list_nrs[2];
			for (i = 0; i < list_nrs[3]; ++i) {
				uppstr[i].z = list[i + offset];
				uppstr[i].slope = list[i + offset + list_nrs[3]];
				uppstr[i].absc = list[i + offset + 2*list_nrs[3]];
				uppstr[i].center = list[i + offset + 3*list_nrs[3]];
			}
			temp_arch1.upperstore.assign(uppstr.begin(), uppstr.end());
			offset += 4*list_nrs[3];
			temp_arch1.startstore = list[offset++];
			temp_arch1.scalestore = list[offset++];
			temp_arch1.normstore = list[offset++];
			for (i = 0; i < list_nrs[4]; ++i) {
				sstr[i] = list[i + offset];
			}
			temp_arch1.sstore.assign(sstr.begin(), sstr.end());
			offset += list_nrs[4];

			if (list_nrs[0] == 2) {

				ars_archiv temp_arch2;
				std::vector<point> hstr2(list_nrs[5]);
				std::vector<piece> lowstr2(list_nrs[6]);
				std::vector<piece> uppstr2(list_nrs[7]);
				std::vector<double> sstr2(list_nrs[8]);
				int i;
				for (i = 0; i < list_nrs[5]; ++i) {
					hstr2[i].x = list[i + offset];
					hstr2[i].h = list[i + offset + list_nrs[5]];
					hstr2[i].dh = list[i + offset + list_nrs[5] + list_nrs[5]];
				}
				temp_arch2.hstore.assign(hstr2.begin(), hstr2.end());
				offset += 3 * list_nrs[5];
				for (i = 0; i < list_nrs[6]; ++i) {
					lowstr2[i].z = list[i + offset];
					lowstr2[i].slope = list[i + offset + list_nrs[6]];
					lowstr2[i].absc = list[i + offset + 2*list_nrs[6]];
					lowstr2[i].center = list[i + offset + 3*list_nrs[6]];
				}
				temp_arch2.lowerstore.assign(lowstr2.begin(), lowstr2.end());
				offset += 4*list_nrs[6];
				for (i = 0; i < list_nrs[7]; ++i) {
					uppstr2[i].z = list[i + offset];
					uppstr2[i].slope = list[i + offset + list_nrs[7]];
					uppstr2[i].absc = list[i + offset + 2*list_nrs[7]];
					uppstr2[i].center = list[i + offset + 3*list_nrs[7]];
				}
				temp_arch2.upperstore.assign(uppstr2.begin(), uppstr2.end());
				offset += 4*list_nrs[7];
				temp_arch2.startstore = list[offset++];
				temp_arch2.scalestore = list[offset++];
				temp_arch2.normstore = list[offset++];
				for (i = 0; i < list_nrs[8]; ++i) {
					sstr2[i] = list[i + offset];
				}
				temp_arch2.sstore.assign(sstr2.begin(), sstr2.end());

				if (R == 0) {
					ars_arch1 = temp_arch1;
					ars_arch2 = temp_arch2;
					use_store = 1;
				} else if (R == 1) {
					ars_arch1 = temp_arch2;
					use_store = 1;
				} else if (R == 2) {
					ars_arch1 = temp_arch1;
					use_store = 1;
				}

			} else {
				if (R != 0) {
					ars_arch1 = temp_arch1;
					use_store = 1;
				}
			}

		}


		/* sampling */
		run_make_rwiener(choice, N, a, v, w, sv, sw, R, bound, eps, K, epsFLAG, NThreads, Rq, Rresp, &ars_arch1, &ars_arch2, use_store);


		/* fill list store */
		if (storeFLAG) {
			int hstore1_size = static_cast<int>(ars_arch1.hstore.size());
			int lowerstore1_size = static_cast<int>(ars_arch1.lowerstore.size());
			int upperstore1_size = static_cast<int>(ars_arch1.upperstore.size());
			int sstore1_size = static_cast<int>(ars_arch1.sstore.size());

			SEXP hstore1_x = PROTECT(Rf_allocVector(REALSXP, hstore1_size));
			SEXP hstore1_h = PROTECT(Rf_allocVector(REALSXP, hstore1_size));
			SEXP hstore1_dh = PROTECT(Rf_allocVector(REALSXP, hstore1_size));
			SEXP lowerstore1_z = PROTECT(Rf_allocVector(REALSXP, lowerstore1_size));
			SEXP lowerstore1_slope = PROTECT(Rf_allocVector(REALSXP, lowerstore1_size));
			SEXP lowerstore1_absc = PROTECT(Rf_allocVector(REALSXP, lowerstore1_size));
			SEXP lowerstore1_center = PROTECT(Rf_allocVector(REALSXP, lowerstore1_size));
			SEXP upperstore1_z = PROTECT(Rf_allocVector(REALSXP, upperstore1_size));
			SEXP upperstore1_slope = PROTECT(Rf_allocVector(REALSXP, upperstore1_size));
			SEXP upperstore1_absc = PROTECT(Rf_allocVector(REALSXP, upperstore1_size));
			SEXP upperstore1_center = PROTECT(Rf_allocVector(REALSXP, upperstore1_size));
			SEXP startstore1 = PROTECT(Rf_allocVector(REALSXP, 1));
			SEXP scalestore1 = PROTECT(Rf_allocVector(REALSXP, 1));
			SEXP normstore1 = PROTECT(Rf_allocVector(REALSXP, 1));
			SEXP sstore1 = PROTECT(Rf_allocVector(REALSXP, sstore1_size));
			prtCnt += 15;

			double *Rhstore1_x = REAL(hstore1_x);
			double *Rhstore1_h = REAL(hstore1_h);
			double *Rhstore1_dh = REAL(hstore1_dh);
			double *Rlowerstore1_z = REAL(lowerstore1_z);
			double *Rlowerstore1_slope = REAL(lowerstore1_slope);
			double *Rlowerstore1_absc = REAL(lowerstore1_absc);
			double *Rlowerstore1_center = REAL(lowerstore1_center);
			double *Rupperstore1_z = REAL(upperstore1_z);
			double *Rupperstore1_slope = REAL(upperstore1_slope);
			double *Rupperstore1_absc = REAL(upperstore1_absc);
			double *Rupperstore1_center = REAL(upperstore1_center);
			double *Rstartstore1 = REAL(startstore1);
			double *Rscalestore1 = REAL(scalestore1);
			double *Rnormstore1 = REAL(normstore1);
			double *Rsstore1 = REAL(sstore1);

			SET_VECTOR_ELT(store1, 0,hstore1_x);
			SET_VECTOR_ELT(store1, 1,hstore1_h);
			SET_VECTOR_ELT(store1, 2,hstore1_dh);
			SET_VECTOR_ELT(store1, 3,lowerstore1_z);
			SET_VECTOR_ELT(store1, 4,lowerstore1_slope);
			SET_VECTOR_ELT(store1, 5,lowerstore1_absc);
			SET_VECTOR_ELT(store1, 6,lowerstore1_center);
			SET_VECTOR_ELT(store1, 7,upperstore1_z);
			SET_VECTOR_ELT(store1, 8,upperstore1_slope);
			SET_VECTOR_ELT(store1, 9,upperstore1_absc);
			SET_VECTOR_ELT(store1,10,upperstore1_center);
			SET_VECTOR_ELT(store1,11,startstore1);
			SET_VECTOR_ELT(store1,12,scalestore1);
			SET_VECTOR_ELT(store1,13,normstore1);
			SET_VECTOR_ELT(store1,14,sstore1);

			for (int i=0; i < hstore1_size; ++i) {
				Rhstore1_x[i] = ars_arch1.hstore[i].x;
				Rhstore1_h[i] = ars_arch1.hstore[i].h;
				Rhstore1_dh[i] = ars_arch1.hstore[i].dh;
			}
			for (int i=0; i < lowerstore1_size; ++i) {
				Rlowerstore1_z[i] = ars_arch1.lowerstore[i].z;
				Rlowerstore1_slope[i] = ars_arch1.lowerstore[i].slope;
				Rlowerstore1_absc[i] = ars_arch1.lowerstore[i].absc;
				Rlowerstore1_center[i] = ars_arch1.lowerstore[i].center;
			}
			for (int i=0; i < upperstore1_size; ++i) {
				Rupperstore1_z[i] = ars_arch1.upperstore[i].z;
				Rupperstore1_slope[i] = ars_arch1.upperstore[i].slope;
				Rupperstore1_absc[i] = ars_arch1.upperstore[i].absc;
				Rupperstore1_center[i] = ars_arch1.upperstore[i].center;
			}
			Rstartstore1[0] = ars_arch1.startstore;
			Rscalestore1[0] = ars_arch1.scalestore;
			Rnormstore1[0] = ars_arch1.normstore;
			for (int i=0; i < sstore1_size; ++i) {
				Rsstore1[i] = ars_arch1.sstore[i];
			}

			if (R == 0) {
				int hstore2_size = static_cast<int>(ars_arch2.hstore.size());
				int lowerstore2_size = static_cast<int>(ars_arch2.lowerstore.size());
				int upperstore2_size = static_cast<int>(ars_arch2.upperstore.size());
				int sstore2_size = static_cast<int>(ars_arch2.sstore.size());

				SEXP hstore2_x = PROTECT(Rf_allocVector(REALSXP, hstore2_size));
				SEXP hstore2_h = PROTECT(Rf_allocVector(REALSXP, hstore2_size));
				SEXP hstore2_dh = PROTECT(Rf_allocVector(REALSXP, hstore2_size));
				SEXP lowerstore2_z = PROTECT(Rf_allocVector(REALSXP, lowerstore2_size));
				SEXP lowerstore2_slope = PROTECT(Rf_allocVector(REALSXP, lowerstore2_size));
				SEXP lowerstore2_absc = PROTECT(Rf_allocVector(REALSXP, lowerstore2_size));
				SEXP lowerstore2_center = PROTECT(Rf_allocVector(REALSXP, lowerstore2_size));
				SEXP upperstore2_z = PROTECT(Rf_allocVector(REALSXP, upperstore2_size));
				SEXP upperstore2_slope = PROTECT(Rf_allocVector(REALSXP, upperstore2_size));
				SEXP upperstore2_absc = PROTECT(Rf_allocVector(REALSXP, upperstore2_size));
				SEXP upperstore2_center = PROTECT(Rf_allocVector(REALSXP, upperstore2_size));
				SEXP startstore2 = PROTECT(Rf_allocVector(REALSXP, 1));
				SEXP scalestore2 = PROTECT(Rf_allocVector(REALSXP, 1));
				SEXP normstore2 = PROTECT(Rf_allocVector(REALSXP, 1));
				SEXP sstore2 = PROTECT(Rf_allocVector(REALSXP, sstore2_size));
				prtCnt += 15;

				double *Rhstore2_x = REAL(hstore2_x);
				double *Rhstore2_h = REAL(hstore2_h);
				double *Rhstore2_dh = REAL(hstore2_dh);
				double *Rlowerstore2_z = REAL(lowerstore2_z);
				double *Rlowerstore2_slope = REAL(lowerstore2_slope);
				double *Rlowerstore2_absc = REAL(lowerstore2_absc);
				double *Rlowerstore2_center = REAL(lowerstore2_center);
				double *Rupperstore2_z = REAL(upperstore2_z);
				double *Rupperstore2_slope = REAL(upperstore2_slope);
				double *Rupperstore2_absc = REAL(upperstore2_absc);
				double *Rupperstore2_center = REAL(upperstore2_center);
				double *Rstartstore2 = REAL(startstore2);
				double *Rscalestore2 = REAL(scalestore2);
				double *Rnormstore2 = REAL(normstore2);
				double *Rsstore2 = REAL(sstore2);

				SET_VECTOR_ELT(store2, 0,hstore2_x);
				SET_VECTOR_ELT(store2, 1,hstore2_h);
				SET_VECTOR_ELT(store2, 2,hstore2_dh);
				SET_VECTOR_ELT(store2, 3,lowerstore2_z);
				SET_VECTOR_ELT(store2, 4,lowerstore2_slope);
				SET_VECTOR_ELT(store2, 5,lowerstore2_absc);
				SET_VECTOR_ELT(store2, 6,lowerstore2_center);
				SET_VECTOR_ELT(store2, 7,upperstore2_z);
				SET_VECTOR_ELT(store2, 8,upperstore2_slope);
				SET_VECTOR_ELT(store2, 9,upperstore2_absc);
				SET_VECTOR_ELT(store2,10,upperstore2_center);
				SET_VECTOR_ELT(store2,11,startstore2);
				SET_VECTOR_ELT(store2,12,scalestore2);
				SET_VECTOR_ELT(store2,13,normstore2);
				SET_VECTOR_ELT(store2,14,sstore2);

				for (int i=0; i < hstore2_size; ++i) {
					Rhstore2_x[i] = ars_arch2.hstore[i].x;
					Rhstore2_h[i] = ars_arch2.hstore[i].h;
					Rhstore2_dh[i] = ars_arch2.hstore[i].dh;
				}
				for (int i=0; i < lowerstore2_size; ++i) {
					Rlowerstore2_z[i] = ars_arch2.lowerstore[i].z;
					Rlowerstore2_slope[i] = ars_arch2.lowerstore[i].slope;
					Rlowerstore2_absc[i] = ars_arch2.lowerstore[i].absc;
					Rlowerstore2_center[i] = ars_arch2.lowerstore[i].center;
				}
				for (int i=0; i < upperstore2_size; ++i) {
					Rupperstore2_z[i] = ars_arch2.upperstore[i].z;
					Rupperstore2_slope[i] = ars_arch2.upperstore[i].slope;
					Rupperstore2_absc[i] = ars_arch2.upperstore[i].absc;
					Rupperstore2_center[i] = ars_arch2.upperstore[i].center;
				}
				Rstartstore2[0] = ars_arch2.startstore;
				Rscalestore2[0] = ars_arch2.scalestore;
				Rnormstore2[0] = ars_arch2.normstore;
				for (int i=0; i < sstore2_size; ++i) {
					Rsstore2[i] = ars_arch2.sstore[i];
				}
			}

			SEXP store_names = PROTECT(Rf_allocVector(STRSXP, 15));
			prtCnt++;
			SET_STRING_ELT(store_names, 0,Rf_mkChar("hstore_x"));
			SET_STRING_ELT(store_names, 1,Rf_mkChar("hstore_h"));
			SET_STRING_ELT(store_names, 2,Rf_mkChar("hstore_dh"));
			SET_STRING_ELT(store_names, 3,Rf_mkChar("lowerstore_z"));
			SET_STRING_ELT(store_names, 4,Rf_mkChar("lowerstore_slope"));
			SET_STRING_ELT(store_names, 5,Rf_mkChar("lowerstore_absc"));
			SET_STRING_ELT(store_names, 6,Rf_mkChar("lowerstore_center"));
			SET_STRING_ELT(store_names, 7,Rf_mkChar("upperstore_z"));
			SET_STRING_ELT(store_names, 8,Rf_mkChar("upperstore_slope"));
			SET_STRING_ELT(store_names, 9,Rf_mkChar("upperstore_absc"));
			SET_STRING_ELT(store_names,10,Rf_mkChar("upperstore_center"));
			SET_STRING_ELT(store_names,11,Rf_mkChar("startstore"));
			SET_STRING_ELT(store_names,12,Rf_mkChar("scalestore"));
			SET_STRING_ELT(store_names,13,Rf_mkChar("normstore"));
			SET_STRING_ELT(store_names,14,Rf_mkChar("sstore"));

			Rf_setAttrib(store1,R_NamesSymbol,store_names);
			if (R == 0) Rf_setAttrib(store2,R_NamesSymbol,store_names);
		}


		/* set elements of list(s) store and give names */


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,q);
		SET_VECTOR_ELT(out,1,resp);
		if (storeFLAG) {
			SET_VECTOR_ELT(out,2,store1);
			if (R == 0) SET_VECTOR_ELT(out,3,store2);
		}


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt+add));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("q"));
		SET_STRING_ELT(names,1,Rf_mkChar("resp"));
		if (storeFLAG) {
			if (R != 0) SET_STRING_ELT(names,2,Rf_mkChar("ars.store"));
			if (R == 0) {
				SET_STRING_ELT(names,2,Rf_mkChar("ars.store.upp"));
				SET_STRING_ELT(names,3,Rf_mkChar("ars.store.low"));
			}
		}
		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}
