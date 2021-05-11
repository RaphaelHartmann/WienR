
// Chair of Social Psychology, University of Freiburg
// Authors: Raphael Hartmann and Christoph Klauer

#include "derivs.h"
#include <Rinternals.h>
#include "rwiener.h"

extern "C" {

	SEXP dWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

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
		PDF(t, a, v, w, eps, resp, K, N, epsFLAG, Rpdf, Rlogpdf, NThreads);


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

	SEXP dtdWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

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
		dtPDF(t, a, v, w, eps, resp, K, N, epsFLAG, Rderiv, NThreads);


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

	SEXP dadWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

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
		daPDF(t, a, v, w, eps, resp, K, N, epsFLAG, Rderiv, NThreads);


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

	SEXP dvdWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

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
		dvPDF(t, a, v, w, eps, resp, K, N, epsFLAG, Rderiv, NThreads);


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

	SEXP dwdWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

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
		dwPDF(t, a, v, w, eps, resp, K, N, epsFLAG, Rderiv, NThreads);


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

	SEXP dxdWiener(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

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
		dxPDF(t, a, v, w, eps, resp, K, N, epsFLAG, Rda, Rdv, Rdw, NThreads);


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

	SEXP dDiffusion7(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP re7, SEXP re8, SEXP re9, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP in5, SEXP bo1) {

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

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP value = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP value_ln = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rvalue = REAL(value);
		double *Rvalue_ln = REAL(value_ln);


		/* calculate the PDF or derivatives */
		PDF7(choice, t, resp, a, v, to, w, sw, sv, st, eps, K, N, epsFLAG, Rvalue, Rvalue_ln, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,value);
		SET_VECTOR_ELT(out,1,value_ln);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		if (choice == 0) {
			SET_STRING_ELT(names,0,Rf_mkChar("pdf"));
			SET_STRING_ELT(names,1,Rf_mkChar("logpdf"));
		} else {
			Rvalue_ln = nullptr;
			SET_STRING_ELT(names,0,Rf_mkChar("deriv"));
			SET_STRING_ELT(names,1,Rf_mkChar("deriv_ln"));
		}


		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}


extern "C" {

	SEXP pDiffusion7(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP re7, SEXP re8, SEXP re9, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP in5, SEXP bo1) {

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

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP value = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP value_ln = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rvalue = REAL(value);
		double *Rvalue_ln = REAL(value_ln);


		/* calculate the PDF or derivatives */
		CDF7(choice, t, resp, a, v, to, w, sw, sv, st, eps, K, N, epsFLAG, Rvalue, Rvalue_ln, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,value);
		SET_VECTOR_ELT(out,1,value_ln);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		if (choice == 0) {
			SET_STRING_ELT(names,0,Rf_mkChar("cdf"));
			SET_STRING_ELT(names,1,Rf_mkChar("logcdf"));
		} else {
			Rvalue_ln = nullptr;
			SET_STRING_ELT(names,0,Rf_mkChar("deriv"));
			SET_STRING_ELT(names,1,Rf_mkChar("deriv_ln"));
		}

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}





extern "C" {

	SEXP dxdDiffusion7(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP re7, SEXP re8, SEXP re9, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

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


		/* calculate the derivatives */
		dxPDF7(t, resp, a, v, t0, w, sw, sv, st, eps, K, N, epsFLAG, Rda, Rdv, Rdt0, Rdw, Rdsw, Rdsv, Rdst, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,da);
		SET_VECTOR_ELT(out,1,dv);
		SET_VECTOR_ELT(out,2,dt0);
		SET_VECTOR_ELT(out,3,dw);
		SET_VECTOR_ELT(out,4,dsw);
		SET_VECTOR_ELT(out,5,dsv);
		SET_VECTOR_ELT(out,6,dst);


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

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}




extern "C" {

	SEXP dxpDiffusion7(SEXP re1, SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP re7, SEXP re8, SEXP re9, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP bo1) {

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


		/* calculate the derivatives */
		dxCDF7(t, resp, a, v, t0, w, sw, sv, st, eps, K, N, epsFLAG, Rda, Rdv, Rdt0, Rdw, Rdsw, Rdsv, Rdst, NThreads);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,da);
		SET_VECTOR_ELT(out,1,dv);
		SET_VECTOR_ELT(out,2,dt0);
		SET_VECTOR_ELT(out,3,dw);
		SET_VECTOR_ELT(out,4,dsw);
		SET_VECTOR_ELT(out,5,dsv);
		SET_VECTOR_ELT(out,6,dst);


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

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}



/* SAMPLING */

extern "C" {

	SEXP randWiener(SEXP re2, SEXP re3, SEXP re4, SEXP re5, SEXP re6, SEXP re7, SEXP re8, SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP in5, SEXP bo1) {

		/* define input variables */
		double a = REAL(re2)[0];
		double v = REAL(re3)[0];
		double w = REAL(re4)[0];
		double sv = REAL(re5)[0];
		double sw = REAL(re6)[0];
		double eps = REAL(re7)[0];
		double bound = REAL(re8)[0];

		int R = INTEGER(in1)[0];
		int K = INTEGER(in2)[0];
		int N = INTEGER(in3)[0];
		int NThreads = INTEGER(in4)[0];
		int choice = INTEGER(in5)[0];

		int epsFLAG = INTEGER(bo1)[0];


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP q = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP resp = PROTECT(Rf_allocVector(INTSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rq = REAL(q);
		int *Rresp = INTEGER(resp);

		/* sampling */
		run_make_rwiener(choice, N, a, v, w, sv, sw, R, bound, eps, K, epsFLAG, Rq, Rresp);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,q);
		SET_VECTOR_ELT(out,1,resp);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("q"));
		SET_STRING_ELT(names,1,Rf_mkChar("resp"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		return(out);
	}

}
