
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#include "tools.h"
#include "pdf_fncs.h"
#include "cdf_fncs.h"
#include "fncs_seven.h"
#include "cubature.h"



/* DENSITY */

/* dependencies */

/* integrand density */
int int_ddiff(unsigned dim, const double *x, void *p, unsigned fdim, double *retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	// double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t - tau), a, nu, omega, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = -0.5 * pow(y, 2) - M_LN_SQRT_PI - 0.5 * M_LN2 + log1p(temp) - 2 * log1p(-temp);

		double integrand = exp(ldW + temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/da */
int int_daddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, nu, omega, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dadwiener(low_or_up * (t-tau), a, nu, omega, ldW, val_ptr, errorW, K, epsFLAG);
		double dda = val_ptr[0];

		double integrand = dda * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dv */
int int_dvddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, nu, omega, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dvdwiener(low_or_up * (t-tau), a, nu, omega, ldW, val_ptr);

		double integrand = val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dt0 */
int int_dt0ddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, nu, omega, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		dtdwiener(t-tau, a, -low_or_up*nu, wn, ldW, val_ptr, errorW, K, epsFLAG);

		double integrand = -val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dz */
int int_dwddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, nu, omega, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dwdwiener(low_or_up * (t-tau), a, nu, omega, ldW, val_ptr, errorW, K, epsFLAG);

		double integrand = val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dsz */
int int_dswddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? w + sw * (x[1] - 0.5) : w + sw * (x[0] - 0.5);
	double tau = sv ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0);
	double temp_sw = sv ? (x[1]-0.5) : (x[0]-0.5);

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, nu, omega, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dwdwiener(low_or_up * (t-tau), a, nu, omega, ldW, val_ptr, errorW, K, epsFLAG);

		double integrand = temp_sw * val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dsv */
int int_dsvddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = pow(x[0], 2);
	double y = x[0] / (1 - temp);
	double nu = v + sv * y;
	double omega = sw ? w + sw * (x[1] - 0.5) : w;
	double tau = sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0);

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, nu, omega, errorW, K, epsFLAG);

		double temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dvdwiener(low_or_up * (t-tau), a, nu, omega, ldW, val_ptr);

		double integrand = y * val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dst0 */
int int_dst0ddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );
	double temp_st0 = sw ? (sv ? -x[2] : -x[1]) : (sv ? -x[1] : -x[0]);

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, nu, omega, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		dtdwiener(t-tau, a, -low_or_up*nu, wn, ldW, val_ptr, errorW, K, epsFLAG);

		double integrand = temp_st0 * val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dt */
int int_dtddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, nu, omega, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		dtdwiener(t-tau, a, -low_or_up*nu, wn, ldW, val_ptr, errorW, K, epsFLAG);

		double integrand = val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* calculate density for 7-param diffusion */
void ddiff(int choice, double t, int low_or_up, double a, double v, double t0, double w, double sw, double sv, double st, double myerr, int K, int epsFLAG, int Neval, double *derivF, double *Rerr) {

	//double result;
	//double error;

	double value;
	// double valueln;

	double *val_ptr = &value;
	double errorW = myerr ? myerr*0.1 : 1.e-12*0.1;

	my_params params = {t, low_or_up, a, v, t0, w, sw, sv, st, errorW, K, epsFLAG, val_ptr};

	int dim = (sw!=0)+(sv!=0)+(st!=0);

	double *xmin = (double*)malloc(dim * sizeof(double));
	double *xmax = (double*)malloc(dim * sizeof(double));

	// 0  = s (v); 1 = u (w), 2 = v (w)
	if(sv) {
		xmin[0] = -1; xmax[0] = 1;
		for (int i = 1; i < dim; i++) {
			xmin[i] = 0;
			xmax[i] = 1;
		}
	} else {
		for (int i = 0; i < dim; i++) {
			xmin[i] = 0;
			xmax[i] = 1;
		}
	}
	if (st) xmax[dim-1] = fmin(1.0, (t-t0)/st);

	double reltol = 0.0;
	double abstol = myerr ? myerr*0.9 : 1.e-12*0.9;

	double val, err;

	int Meval = Neval;

	//	printf("%u-dim integral, tolerance = %g\n", dim, tol);
	switch (choice) {
		case 0: hcubature(1, int_ddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
	  case 1: hcubature(1, int_daddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 2: hcubature(1, int_dvddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 3: hcubature(1, int_dt0ddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 4: hcubature(1, int_dwddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 5: hcubature(1, int_dswddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 6: hcubature(1, int_dsvddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 7: hcubature(1, int_dst0ddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 8: hcubature(1, int_dtddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
	}
	//if(err > abstol) Rprintf("absolute error not achieved: %g < %g\n", abstol, err);

	free(xmin); free(xmax);
	*derivF = val;
	if (*Rerr < err+errorW) *Rerr = err+errorW;

}

/*-----------------------------------------------*/

/* DISTRIBUTION */

/* dependencies */

/* integrand distribution */
int int_pdiff(unsigned dim, const double *x, void *p, unsigned fdim, double *retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	// double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		double lpW = pwiener((t-tau), a, -low_or_up*nu, wn, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		double integrand = exp(lpW + temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/da */
int int_dapdiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		double lpW = pwiener((t-tau), a, -low_or_up*nu, wn, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dapwiener(low_or_up, (t-tau), a, nu, omega, lpW, val_ptr, errorW, K, epsFLAG);

		double integrand = val_ptr[0] * exp(temp2);
//if(std::isnan(integrand)) Rprintf("t-tau, a, v, w, lp:%8g%8g%8g%8g%8g\n", t-tau, a, nu, omega, lpW);
		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dv */
int int_dvpdiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		double lpW = pwiener((t-tau), a, -low_or_up*nu, wn, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dvpwiener(low_or_up, (t-tau), a, nu, omega, lpW, val_ptr, errorW, K, epsFLAG);

		double integrand = val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dt0 */
int int_dt0pdiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	// double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, nu, omega, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		double integrand = -exp(ldW + temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dz */
int int_dwpdiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		double lpW = pwiener((t-tau), a, -low_or_up*nu, wn, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dwpwiener(low_or_up, (t-tau), a, nu, omega, lpW, val_ptr, errorW, K, epsFLAG);

		double integrand = val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dsz */
int int_dswpdiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );
	double temp_sw = sv ? (x[1]-0.5) : (x[0]-0.5);

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		double lpW = pwiener((t-tau), a, -low_or_up*nu, wn, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dwpwiener(low_or_up, (t-tau), a, nu, omega, lpW, val_ptr, errorW, K, epsFLAG);

		double integrand = temp_sw * val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dsv */
int int_dsvpdiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );


	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		double lpW = pwiener((t-tau), a, -low_or_up*nu, wn, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dvpwiener(low_or_up, (t-tau), a, nu, omega, lpW, val_ptr, errorW, K, epsFLAG);

		double integrand = y * val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dst0 */
int int_dst0pdiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
	my_params *params = static_cast<my_params*>(p);
	double t = (params->t);
	int low_or_up = (params->low_or_up);
	double a = (params->a);
	double v = (params->v);
	double t0 = (params->t0);
	double w = (params->w);
	double sw = (params->sw);
	double sv = (params->sv);
	double st = (params->st);
	double errorW = (params->errorW);
	int K = (params->K);
	int epsFLAG = (params->epsFLAG);
	// double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );
	double temp_st0 = sw ? (sv ? -x[2] : -x[1]) : (sv ? -x[1] : -x[0]);

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, nu, omega, errorW, K, epsFLAG);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		double integrand = temp_st0 * exp(ldW + temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* calculate distribution for 7-param diffusion */
void pdiff(int choice, double t, int low_or_up, double a, double v, double t0, double w, double sw, double sv, double st, double myerr, int K, int epsFLAG, int Neval, double *derivF, double *Rerr) {

	//ouble result;
	//double error;

	double value;
	// double valueln;

	double *val_ptr = &value;
	double errorW = myerr ? myerr*0.1 : 1.e-12*0.1;

	my_params params = {t, low_or_up, a, v, t0, w, sw, sv, st, errorW, K, epsFLAG, val_ptr};

	int dim = (sw!=0)+(sv!=0)+(st!=0);

	double *xmin = (double*)malloc(dim * sizeof(double));
	double *xmax = (double*)malloc(dim * sizeof(double));

	// 0  = s (v); 1 = u (w), 2 = v (w)
	if(sv) {
		xmin[0] = -1; xmax[0] = 1;
		for (int i = 1; i < dim; i++) {
			xmin[i] = 0;
			xmax[i] = 1;
		}
	} else {
		for (int i = 0; i < dim; i++) {
			xmin[i] = 0;
			xmax[i] = 1;
		}
	}
	if (st) xmax[dim-1] = fmin(1.0, (t-t0)/st);

	double reltol = 0.0;
	double abstol = myerr ? myerr*0.9 : 1.e-12*0.9;

	double val, err;

	int Meval = Neval;

	//	printf("%u-dim integral, tolerance = %g\n", dim, tol);
	switch (choice) {
		case 0: hcubature(1, int_pdiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
	  case 1: hcubature(1, int_dapdiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 2: hcubature(1, int_dvpdiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 3: hcubature(1, int_dt0pdiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 4: hcubature(1, int_dwpdiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 5: hcubature(1, int_dswpdiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 6: hcubature(1, int_dsvpdiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
		case 7: hcubature(1, int_dst0pdiff, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err); break;
	}
	//if(err > abstol) Rprintf("absolute error not achieved: %g < %g\n", abstol, err);

	free(xmin); free(xmax);
	*derivF = val;
	if (*Rerr < err+errorW) *Rerr = err+errorW;

}
