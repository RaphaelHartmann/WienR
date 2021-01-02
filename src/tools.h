
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#ifndef TOOLS_H
#define TOOLS_H
#include "cmath"

/*------------------------------------------------*/
/*MATHEMATICAL CONSTANTS*/

#ifndef M_PISQ
#define M_PISQ        9.86960440108935861883449099987615113531
#endif

#ifndef LNNORM_MAX_X
#define LNNORM_MAX_X  38.0
#endif

#ifndef LNNORM_MIN_X
#define LNNORM_MIN_X  -1.00e9
#endif

#ifndef SQR2PI
#define SQR2PI        2.506628274631000502415765e0
#endif
  /*-------------------------*/

  /* The following mathematical constants
  are copied from GSL version 2.6 */
#ifndef M_LNPI
#define M_LNPI        1.14472988584940017414342735135      /* ln(pi) */
#endif
  /*-------------------------*/

  /* The following mathematical constants
  are copied from Rmath.h */
#ifndef M_PI
#define M_PI          3.141592653589793238462643383280        /* pi */
#endif

#ifndef M_LN2
#define M_LN2         0.693147180559945309417232121458        /* ln(2) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI  0.572364942924700087071713675677        /* log(sqrt(pi)) == log(pi)/2 */
#endif
  /*-------------------------*/
/*------------------------------------------------*/



/* USEFUL FUNCTIONS */
double lnnorm(double);
double logsum(double, double);
double logdiff(double, double);
double rexp(double);
double lognormal(double);
double logMill(double);

/* rat_eval, small, intermediate, tail, gsl_cdf_ugaussian_Pinv
  are copied from the GNU scientific library version 2.6 */
double rat_eval(const double, const size_t, const double, const size_t, const double);
double small(double);
double intermediate(double);
double tail(double);
double gsl_cdf_ugaussian_Pinv(const double);
/* -------------------------------------------------------- */

struct my_params {
	double t;
	int low_or_up;
	double a;
	double v;
	double t0;
	double w;
	double sw;
	double sv;
	double st;
	double errorW;
	int K;
	int epsFLAG;
	double *val_ptr;
};



#endif
