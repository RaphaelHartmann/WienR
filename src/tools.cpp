
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#include "tools.h"



double logsum(double xa, double xb) {
	double temp;
	if (xa <= -INFINITY) return xb;
	if (xb <= -INFINITY) return xa;
	if (xa > xb) temp = xa + log1p(exp(xb - xa));
	else temp = xb + log1p(exp(xa - xb));
	return temp;
}

double logdiff(double xa, double xb) {
	double result;
	if (xb <= -INFINITY) return(xa);
	if (xa <= -INFINITY) return(xb);
	if (xa > xb) result = (xa + log1p(-exp(xb - xa)));
	else
		if (xb > xa) result = (xb + log1p(-exp(xa - xb)));
		else result = -INFINITY;
	return result;
}

double rexp(double x) {
	double result;
	if (x <= 700.0)
		result = exp(x);
	else
		result = exp(700.00);
	return result;
}

double lognormal(double x) {
	return	-0.5*x*x - M_LN_SQRT_PI - 0.5*M_LN2;
}

/*
The original function lnnorm was written by
	Jean Marie Linhart
	StataCorp LP
	jlinhart@stata.com
	January 4, 2008
and later modified
*/
double lnnorm(double z)
{
        int lower ;

        double z2, y, s, p1, q1, p2, q2, t, a1, a2 ;
        double n, m ;

	if (z==0.0e0) return(log(0.50e0)) ;

	if (z > LNNORM_MAX_X) return(0.0e0);
        if (z <= LNNORM_MIN_X) return(-0.5e0*z*z);

        if (z<0.0e0) {
                z= -z ;
                lower=1 ;
        }
        else    lower=0 ;
        //upper = !lower ;

        z2 = z*z ;

        y = exp(-0.5*z2) / SQR2PI ;
        n = y/z ;

        if (!( ( z>2.00e0))) {
                z *= y ;
                s=z ;
                t=0.0e0 ;

                for (n=3.0e0;s!=t;n+=2.0e0) {
                        t=s ;
                        z *= ((z2)/ (n)) ;
                        s+=z ;
                }
                if (lower) return(log(0.50e0-s)) ;
                return(log(0.50e0+s)) ;
        }

        a1=2.0e0 ;
        a2=0.0e0 ;
        n=z2+3.0e0 ;
        p1=1.0e0 ;
        q1=z ;
        p2=(n-1.0e0) ;
        q2=n*z ;
        m = p1/q1 ;
        t = p2/q2 ;
        s = m ;

        for (n+=4.0e0; m!=t && s!=t; n+=4.0e0) {
                a1 -= 8.0 ;
                a2 += a1 ;
                s = a2*p1 + n*p2 ;
                p1=p2 ;
                p2=s ;
                s = a2*q1 + n*q2 ;
                q1=q2 ;
                q2=s ;
                s=m ;
                m=t ;
                if (q2>1.0e30) {
                        p1 /= 1.0e30 ;
                        p2 /= 1.0e30 ;
                        q1 /= 1.0e30 ;
                        q2 /= 1.0e30 ;
                }
 		            t = p2/q2;
        }
        t = lower ? log(t) - 0.5*z2 - log(SQR2PI) : log1p(-y*t);
        return(t) ;
}
/* ----------------------------------- */

double logMill(double x) {
	double m;
	if (x > 1.0e5) return -log(x);
	m = lnnorm(-x) - lognormal(x);
	return m;
}

/* rat_eval, small, intermediate, tail, gsl_cdf_ugaussian_Pinv
  are copied from the GNU scientific library version 2.6 */
double rat_eval(const double a[], const size_t na,
          const double b[], const size_t nb, const double x) {
  size_t i, j;
  double u, v, r;
  u = a[na - 1];
  for (i = na - 1; i > 0; i--)
    {
      u = x * u + a[i - 1];
    }
  v = b[nb - 1];
  for (j = nb - 1; j > 0; j--)
    {
      v = x * v + b[j - 1];
    }
  r = u / v;
  return r;
}

double small(double q) {
  const double a[8] = { 3.387132872796366608, 133.14166789178437745,
    1971.5909503065514427, 13731.693765509461125,
    45921.953931549871457, 67265.770927008700853,
    33430.575583588128105, 2509.0809287301226727
  };
  const double b[8] = { 1.0, 42.313330701600911252,
    687.1870074920579083, 5394.1960214247511077,
    21213.794301586595867, 39307.89580009271061,
    28729.085735721942674, 5226.495278852854561
  };
  double r = 0.180625 - q * q;
  double x = q * rat_eval (a, 8, b, 8, r);
  return x;
}

double intermediate(double r) {
  const double a[] = { 1.42343711074968357734, 4.6303378461565452959,
    5.7694972214606914055, 3.64784832476320460504,
    1.27045825245236838258, 0.24178072517745061177,
    0.0227238449892691845833, 7.7454501427834140764e-4
  };
  const double b[] = { 1.0, 2.05319162663775882187,
    1.6763848301838038494, 0.68976733498510000455,
    0.14810397642748007459, 0.0151986665636164571966,
    5.475938084995344946e-4, 1.05075007164441684324e-9
  };
  double x = rat_eval (a, 8, b, 8, (r - 1.6));
  return x;
}

double tail(double r) {
  const double a[] = { 6.6579046435011037772, 5.4637849111641143699,
    1.7848265399172913358, 0.29656057182850489123,
    0.026532189526576123093, 0.0012426609473880784386,
    2.71155556874348757815e-5, 2.01033439929228813265e-7
  };
  const double b[] = { 1.0, 0.59983220655588793769,
    0.13692988092273580531, 0.0148753612908506148525,
    7.868691311456132591e-4, 1.8463183175100546818e-5,
    1.4215117583164458887e-7, 2.04426310338993978564e-15
  };
  double x = rat_eval (a, 8, b, 8, (r - 5.0));
  return x;
}

double gsl_cdf_ugaussian_Pinv(const double P) {
  double r, x, pp;
  double dP = P - 0.5;
  if (P == 1.0)
    {
      return INFINITY;
    }
  else if (P == 0.0)
    {
      return -INFINITY;
    }

  if (fabs (dP) <= 0.425)
    {
      x = small (dP);

      return x;
    }

  pp = (P < 0.5) ? P : 1.0 - P;
  r = sqrt (-log (pp));
  if (r <= 5.0)
    {
      x = intermediate (r);
    }
  else
    {
      x = tail (r);
    }

  if (P < 0.5)
    {
      return -x;
    }
  else
    {
      return x;
    }
}
/* ----------------------------------- */
