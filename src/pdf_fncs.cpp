
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#include "cstdio"
#include "pdf_fncs.h"
#include "tools.h"
#include <Rinternals.h>



/* DENSITY */

/* calculate number of terms needed for short t */
double ks(double t, double w, double eps) {
	double K1 = (sqrt(2.0 * t) + w) / 2.0;
	double u_eps = fmin(-1.0, M_LN2 + M_LNPI + 2.0 * log(t) + 2.0 * (eps));
	double	arg = -t * (u_eps - sqrt(-2.0 * u_eps - 2.0));
	double 	K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
	return ceil(fmax(K1, K2));
}

/* calculate number of terms needed for large t */
double kl(double q, double v, double w, double err) {
	double K1 = 1.0 / (M_PI * sqrt(q)), K2=0.0;
	double temp = -2.0 * (log(M_PI * q) + err);
	if (temp>=0) K2 = sqrt(temp/(pow(M_PI, 2) * q));
	return ceil(fmax(K1,K2));
}

/* calculate terms of the sum for short t */
double logfs(double t, double w, int K) {
	if (w == 0) return -INFINITY;
	double	fplus = -INFINITY, fminus = -INFINITY, twot = 2.0 * t;
	if (K > 0)
		for (int k = K; k >= 1; k--) {
			double temp1 = w + 2.0 * k, temp2 = w - 2.0 * k;

			fplus = logsum(log(temp1) - pow(temp1, 2) / twot, fplus);
			fminus = logsum(log(-temp2) - pow(temp2, 2) / twot, fminus);
		}
	fplus = logsum(log(w) - pow(w, 2) / twot, fplus);
	return  -0.5 * M_LN2 - M_LN_SQRT_PI - 1.5 * log(t) + logdiff(fplus, fminus);
}

/* calculate terms of the sum for large t */
double logfl(double q, double v, double w, int K) {
	if (w == 0) return -INFINITY;
	double fplus = -INFINITY, fminus = -INFINITY;
	double halfq = q / 2.0;
	for (int k = K; k >= 1; k--) {
		double temp = k * M_PI;
		double check = sin(temp * w);
		if (check > 0) fplus = logsum(log(k) - pow(temp, 2) * halfq + log(check), fplus);
		else fminus = logsum(log(k) - pow(temp, 2) * halfq + log(-check), fminus);
	}
	return	logdiff(fplus, fminus) + M_LNPI;
}

/* calculate density */
double dwiener(double q, double a, double vn, double wn, double err, int K, int epsFLAG) {
	if (q == 0.0) {
		return -INFINITY;
	}
	double kll, kss, ans, v, w;
	if(!epsFLAG && K==0) {
		err = -27.63102;  // exp(err) = 1.e-12
		epsFLAG = 1;
	}
	else if(!epsFLAG && K>0) err = -27.63102;  // exp(err) = 1.e-12
	else if(epsFLAG) err = log(err);

	if (q >= 0) {
		w = 1.0 - wn;
		v = -vn;
	}
	else {
		q = fabs(q);
		w = wn;
		v = vn;
	}

	double q_asq = q / pow(a, 2);
	ans = 0.0;

	/* calculate the number of terms needed for short t*/
	double lg1 = (-v * a * w - (pow(v, 2)) * q / 2.0) - 2.0*log(a);
	double es = (err - lg1);
	kss = ks(q_asq, w, es);
	/* calculate the number of terms needed for large t*/
	double el = es;
	kll = kl(q_asq, v, w, el);

	// if small t is better
	if (2 * kss <= kll) {
		if((epsFLAG && kss<K) || !epsFLAG) kss = K;
		ans = lg1 + logfs(q_asq, w, static_cast<int>(kss));
	}
	// if large t is better
	else {
		if((epsFLAG && kll<K) || !epsFLAG) kll = K;
		ans = lg1 + logfl(q_asq, v, w, static_cast<int>(kll));
	}

	return ans;
}

/*-----------------------------------------------*/



/* d/dt DENSITY */

/* calculate number of terms needed for short t */
double dtaks(double t, double w, double eps) {
	double K1 = (sqrt(3.0*t) + w) / 2.0;
	double u_eps = fmin(-1.0, (log(8.0 / 27.0) + M_LNPI + 4.0*log(t) + 2.0*eps)/3.0);
	double	arg = -3.0 * t * (u_eps - sqrt(-2.0 * u_eps - 2.0));
	double 	K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;

	return ceil(fmax(K1, K2));
}

/* calculate number of terms needed for large t */
double dtakl(double q, double v, double a, double err) {
	double K1 = sqrt(3.0 / q) / M_PI;
	double u_eps = fmin(-1.0, err + log(0.6) + M_LNPI + 2.0 * log(q) );
	double	arg = -2.0/M_PISQ/q*(u_eps - sqrt(-2.0 * u_eps - 2.0));
	double 	kl = (arg > 0) ? sqrt(arg) : K1;

	return ceil(fmax(kl,K1));
}

/* calculate terms of the sum for short t */
void logdtfs(double t, double w, int K, double &erg, int &newsign) {
  double	fplus = -INFINITY, fminus = -INFINITY, twot=2.0*t;
  {
    if (K > 0)
    for (int k = K; k >= 1; k--) {
      double temp1 = w + 2.0 * k, temp2 = w - 2.0 * k;
      fplus = logsum(3.0*log(temp1) - temp1 * temp1 / twot, fplus);
      fminus = logsum(3.0*log(-temp2) - temp2 * temp2 / twot, fminus);
    }
  }
  newsign = 1;
  fplus = logsum(3.0*log(w) - w * w / twot, fplus);
  erg = logdiff(fplus, fminus);
  if (fplus < fminus) newsign = -1;
}

/* calculate terms of the sum for large t */
void logdtfl(double q, double w, int K, double &erg, int &newsign) {
  double fplus = -INFINITY, fminus = -INFINITY;
  double  halfq = q / 2.0;
  for (int k = K; k >= 1; k--) {
    double temp = M_PI * k, zwi = sin(temp * w);
    if (zwi > 0) {
      fplus = logsum(3.0*log(k) - temp * temp * halfq + log(zwi), fplus);

    }
    if (zwi < 0) {
      fminus = logsum(3.0*log(k) - temp * temp * halfq + log(-zwi), fminus);

    }
  }
  erg = logdiff(fplus, fminus); newsign = (fplus > fminus) ? 1 : -1;
}

/* calculate derivative of density with respect to t */
void dtdwiener(double q, double a, double v, double w, double ld, double *derivF, double *derivlnF, double err, int K, int epsFLAG) {
	if (q == 0.0) {
		*derivF = 0.0;
		*derivlnF = 0.0/0.0;
	} else {
		double kll, kss, ans;
		if(!epsFLAG && K==0) {
			err = -27.63102; // exp(err) = 1.e-12
			epsFLAG = 1;
		}
		else if(!epsFLAG && K>0) err = -27.63102; // exp(err) = 1.e-12
		else if(epsFLAG) err = log(err);

		/* prepare some variables */
	  double q_asq = q / pow(a, 2);
	  double ans0 = -pow(v, 2) / 2.0;
	  double la = 2.0*log(a);
	  double lg1 = -v * a * w - pow(v, 2) * q / 2.0 - la;
	  double factor = lg1 - la;


		/* calculate the number of terms needed for short t */
	  double es = err - lg1 + ld;
	  es = es + la;
		kss = dtaks(q_asq, w, es);
		/* calculate the number of terms needed for large t */
	  double el = (err - lg1 + ld) ;
	  el = el +la;
	  kll = dtakl(q_asq, v, a, el);

	  // if small t is better
	  if (2*kss < kll) {
			if((epsFLAG && kss<K) || !epsFLAG) kss = K;
	    double erg; int newsign;
	    logdtfs(q_asq, w, static_cast<int>(kss), erg, newsign);
	    ans = ans0 - 1.5 / q + newsign * exp(factor - 1.5 * M_LN2 - M_LN_SQRT_PI - 3.5 * log(q_asq) + erg - ld);
	  }
	  // if large t is better...
	  else {
			if((epsFLAG && kll<K) || !epsFLAG) kll = K;
	    double erg; int newsign;
	    logdtfl(q_asq, w, static_cast<int>(kll), erg, newsign);
	    ans = ans0 - newsign * exp(factor + 3.0*M_LNPI - M_LN2 + erg - ld);
	  }
	  //return ans;
		*derivlnF = ans;
		*derivF = ans*exp(ld);
	}
}

/*-----------------------------------------------*/



/* d/da DENSITY */

/* calculate derivative of density with respect to a */
void dadwiener(double q, double a, double vn, double wn, double ld, double *derivF, double *derivlnF, double err, int K, int epsFLAG) {
	if (q == 0.0) {
		*derivF = 0.0;
		*derivlnF = 0.0/0.0;
	} else {
		double kll, kss, ans, v, w;
		if(!epsFLAG && K==0) {
			err = -27.63102; // exp(err) = 1.e-12
			epsFLAG = 1;
		}
		else if(!epsFLAG && K>0) err = -27.63102; // exp(err) = 1.e-12
		else if(epsFLAG) err = log(err);

		if (q >= 0) {
			w = 1.0 - wn;
			v = -vn;
		}
		else {
			q = fabs(q);
			w = wn;
			v = vn;
		}

		double la = log(a), lq = log(q);

		/* prepare some variables */
		double q_asq = q / pow(a, 2);
		double ans0 =  - v * w;
		double lg1 = -v * a*w - pow(v, 2)*q / 2.0 - 2.0*la;
		double factor = lg1 - 3.0*la;

		/* calculate the number of terms needed for short t */
		double es = err - lg1 + ld;
		es = es + la;
		es = es -  M_LN2 + 2.0*la - lq;
		kss = dtaks(q_asq, w, es);
		/* calculate the number of terms needed for large t */
		double el = err - lg1 + ld;
		el = el + la;
		el = el - M_LN2 + 2.0*la - lq;
		kll = dtakl(q_asq, v, a, el);

		// if small t is better
		if (2 * kss < kll) {
			if((epsFLAG && kss<K) || !epsFLAG) kss = K;
			double erg; int newsign;
			logdtfs(q_asq, w, static_cast<int>(kss), erg, newsign);
			ans = ans0 + 1.0 / a - newsign * exp(-0.5*M_LN2 - M_LN_SQRT_PI - 2.5*lq + 4.0*la + lg1 + erg - ld);
		}
		// if large t is better...
		else {
			if((epsFLAG && kll<K) || !epsFLAG) kll = K;
			double erg; int newsign;
			logdtfl(q_asq, w, static_cast<int>(kll), erg, newsign);
			ans = ans0 - 2.0 / a + newsign * exp(lq + factor + 3.0*M_LNPI + erg - ld);
		}
		//return ans;
		*derivlnF = ans;
		*derivF = ans*exp(ld);
	}
}
/*-----------------------------------------------*/



/* d/dv DENSITY */

/* calculate derivative of density with respect to v */
void dvdwiener(double q, double a, double vn, double wn, double ld, double *derivF, double *derivlnF) {
	if (q == 0.0) {
		*derivF = 0.0;
		*derivlnF = 0.0/0.0;
	} else {
		double ans, v, w;
		int sign = 1;

		if (q >= 0) {
			w = 1.0 - wn;
			v = -vn;
			sign = -1;
		}
		else {
			q = fabs(q);
			w = wn;
			v = vn;
		}

		ans =  sign*(- a * w - v * q);

		*derivlnF = ans;
		*derivF = ans*exp(ld);
	}
}
/*-----------------------------------------------*/



/* d/dw DENSITY */

/* calculate number of terms needed for short t */
double dwks(double t, double w, double eps) {
	double K1 = (sqrt(3.0*t) + w) / 2.0;
	double u_eps = fmin(-1.0, 2.0*(eps) + M_LN2 + M_LNPI + 2.0*log(t));
	double	arg = -t * (u_eps - sqrt(-2.0 * u_eps - 2.0));
	double 	K2 = (arg > 0) ? (sqrt(arg) + w) / 2.0 : K1;
	return ceil(fmax(K1, K2));
}

/* calculate number of terms needed for large t */
double dwkl(double q, double v, double err) {
	double K1 = sqrt(2.0 / q) / M_PI;
	double u_eps = fmin(-1.0, log(4.0 / 9.0) + 2.0*M_LNPI + 3.0*log(q) + 2.0*(err));
	double	arg = -(u_eps - sqrt(-2.0 * u_eps - 2.0));
	double 	K2 = (arg > 0) ? 1.0 / M_PI * sqrt(arg / q) : K1;
	return ceil(fmax(K1, K2));
}

/* calculate terms of the sum for short t */
void logdwfs(double t, double w, int K, double &erg, int &sign) {
	double	fplus = -INFINITY, fminus = -INFINITY, twot = 2 * t;
	for (int k = K; k >= 1; k--) {
		double temp1 = pow((w + 2 * k), 2), temp2 = pow((w - 2 * k), 2), temp3 = temp1 - t, temp4 = temp2 - t;
		if (temp3 > 0) fplus = logsum(log(temp3) - temp1 / twot, fplus);
		else if (temp3 < 0) fminus = logsum(log(-(temp3)) - temp1 / twot, fminus);
		if (temp4 > 0) fplus = logsum(log(temp4) - temp2 / twot, fplus);
		else if (temp4 < 0) fminus = logsum(log(-(temp4)) - temp2 / twot, fminus);
	}
	double temp = pow(w, 2), temp1 = temp - t;
	if (temp1 > 0) fplus = logsum(log(temp1) - temp / twot, fplus);
	else if (temp1 < 0) fminus = logsum(log(-(temp1)) - temp / twot, fminus);
	erg = logdiff(fplus, fminus);
	sign = (fplus < fminus)? -1: 1;
	// 1/sqrt(2*pi*t^5) added outside
}

/* calculate terms of the sum for large t */
void logdwfl(double q, double v,  double w, int K, double &erg, int &sign) {
	double fplus=-INFINITY, fminus = -INFINITY;
	double  halfq=q/2.0;

	for (int k = K; k >= 1; k--) {
		double temp = M_PI * k;
		double x = cos(temp * w);
		if (x > 0)
			fplus = logsum(2.0*log(k) - pow(temp, 2) * halfq + log(x), fplus);
		else if (x < 0)
			fminus = logsum(2.0*log(k) - pow(temp, 2) * halfq + log(-x), fminus);
	}
	erg = logdiff(fplus, fminus);
	sign = (fplus < fminus)? -1: 1;
// pi^2 added outside
}

/* calculate derivative of density with respect to w */
void dwdwiener(double q, double a, double vn, double wn, double ld, double *derivF, double *derivlnF, double err, int K, int epsFLAG) {
	if (q == 0.0) {
		*derivF = 0.0;
		*derivlnF = 0.0/0.0;
	} else {
		double kll, kss, ans, v, w;
		if(!epsFLAG && K==0) {
			err = -27.63102; // exp(err) = 1.e-12
			epsFLAG = 1;
		}
		else if(!epsFLAG && K>0) err = -27.63102; // exp(err) = 1.e-12
		else if(epsFLAG) err = log(err);

		int sign = 1;

		if (q>=0) {
			w = 1.0 - wn;
			v = -vn;
			sign = -1;
		}
		else {
			q = fabs(q);
			w = wn;
			v = vn;
		}

		/* prepare some variables */
		double q_asq = q / pow(a, 2);
		double ans0 = -v * a;
		double lg1 = (-v * a*w - pow(v, 2)*(q) / 2.0) - 2.0*log(a);
		double ls = -lg1 + ld;
		double ll = -lg1 + ld;

		/* calculate the number of terms needed for short t and large t */
		kss = dwks(q_asq, w, err+ls);
		kll = dwkl(q_asq, v, err+ll);

		// if small t is better
		if (2 * kss < kll) {
			if((epsFLAG && kss<K) || !epsFLAG) kss = K;
			double erg; int newsign;
			logdwfs(q_asq, w, static_cast<int>(kss), erg, newsign);
			ans = sign*(ans0 - newsign * exp(erg - ls - 2.5 * log(q_asq) - 0.5 * M_LN2 - 0.5 * M_LNPI));
		}
		// if large t is better...
		else {
			if((epsFLAG && kll<K) || !epsFLAG) kll = K;
	  	double erg; int newsign;
			logdwfl(q_asq, v, w, static_cast<int>(kll), erg, newsign);
			ans = sign*(ans0 + newsign * exp(erg - ll + 2.0 * M_LNPI));
		}
		//return ans*sign;
		*derivlnF = ans;
		*derivF = ans*exp(ld);
	}
}
/*-----------------------------------------------*/



/* d/dw DENSITY */

/* calculate number of terms needed for short t */
void dxks(double q, double t, double w, double a, double eps, double &Kas, double &Kws) {
	double la = log(a);
	double lq = log(t);
	double es_a = eps + la;
	es_a = es_a - M_LN2 +  2*la - lq;

	double K1_a = (sqrt(3.0*q) + w) / 2.0;
	double K1_w = K1_a;

	double u_eps_a = fmin(-1.0, (log(8.0 / 27.0) + M_LNPI + 4.0*log(q) + 2.0*es_a)/3.0);
	double arg_a = -3.0 * q * (u_eps_a - sqrt(-2.0 * u_eps_a - 2.0));
	double u_eps_w = fmin(-1, 2*(eps) +M_LN2+ M_LNPI + 2.0*log(q));
	double arg_w = -q * (u_eps_w - sqrt(-2 * u_eps_w - 2));

	double K2_a = (arg_a > 0) ? 0.5 * (sqrt(arg_a) - w) : K1_a;
	double K2_w = (arg_w > 0) ? 0.5 * (sqrt(arg_w) + w) : K1_w;

	Kas = ceil(fmax(K1_a, K2_a));
	Kws = ceil(fmax(K1_w, K2_w));
}

/* calculate number of terms needed for large t */
void dxkl(double q, double t, double v, double a, double err, double &Kal, double &Kwl) {
	double la = log(a);
	double lq = log(t);
	double el_a = err + la;
	el_a = el_a - M_LN2 + 2*la - lq;

	double K1_a = sqrt(3 / q)/M_PI;
	double K1_w = sqrt(2 / q) / M_PI;

	double u_eps_a = fmin(-1.0, el_a + log(0.6) + M_LNPI + 2.0 * log(q) );
	double arg_a = -2.0/M_PISQ/q*(u_eps_a - sqrt(-2.0 * u_eps_a - 2.0));
	double u_eps_w = fmin(-1, log(4.0 / 9.0) + 2 * M_LNPI + 3.0*log(q) + 2.0*(err));
	double arg_w = -(u_eps_w - sqrt(-2 * u_eps_w - 2));

	double K2_a = (arg_a > 0) ? sqrt(arg_a) : K1_a;
	double K2_w = (arg_w > 0) ? 1.0 / M_PI * sqrt(arg_w / q) : K1_w;

	Kal = ceil(fmax(K1_a, K2_a));
	Kwl = ceil(fmax(K1_w, K2_w));
}

/* calculate terms of the sum for short t */
void logdxfs(double q, double w, int Kas, int Kws, double &erg_a, double &erg_w, int &sign_a, int &sign_w) {
	double fplus_a, fminus_a, fplus_w, fminus_w;
	fplus_a = fminus_a = fplus_w = fminus_w = -INFINITY;
	double twot=2.0*q;

	int K = fmax(Kas, Kws);
  {
    if (K > 0) {
			for (int k = K; k >= 1; k--) {
	      double temp1 = w + 2.0 * k, temp2 = w - 2.0 * k;
				if (k <= Kas) {
					fplus_a = logsum(3.0*log(temp1) - temp1 * temp1 / twot, fplus_a);
		      fminus_a = logsum(3.0*log(-temp2) - temp2 * temp2 / twot, fminus_a);
				}
				if (k <= Kws) {
					double temp1sq = pow(temp1, 2), temp2sq = pow(temp2, 2);
					double temp3 = temp1sq - q, temp4 = temp2sq - q;
					if (temp3 > 0) fplus_w = logsum(log(temp3) - temp1sq / twot, fplus_w);
					else if (temp3 < 0) fminus_w = logsum(log(-(temp3)) - temp1sq / twot, fminus_w);
					if (temp4 > 0) fplus_w = logsum(log(temp4) - temp2sq / twot, fplus_w);
					else if (temp4 < 0) fminus_w = logsum(log(-(temp4)) - temp2sq / twot, fminus_w);
				}
	    }
		}
		double temp = pow(w, 2), temp1 = temp - q;
		if (temp1 > 0) fplus_w = logsum(log(temp1) - temp / twot, fplus_w);
		else if (temp1 < 0) fminus_w = logsum(log(-(temp1)) - temp / twot, fminus_w);
  }
	fplus_a = logsum(3.0*log(w) - w * w / twot, fplus_a);

	erg_a = logdiff(fplus_a, fminus_a); sign_a = (fplus_a > fminus_a) ? 1 : -1;
	erg_w = logdiff(fplus_w, fminus_w); sign_w = (fplus_w > fminus_w) ? 1 : -1;
}

/* calculate terms of the sum for large t */
void logdxfl(double q, double w, int Kal, int Kwl, double &erg_a, double &erg_w, int &sign_a, int &sign_w) {
	double fplus_a, fminus_a, fplus_w, fminus_w;
	fplus_a = fminus_a = fplus_w = fminus_w = -INFINITY;
	double halfq = q / 2.0;

	int K = fmax(Kal, Kwl);
	for (int k = K; k >= 1; k--) {
		double temp = M_PI * k, zwi = sin(temp * w), x = cos(temp * w);
		if (zwi > 0 && k <= Kal) {
			if (k <= Kal) fplus_a = logsum(3.0*log(k) - temp * temp * halfq + log(zwi), fplus_a);
		}
		if (zwi < 0 && k <= Kal) {
			if (k <= Kal) fminus_a = logsum(3.0*log(k) - temp * temp * halfq + log(-zwi), fminus_a);
		}
		if (x > 0 && k <= Kwl) {
			if (k <= Kwl) fplus_w = logsum(2.0 * log(k) - pow(temp, 2) * halfq + log(x), fplus_w);
		}
		if (x < 0 && k <= Kwl) {
			if (k <= Kwl) fminus_w = logsum(2.0 * log(k) - pow(temp, 2) * halfq + log(-x), fminus_w);
		}
	}
	erg_a = logdiff(fplus_a, fminus_a); sign_a = (fplus_a > fminus_a) ? 1 : -1;
	erg_w = logdiff(fplus_w, fminus_w); sign_w = (fplus_w > fminus_w) ? 1 : -1;
}

/* calculate derivative of density with respect to w */
void dxdwiener(double q, double a, double vn, double wn, double ld, double err, int K, int epsFLAG, double *da, double *da_ln, double *dv, double *dv_ln, double *dw, double *dw_ln) {

	if (q == 0.0) {
		// -----------d/da------------
		*da_ln = 0.0/0.0;
		*da = 0.0;

		// -----------d/dv------------
		*dv_ln = 0.0/0.0;
		*dv = 0.0;

		// -----------d/dw------------
		*dw_ln = 0.0/0.0;
		*dw = 0.0;
	} else {
		double v, w, kal, kas, kwl, kws;
		int sign = 1;
		if(!epsFLAG && K==0) {
			err = -27.63102; // exp(err) = 1.e-12
			epsFLAG = 1;
		}
		else if(!epsFLAG && K>0) err = -27.63102; // exp(err) = 1.e-12
		else if(epsFLAG) err = log(err);

		if (q >= 0) {
			w = 1.0 - wn;
			v = -vn;
			sign = -1;
		}
		else {
			q = fabs(q);
			w = wn;
			v = vn;
		}

		double q_asq = q / pow(a, 2);

		double ans_a, ans0_a =  - v * w;
		double ans_w, ans0_w = -v * a;
		double ans_v = sign*(- a * w - v * q);

		double la = log(a);
		double lq = log(q);
		double lg1 = -v * a * w - pow(v, 2) * q / 2.0 - 2.0*la;
		double es = err - lg1 + ld;
		double el = es;

		double factor_a = lg1 - 3.0*la;
		double factor_w = -lg1 + ld;

		dxks(q_asq, q, w, a, es, kas, kws);
		dxkl(q_asq, q, v, a, el, kal, kwl);

		double S = kas + kws;
		double L = kal + kwl;

		if (2*S < L) {
			if((epsFLAG && S < 2*K) || !epsFLAG) kas = kws = K;
			double erg_a, erg_w; int newsign_a, newsign_w;
			logdxfs(q_asq, w, static_cast<int>(kas), static_cast<int>(kws), erg_a, erg_w, newsign_a, newsign_w);
			ans_a = ans0_a + 1.0 / a - newsign_a * exp(-0.5*M_LN2 - M_LN_SQRT_PI - 2.5*lq + 4 * la + lg1 + erg_a - ld);
			ans_w = ans0_w - newsign_w * exp(erg_w - factor_w - 2.5 * log(q_asq) - 0.5 * M_LN2 - 0.5 * M_LNPI);
		}
		// if large t is better...
		else {
			if((epsFLAG && L < 2*K) || !epsFLAG) kal = kwl = K;
			double erg_a, erg_w; int newsign_a, newsign_w;
			logdxfl(q_asq, w, static_cast<int>(kal), static_cast<int>(kwl), erg_a, erg_w, newsign_a, newsign_w);
			ans_a = ans0_a - 2.0 / a + newsign_a * exp(lq + factor_a + 3 * M_LNPI + erg_a - ld);
			ans_w = ans0_w + newsign_w * exp(erg_w - factor_w + 2.0 * M_LNPI);
		}

		// -----------d/da------------
		*da_ln = ans_a;
		*da = ans_a*exp(ld);

		// -----------d/dv------------
		*dv_ln = ans_v;
		*dv = ans_v*exp(ld);

		// -----------d/dw------------
		*dw_ln = ans_w*sign;
		*dw = (ans_w*sign)*exp(ld);
	}

}
/*-----------------------------------------------*/
