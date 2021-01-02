
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#include "cstdio"
#include "cdf_fncs.h"
#include "tools.h"
#include <Rinternals.h>
#include <limits>


/* DISTRIBUTION */

/* P term in the distribution function */
double logP(int pm, double a, double v, double w) {
	if (pm == 1) { v = -v; w = 1.0 - w; }
	if (fabs(v) == 0.0) return log1p(-w);

	double prob;
	double e = (-2.0 * v * a * (1.0 - w));
	if (e < 0) {
		double tt = exp(e);
		tt = log1p(-tt) - logdiff(2 * v * a * w, e);
		prob = tt;
	}
	else {
		double tt = log1p(-exp(-e)) - log1p(-exp(2 * v * a));
		prob = tt;
	}
	return prob;
}

/* calculate number of terms needed for short t */
double Ks(double t, double v, double a, double w, double eps)
{
	double K1 = 0.5 * (fabs(v) / a * t - w);
	double arg = fmax(0, fmin(1, exp(v*a*w + pow(v, 2)*t / 2 + (eps)) / 2));
	double K2 = (arg==0) ? INFINITY : (arg==1) ? -INFINITY : -sqrt(t) / 2 / a * gsl_cdf_ugaussian_Pinv(arg);
	return ceil(fmax(K1, K1 + K2));
}

/* calculate number of terms needed for large t */
double Kl(double t, double v, double a, double w, double err) {
	double api = a / M_PI, vsq = pow(v, 2);
	double sqrtL1 = sqrt(1 / t) * api;
	double sqrtL2 = sqrt(fmax(1.0, -2 / t * pow(api, 2) * (err+log(M_PI*t / 2 * (vsq + pow((M_PI / a), 2))) + v * a * w + vsq * t / 2)));
	return ceil(fmax(sqrtL1, sqrtL2));
}

/* calculate terms of the sum for short t */
double logFs(double t, double v, double a, double w, int K)
{
	double fplus = -INFINITY, fminus = -INFINITY;
	double sqt = sqrt(t), temp = -v * a*w - pow(v, 2)*t / 2;
	double vt = v * t;

	for (int k = K; k >= 0; k--)
	{
		double rj = a*(2 * k + w);
		double dj = lognormal(rj / sqt);
		double pos1 = dj + logMill((rj - vt) / sqt);
		double pos2 = dj + logMill((rj + vt) / sqt);
		fplus = logsum(logsum(pos1, pos2), fplus);
		rj = a*(2.0 * k + 2.0 - w);
		dj =  lognormal(rj / sqt);
		double neg1 = dj + logMill((rj - vt) / sqt);
		double neg2 = dj + logMill((rj + vt) / sqt);
		fminus = logsum(logsum(neg1, neg2), fminus);
	}

	return logdiff(fplus, fminus)+temp;
}

/* calculate terms of the sum for large t */
double logFl(double q, double v, double a, double w, int K)
{
	double fplus = -INFINITY, fminus = -INFINITY;
	double la = log(a), lv = log(fabs(v));
	double F = -INFINITY;
	for (int k = K; k >= 1; k--) {
		double temp0 = log(k * 1.0), temp1 = k * M_PI, temp2 = temp1 * w;
		double check = sin(temp2);
		if (check > 0) {
			double temp = temp0 - logsum(2 * lv, 2 * (temp0 + M_LNPI - la)) - 0.5 * pow((temp1 / a), 2) * q + log(check);
			fplus = logsum(temp, fplus);
		}
		else if (check < 0)
		{
			double temp = temp0 - logsum(2 * lv, 2 * (temp0 + M_LNPI - la)) - 0.5 * pow((temp1 / a), 2) * q + log(-check);
			fminus = logsum(temp, fminus);
		}
	}
	F = logdiff(fplus, fminus);
	return (F - v * a * w - 0.5 * pow(v, 2) * q);
}

/* calculate distribution */
double pwiener(double q, double a, double v, double w, double err, int K, int epsFLAG) {
	//if (q == 0) return(GSL_NEGINF);
	double Kll, Kss, ans;
	if(!epsFLAG && K==0) {
		err = -27.63102; // exp(err) = 1.e-12
		epsFLAG = 1;
	}
	else if(!epsFLAG && K>0) err = -27.63102; // exp(err) = 1.e-12
	else if(epsFLAG) err = log(err);

	if (std::isinf(q)) return logP(0, a, v, w);

	Kss = Ks(q, v, a, w, err);
	Kll = Kl(q, v, a, w, err);
	double lg = M_LN2 + M_LNPI - 2.0 * log(a);

  if (3 * Kss < Kll) {
		if((epsFLAG && Kss<K) || !epsFLAG) Kss = K;
		ans = logFs(q, v, a, w, static_cast<int>(Kss));
		//Rprintf("short hier\n");
	}
	else {
		if((epsFLAG && Kll<K) || !epsFLAG) Kll = K;
		ans = logdiff(logP(0, a, v, w), lg + logFl(q, v, a, w, static_cast<int>(Kll)));
		//Rprintf("large\n");
	}
	return ans;
}

/*-----------------------------------------------*/



/* d/da DISTRIBUTION */

/* P term in the derivative with respect to d/da and d/dv */
double davlogP(int pm, double a, double v, double w) {

	double tt;
	if (pm == 1) {
		w = 1.0 - w;
		v = -v;
	}

	if (fabs(v)==0.0) return(-w);

	if (v < 0) {
		double emw = (2.0 * v * a * (1.0 - w)), ew = 2 * a * v * w, e = 2 * a * v;
		tt = M_LN2 + emw - log1p(-exp(emw));
		double temp = log1p(-exp(ew)) - log1p(-exp(e));
		if (log(w) > temp) {
			tt += logdiff(log(w), temp);
			tt = rexp(tt);
		}
		else {
			tt += logdiff(temp, log(w));
			tt = -rexp(tt);
		}


	}
	else {
		double emw = (-2.0 * v * a * (1.0 - w)), e = (-2 * a * v);
		tt = M_LN2 - log1p(-exp(emw));
		double temp = logdiff(emw, e) - log1p(-exp(e));
		if (log(w) > temp) {
			tt += logdiff(log(w), temp);
			tt = -rexp(tt);
		}
		else {
			tt += logdiff(temp, log(w));
			tt = rexp(tt);
		}
	}
	if (std::isfinite(tt)) return(tt);
	else {
		Rprintf("dalogprob %20g%20g%20g\n", a, v, w);
		//std::cout << "dalogprob " << setw(20) << a << setw(20) << v << setw(20) << w << std::endl;
		return(-INFINITY);
	}
	return tt;
}

/* extension for d/da */
double dalogP(int pm, double a, double v, double w, double dav) {
	if (fabs(v) == 0.0) return 0.0;
	double tt;
	tt = dav * v;
	tt = (pm == 1) ? -tt : tt;
	if (std::isfinite(tt)) return(tt);
	else {
		Rprintf("dalogprob %20g%20g%20g\n", a, v, w);
		//std::cout << "dalogprob " << setw(20) << a << setw(20) << v << setw(20) << w << std::endl;
		return(-INFINITY);
	}
	return tt;
}

/* calculate number of terms needed for short t */
void dakS(double q, double a, double v, double w, double err, double &Kas) {
	double la = log(a), sqt = sqrt(q), lv = log1p(pow(v, 2) * q);
	double factor = v * a * w + pow(v, 2) * q / 2 + err;
	double wdash = fmin(w, 1.0 - w);

  double K1a = sqt / a - wdash;

	double ueps = fmin(-1, 2 * (factor + la - lv) + M_LNPI);
	Kas = (sqt * sqrt(-(ueps - sqrt(-2 * ueps - 2))) - a * wdash) / a;
	Kas = ceil(fmax(fmax(Kas, K1a), 1.0));
}

/* calculate number of terms needed for large t */
void dakL(double q, double a, double v, double w, double err, double &Kal) {
	double lt = log(q), la = log(a);
	double factor = v * a * w + pow(v, 2) * q / 2 + err;
  double K1 = a / M_PI / sqrt(q);
  double C1 = M_LN2 - logsum(2 * log(fabs(v)), 2 * (M_LNPI - la)); C1 = logsum(C1, lt);
	double alphka = fmin(factor + M_LNPI + lt + la - M_LN2 - C1, 0.0);
	Kal = ceil(fmax(fmax(sqrt(-2 * alphka / q) * a / M_PI, K1),1.0));
}

/* calculate terms of the sum for short t */
void logdaFs(int pm, int Ksa, double t, double a, double v, double w, double &derF, double lp) {
	if (pm == 1) {
		v = -v;
		w = 1 - w;
	}

	double Fj = 0.0;
	derF = 0.0;

	double sqt = sqrt(t), vt = v * t;

	for (int k = Ksa; k >= 0; k--)
	{
		double ta1, ta2, ta3, ta4;
		double temp3;

		double rj = 2 * k * a + a * w;
		double dj = lognormal(rj / sqt);

		double x = rj - vt, xsqt = x / sqt;
		double temp = rexp(dj + logMill(xsqt)), temp2 = exp(dj);
		temp3 = temp * (-vt) - sqt * temp2;
		ta1 =  temp3 * (2 * k + w);

		x = rj + vt, xsqt = x / sqt;
		temp = rexp(dj + logMill(xsqt));
		temp3 = temp * vt - sqt * temp2;
		ta2 = temp3 * (2 * k + w);

		rj = (2 * k + 1) * a + a * (1 - w);
		dj = lognormal(rj / sqt);
		x = rj - vt, xsqt = x / sqt;
		temp = rexp(dj + logMill(xsqt)), temp2 = exp(dj);
		temp3 = temp * (-vt) - sqt * temp2;
		ta3 = -temp3 * (2 * k + 2.0 - w);

		x = rj + vt, xsqt = x / sqt;
		temp = rexp(dj + logMill(xsqt));
		temp3 = temp * vt - sqt * temp2;
		ta4 = -temp3 * (2 * k + 2.0 - w);

		derF = (ta1 + ta2 + ta3 + ta4) + derF;

	}

//	Fj = rexp(lp + v * a * w + 0.5 * pow(v, 2) * t);
//	derF = -v * w + derF / t / Fj;
	Fj = rexp(v * a * w + 0.5 * pow(v, 2) * t);
	derF = -v * w * exp(lp) + derF / t/ Fj;
/*	if(Fj == 0) {Rprintf("One derivative d/da had to be calculated numerically\n");
		double epsilon = .00000000001;
		double tempL = logdiff(logsum(lp, log(epsilon)), logdiff(lp, log(epsilon)))-log(2*epsilon);
		derF = exp(tempL-lp);
	}*/
}

/* calculate terms of the sum for large t */
void logdaFl(int pm,  int Kal, double q, double a, double v, double w, double &derF, double lp) {
	if (pm == 1) {
		v = -v;
		w = 1 - w;
	}

	derF = 0.0;
	for (int k = Kal; k >= 1; k--) {
		double sin1 = sin(M_PI * k * w), kpi = k * M_PI, kpia2 = pow((kpi / a), 2), ekpia2q = exp(-0.5 * kpia2 * q), denom = 1.0 / (pow(v, 2) + kpia2), denomk = k * denom;
		double last = pow(kpi, 2) / pow(a, 3) * (q + 2.0 * denom) * denomk * ekpia2q;
		derF = derF - last * sin1;
	}

	double evaw = exp(-v * a * w - 0.5 * pow(v, 2) * q);
	double temp = rexp(logP(0, a, v, w));
	double dav = davlogP(0, a, v, w);
	double pia2 = 2 * M_PI / pow(a, 2);

	lp = exp(lp);
	double tempa = dalogP(0, a, v, w, dav) * temp;
	derF = (-2 / a - v * w) * (lp-temp) + derF * pia2 * evaw;
	derF = tempa + derF;
	//derF /= p;

}

/* calculate derivative of distribution with respect to a */
void dapwiener(int pm, double q, double a, double v, double w, double lp, double *derivF, double err, int K, int epsFLAG) {
  if(!epsFLAG && K==0) {
		err = -27.63102; // exp(err) = 1.e-12
		epsFLAG = 1;
	}
	else if(!epsFLAG && K>0) err = -27.63102; // exp(err) = 1.e-12
	else if(epsFLAG) err = log(err);
  err = err;// + lp;

	double Kal, Kas;
	if (pm == -1) dakL(q, a, v, w, err, Kal); else dakL(q, a, -v, 1.0 - w, err, Kal);
	if (pm == -1) dakS(q, a, v, w, err, Kas); else dakS(q, a, -v, 1.0 - w, err, Kas);

	double erg;
	if (Kal < 4 * Kas) {
    if((epsFLAG && Kal<K) || !epsFLAG) Kal = K;
    logdaFl(pm, static_cast<int>(Kal), q, a, v, w, erg, lp);
  }
	else {
    if((epsFLAG && Kas<K) || !epsFLAG) Kas = K;
    logdaFs(pm, static_cast<int>(Kas), q, a, v, w, erg, lp);
  }
  *derivF = erg;
}
/*-----------------------------------------------*/



/* d/dv DISTRIBUTION */

/* extension for d/dv */
double dvlogP(int pm, double a, double v, double w, double dav) {
	int sign = 1; if (pm == 1) sign = -1;
	double tt = dav * a * sign;
	if (std::isfinite(tt)) return(tt); else
	{
		Rprintf("dvlogprob %20g%20g%20g\n", a, v, w);
		//std::cout << "dvlogprob " << setw(20) << a << setw(20) << v << setw(20) << w << std::endl;
		return(-INFINITY);
	}
	return tt;
}

/* calculate number of terms needed for short t */
void dvkS(double q, double a, double v, double w, double err, double &Kvs) {
	double lt = log(q), sqt = sqrt(q);
	double factor = v * a * w + pow(v, 2) * q / 2 + err;
	double wdash = fmin(w, 1.0 - w);

	double K1wv = fabs(v) / a * q - wdash;

	double alphKv = factor + 0.5 * (M_LN2 - lt + M_LNPI);
	Kvs = (alphKv < 0) ? sqt * sqrt(-2 * alphKv)/a -  wdash  : 0;
	Kvs = ceil(fmax(fmax(Kvs, K1wv), 1.0));
}

/* calculate number of terms needed for large t */
void dvkL(double q, double a, double v, double w, double err, double &Kvl) {
	double lt = log(q), la = log(a);
	double temp = -rexp(la - M_LNPI - 0.5*lt);

	double factor = v * a * w + pow(v, 2) * q / 2 + err;

	if (v == 0) Kvl = 1.0;
	else {
		double lv = log(fabs(v));
		double alphKv = rexp(factor + 0.5 * (7*M_LNPI + lt) - 2.5 * M_LN2 - 3 * la - lv);
		alphKv = fmax(0.0, fmin(1.0, alphKv));
		Kvl = fmax(ceil((alphKv == 0) ? INFINITY : (alphKv == 1) ? -INFINITY : temp * gsl_cdf_ugaussian_Pinv(alphKv)), 1.0);
	}
}

/* calculate terms of the sum for short t */
void logdvFs(int pm, int Ksv, double t, double a, double v, double w, double &derF, double lp)
{
	double sign = 1.0;
	if (pm == 1) {
		v = -v;
		w = 1 - w;
		sign = -1.0;
	}

	double Fj = 0.0;
	derF = 0.0;

	double sqt = sqrt(t);

	for (int k = Ksv; k >= 0; k--)
	{
		double tv1, tv2, tv3, tv4;

		double rj = 2 * k * a + a * w;
		double dj = lognormal(rj / sqt);

		double x = rj - v * t, xsqt = x / sqt;
		double temp = rexp(dj + logMill(xsqt));
		tv1 = -temp * x;

		x = rj + v * t, xsqt = x / sqt;
		temp = rexp(dj + logMill(xsqt));
		tv2 = temp * x;

		rj = (2 * k + 1) * a + a * (1 - w);
		dj = lognormal(rj / sqt);
		x = rj - v * t, xsqt = x / sqt;
		temp = rexp(dj + logMill(xsqt));
		tv3 = temp * x;

		x = rj + v * t, xsqt = x / sqt;
		temp = rexp(dj + logMill(xsqt));
		tv4 = -temp * x;

		derF = (tv1 + tv2 + tv3 + tv4) + derF;
	}

//	Fj = rexp(lp + v * a * w + 0.5 * pow(v, 2) * t);
//	derF = ((-w * a - v * t) + derF / Fj) * sign;
	Fj = rexp(v * a * w + 0.5 * pow(v, 2) * t);
	derF = ((-w * a - v * t) * exp(lp) + derF / Fj) * sign;
}

/* calculate terms of the sum for large t */
void logdvFl(int pm,  int Kvl, double q, double a, double v, double w, double &derF, double lp) {
	double sign = 1.0;
	if (pm == 1) {
		v = -v;
		w = 1 - w;
		sign = -1.0;
	}

	derF = 0.0;
	for (int k = Kvl; k >= 1; k--) {
		double sin1 = sin(M_PI * k * w), kpi = k * M_PI, kpia2 = pow((kpi / a), 2), ekpia2q = exp(-0.5 * kpia2 * q), denom = 1.0 / (pow(v, 2) + kpia2), denomk = k * denom;
		double last = denomk * denom * ekpia2q;
		derF = derF - last * sin1;
	}

	double evaw = exp(-v * a * w - 0.5 * pow(v, 2) * q);
	double temp = rexp(logP(0, a, v, w));
	double dav = davlogP(0, a, v, w);
	double pia2 = 2 * M_PI / pow(a, 2);

	lp = exp(lp);
	double tempv = dvlogP(0, a, v, w, dav) * temp;
	derF = (-w * a - v * q) * (lp-temp) + derF * (-2 * v) * pia2 * evaw;
	derF = (tempv + derF) * sign;
//	derF /= lp;
}

/* calculate derivative of distribution with respect to v */
void dvpwiener(int pm, double q, double a, double v, double w, double lp, double *derivF, double err, int K, int epsFLAG) {
	if(!epsFLAG && K==0) {
		err = -27.63102; // exp(err) = 1.e-12
		epsFLAG = 1;
	}
	else if(!epsFLAG && K>0) err = -27.63102; // exp(err) = 1.e-12
	else if(epsFLAG) err = log(err);
  err = err;// + lp;

	double Kvl, Kvs;
	if (pm == -1) dvkL(q, a, v, w, err, Kvl); else dvkL(q, a, -v, 1.0 - w, err, Kvl);
	if (pm == -1) dvkS(q, a, v, w, err, Kvs); else dvkS(q, a, -v, 1.0 - w, err, Kvs);

	double erg;
	if (Kvl < 4 * Kvs) {
		if((epsFLAG && Kvl<K) || !epsFLAG) Kvl = K;
		logdvFl(pm, static_cast<int>(Kvl), q, a, v, w, erg, lp);
	}
	else {
	 	if((epsFLAG && Kvs<K) || !epsFLAG) Kvs = K;
	 	logdvFs(pm, static_cast<int>(Kvs), q, a, v, w, erg, lp);
	}
	*derivF = erg;
}
/*-----------------------------------------------*/



/* d/dw DISTRIBUTION */

/* extension for d/dw */
double dwlogP(int pm, double a, double v, double w) {
	double tt = 1.0;
	if (pm == 1) {
		w = 1.0 - w;
		v = -v;
		tt = -1.0;
	}

	if (fabs(v) == 0.0) return -tt / (1.0 - w);
	if (v < 0) {
		double e = (2.0 * v * a * (1.0 - w));
		double temp = M_LN2 + e + log(fabs(v)) + log(a) - log1p(-exp(e));
		tt *= -exp(temp);
	}
	else
	{
		double e = -(2.0 * v * a * (1.0 - w));
		double temp = M_LN2 + log(fabs(v)) + log(a) - log1p(-exp(e));
		tt *= -exp(temp);
	}
	return tt;
}

/* calculate number of terms needed for short t */
void dwkS(double q, double a, double v, double w, double err, double &Kws) {
	double sqt = sqrt(q), lv = log1p(pow(v, 2) * q);
	double factor = v * a * w + pow(v, 2) * q / 2 + err;
	double wdash = fmin(w, 1.0 - w);

	double K1wv = fabs(v) / a * q - wdash;

	double alphKw = factor  - M_LN2 - lv;
	double arg = fmin(rexp(alphKw), 1.0);
	Kws = (arg == 0) ? INFINITY : (arg == 1) ? -INFINITY : -sqt / a * gsl_cdf_ugaussian_Pinv(arg) - wdash;

	Kws = ceil(fmax(fmax(Kws, K1wv), 1.0));
}

/* calculate number of terms needed for large t */
void dwkL(double q, double a, double v, double w, double err, double &Kwl) {
	double lt = log(q), la = log(a);
	double temp = -rexp(la - M_LNPI - 0.5*lt);
	double factor = v * a * w + pow(v, 2) * q / 2 + err;
	double  alphKw = rexp(factor + 0.5 * (M_LNPI + lt) - 1.5 * M_LN2 - la);
	Kwl = fmax(ceil((alphKw == 0) ? INFINITY : (alphKw == 1) ? -INFINITY : temp * gsl_cdf_ugaussian_Pinv(alphKw)), 1.0);
}

/* calculate terms of the sum for short t */
void logdwFs(int pm, int Ksw, double t, double a, double v, double w, double &derF, double lp)
{
	double sign = 1.0;
	if (pm == 1) {
		v = -v;
		w = 1 - w;
		sign = -1.0;
	}

	double Fj = 0.0;
	derF = 0.0;

	double sqt = sqrt(t), vt = v * t;
	for (int k = Ksw; k >= 0; k--)
	{
		double tw1, tw2, tw3, tw4;
		double temp3;

		double rj = 2 * k * a + a * w;
		double dj = lognormal(rj / sqt);

		double x = rj - v * t, xsqt = x / sqt;
		double temp = rexp(dj + logMill(xsqt)), temp2 = exp(dj);
		temp3 = temp * (-vt) - sqt * temp2;
		tw1 = temp3 * a;

		x = rj + v * t, xsqt = x / sqt;
		temp = rexp(dj + logMill(xsqt));
		temp3 = temp * vt - sqt * temp2;
		tw2 = temp3 * a;

		rj = (2 * k + 1) * a + a * (1 - w);
		dj = lognormal(rj / sqt);
		x = rj - v * t, xsqt = x / sqt;
		temp = rexp(dj + logMill(xsqt)), temp2 = exp(dj);
		temp3 = temp * (-vt) - sqt * temp2;
		tw3 = -temp3 * (-a);

		x = rj + v * t, xsqt = x / sqt;
		temp = rexp(dj + logMill(xsqt));
		temp3 = temp * vt - sqt * temp2;
		tw4 = -temp3 * (-a);

		derF = (tw1 + tw2 + tw3 + tw4) + derF;
	}

//	Fj = rexp(lp + v * a * w + 0.5 * pow(v, 2) * t);
//	derF = (-v * a + derF / t / Fj) * sign;
	Fj = rexp(v * a * w + 0.5 * pow(v, 2) * t);
	derF = (-v * a * exp(lp) + derF / t / Fj) * sign;
}

/* calculate terms of the sum for large t */
void logdwFl(int pm,  int Kwl, double q, double a, double v, double w, double &derF, double lp) {
	double sign = 1.0;
	if (pm == 1) {
		v = -v;
		w = 1 - w;
		sign = -1.0;
	}

	derF = 0.0;
	for (int k = Kwl; k >= 1; k--) {
		double kpi = k * M_PI, kpia2 = pow((kpi / a), 2), ekpia2q = exp(-0.5 * kpia2 * q), denom = 1.0 / (pow(v, 2) + kpia2), denomk = k * denom;
		double last = k * M_PI * denomk * ekpia2q;
		derF = derF - last * cos(kpi * w);
	}

	double evaw = exp(-v * a * w - 0.5 * pow(v, 2) * q);
	double temp = rexp(logP(0, a, v, w));
	double tempw = dwlogP(0, a, v, w) * temp;
	double pia2 = 2 * M_PI / pow(a, 2);

	lp = exp(lp);
	derF = -v * a * (lp-temp) + derF * pia2 * evaw;
	derF = (tempw + derF) * sign;
//	derF /= lp;
}

/* calculate derivative of distribution with respect to w */
void dwpwiener(int pm, double q, double a, double v, double w, double lp, double *derivF, double err, int K, int epsFLAG) {
	if(!epsFLAG && K==0) {
		err = -27.63102; // exp(err) = 1.e-12
		epsFLAG = 1;
	}
	else if(!epsFLAG && K>0) err = -27.63102; // exp(err) = 1.e-12
	else if(epsFLAG) err = log(err);
  err = err;// + lp;

	double Kwl, Kws;
	if (pm == -1) dwkL(q, a, v, w, err, Kwl); else dwkL(q, a, -v, 1.0 - w, err, Kwl);
	if (pm == -1) dwkS(q, a, v, w, err, Kws); else dwkS(q, a, -v, 1.0 - w, err, Kws);

	double erg;
	if (Kwl < 4 * Kws) {
		if((epsFLAG && Kwl<K) || !epsFLAG) Kwl = K;
		logdwFl(pm, static_cast<int>(Kwl), q, a, v, w, erg, lp);
	}
	else {
	 	if((epsFLAG && Kws<K) || !epsFLAG) Kws = K;
	 	logdwFs(pm, static_cast<int>(Kws), q, a, v, w, erg, lp);
	}
	*derivF = erg;
}
/*-----------------------------------------------*/



/* d/dt d/da d/dv d/dw DISTRIBUTION */

/* calculate number of terms needed for short t */
void dxkS(double q, double a, double v, double w, double err, double &Kas, double &Kvs, double &Kws) {
	double lt = log(q), la = log(a), sqt = sqrt(q), lv = log1p(pow(v, 2) * q);
	double factor = v * a * w + pow(v, 2) * q / 2 + err;
	double wdash = fmin(w, 1.0 - w);

	double K1wv = fabs(v) / a * q - wdash;
	double K1a = sqt / a - wdash;


	double ueps = fmin(-1, 2 * (factor + la - lv) + M_LNPI);
	Kas = (sqt * sqrt(-(ueps - sqrt(-2 * ueps - 2))) - a * wdash) / a;

	double alphKw = factor  - M_LN2 - lv;
	double arg = fmin(rexp(alphKw), 1.0);
	Kws = (arg == 0) ? INFINITY : (arg == 1) ? -INFINITY : -sqt / a * gsl_cdf_ugaussian_Pinv(arg) - wdash;

	double alphKv = factor + 0.5 * (M_LN2 - lt + M_LNPI);
	Kvs = (alphKv < 0) ? sqt * sqrt(-2 * alphKv)/a -  wdash  : 0;

	Kws = ceil(fmax(fmax(Kws, K1wv), 1.0));
	Kas = ceil(fmax(fmax(Kas, K1a), 1.0));
	Kvs = ceil(fmax(fmax(Kvs, K1wv), 1.0));
}

/* calculate number of terms needed for large t */
void dxkL(double q, double a, double v, double w, double err, double &Kal, double &Kvl, double &Kwl) {
	double lt = log(q), la = log(a);
	double temp = -rexp(la - M_LNPI - 0.5*lt);

	double factor = v * a * w + pow(v, 2) * q / 2 + err;

	double  alphKw = rexp(factor + 0.5 * (M_LNPI + lt) - 1.5 * M_LN2 - la);
	Kwl = fmax(ceil((alphKw == 0) ? INFINITY : (alphKw == 1) ? -INFINITY : temp * gsl_cdf_ugaussian_Pinv(alphKw)), 1.0);

	alphKw = fmax(0.0, fmin(1.0, alphKw));
	if (v == 0) Kvl = 1.0; else {
		double lv = log(fabs(v));
		double alphKv = rexp(factor + 0.5 * (7*M_LNPI + lt) - 2.5 * M_LN2 - 3 * la - lv);
		alphKv = fmax(0.0, fmin(1.0, alphKv));
		Kvl = fmax(ceil((alphKv == 0) ? INFINITY : (alphKv == 1) ? -INFINITY : temp * gsl_cdf_ugaussian_Pinv(alphKv)), 1.0);
	}

	double K1 = a / M_PI / sqrt(q);
	double C1 = M_LN2 - logsum(2 * log(fabs(v)), 2 * (M_LNPI - la)); C1 = logsum(C1, lt);
	double alphka = fmin(factor + M_LNPI + lt + la - M_LN2 - C1, 0.0);
	Kal = ceil(fmax(fmax(sqrt(-2 * alphka / q) * a / M_PI, K1),1.0));
}

/* calculate terms of the sum for short t */
void logdxFs(int pm, int Ksa, int Ksv, int Ksw, double t, double a, double v, double w, double lp, double &Fa, double &Fv, double &Fw) {
	double sign = 1.0;
	if (pm == 1) {
		v = -v;
		w = 1 - w;
		sign = -1.0;
	}

	double Fj = 0.0;
	Fa = Fv = Fw = 0.0;

	double sqt = sqrt(t), vt = v * t;

	int K = int(fmax(Ksa, Ksv)); K = int(fmax(K, Ksw));

	int Kaw = int(fmax(Ksa, Ksw));

	for (int k = K; k >= 0; k--)
	{
		double ta1, ta2, ta3, ta4;
		double tv1, tv2, tv3, tv4;
		double tw1, tw2, tw3, tw4;
		double temp3;

		double rj = 2 * k * a + a * w;
		double dj = lognormal(rj / sqt);

		double x = rj - v * t, xsqt = x / sqt;
		double temp = rexp(dj + logMill(xsqt)), temp2 = exp(dj);
		if (k <= Ksv) tv1 = -temp * x;
		if (k <= Kaw) temp3 = temp * (-vt) - sqt * temp2;
		if (k <= Ksa) ta1 =  temp3 * (2 * k + w);
		if (k <= Ksw) tw1 = temp3 * a;

		x = rj + v * t, xsqt = x / sqt;
		temp = rexp(dj + logMill(xsqt));
		if (k <= Ksv) tv2 = temp * x;
		if (k <= Kaw) temp3 = temp * vt - sqt * temp2;
		if (k <= Ksa) ta2 = temp3 * (2 * k + w);
		if (k <= Ksw) tw2 = temp3 * a;

		rj = (2 * k + 1) * a + a * (1 - w);
		dj = lognormal(rj / sqt);
		x = rj - v * t, xsqt = x / sqt;
		temp = rexp(dj + logMill(xsqt)), temp2 = exp(dj);
		if (k <= Ksv) tv3 = temp * x;
		if (k <= Kaw) temp3 = temp * (-vt) - sqt * temp2;
		if (k <= Ksa) ta3 = -temp3 * (2 * k + 2.0 - w);
		if (k <= Ksw) tw3 = -temp3 * (-a);

		x = rj + v * t, xsqt = x / sqt;
		temp = rexp(dj + logMill(xsqt));
		if (k <= Ksv) tv4 = -temp * x;
		if (k <= Kaw) temp3 = temp * vt - sqt * temp2;
		if (k <= Ksa) ta4 = -temp3 * (2 * k + 2.0 - w);
		if (k <= Ksw) tw4 = -temp3 * (-a);

		if (k <= Ksv) Fv = (tv1 + tv2 + tv3 + tv4) + Fv;
		if (k <= Ksa) Fa = (ta1 + ta2 + ta3 + ta4) + Fa;
		if (k <= Ksw) Fw = (tw1 + tw2 + tw3 + tw4) + Fw;
	}

//	Fj = rexp(lp + v * a * w + 0.5 * pow(v, 2) * t);
//	Fv = ((-w * a - v * t) + Fv / Fj) * sign;
//	Fa = -v * w + Fa / t / Fj;
//	Fw = (-v * a + Fw / t / Fj) * sign;
	Fj = rexp(v * a * w + 0.5 * pow(v, 2) * t);
	Fv = ((-w * a - v * t) * exp(lp) + Fv / Fj) * sign;
	Fa = -v * w * exp(lp) + Fa / t / Fj;
	Fw = (-v * a * exp(lp) + Fw / t / Fj) * sign;
}

/* calculate terms of the sum for large t */
void logdxFl(int pm, int Kal, int Kvl, int Kwl, double q, double a, double v, double w, double lp, double &Fa, double &Fv, double &Fw) {
	double sign = 1.0;
	if (pm == 1) {
		v = -v;
		w = 1 - w;
		sign = -1.0;
	}

	Fa = Fv = Fw = 0.0;
	int K = int(fmax(Kal, Kvl)); K = int(fmax(K, Kwl));

	for (int k = K; k >= 1; k--) {
		double sin1 = sin(M_PI * k * w), kpi = k * M_PI, kpia2 = pow((kpi / a), 2), ekpia2q = exp(-0.5 * kpia2 * q),
			denom = 1.0 / (pow(v, 2) + kpia2), denomk = k * denom;
		if (k <= Kwl) {
			double last = k * M_PI * denomk * ekpia2q;
			Fw = Fw - last * cos(kpi * w);
		}
		if (k <= Kvl) {
			double last = denomk * denom * ekpia2q;
			Fv = Fv - last * sin1;
		}
		if (k <= Kal) {
			double last = pow(kpi, 2) / pow(a, 3) * (q + 2.0 * denom) * denomk * ekpia2q;
			Fa = Fa - last * sin1;
		}
	}

	double evaw = exp(-v * a * w - 0.5 * pow(v, 2) * q);
	double temp = rexp(logP(0, a, v, w));
	double tempw = dwlogP(0, a, v, w) * temp;
	double dav = davlogP(0, a, v, w);
	double pia2 = 2 * M_PI / pow(a, 2);

	lp = exp(lp);
	Fw = -v * a * (lp-temp) + Fw * pia2 * evaw;
	Fw = (tempw + Fw) * sign;

	double tempv = dvlogP(0, a, v, w, dav) * temp;
	Fv = (-w * a - v * q) * (lp-temp) + Fv * (-2 * v) * pia2 * evaw;
	Fv = (tempv + Fv) * sign;

	double tempa = dalogP(0, a, v, w, dav) * temp;
	Fa = (-2 / a - v * w) * (lp-temp) + Fa * pia2 * evaw;
	Fa = tempa + Fa;

//	Fa /= lp;
//	Fv /= lp;
//	Fw /= lp;
}

/* calculate derivative of distribution with respect to a */
void dxpwiener(int pm, double q, double a, double v, double w, double lp, double err, int K, int epsFLAG, double *da, double *dv, double *dw) {
	if(!epsFLAG && K==0) {
		err = -27.63102; // exp(err) = 1.e-12
		epsFLAG = 1;
	}
	else if(!epsFLAG && K>0) err = -27.63102; // exp(err) = 1.e-12
	else if(epsFLAG) err = log(err);
	err = err;// + lp;

	double Kal, Kvl, Kwl;
	double Kas, Kvs, Kws;
	if (pm == 0) dxkL(q, a, v, w, err, Kal, Kvl, Kwl); else dxkL(q, a, -v, 1.0 - w, err, Kal, Kvl, Kwl);
	if (pm == 0) dxkS(q, a, v, w, err, Kas, Kvs, Kws); else dxkS(q, a, -v, 1.0 - w, err, Kas, Kvs, Kws);

	double S =  Kas + Kvs + Kws;
	double L =  Kal + Kvl + Kwl;

	double Fa = 0.0, Fv = 0.0, Fw = 0.0;
	if (L < 4*S) {
		if((epsFLAG && L < 3*K) || !epsFLAG) Kal = Kvl = Kwl = K;
		logdxFl(pm, static_cast<int>(Kal), static_cast<int>(Kvl), static_cast<int>(Kwl), q, a, v, w, lp, Fa, Fv, Fw);
	}
	else {
		if((epsFLAG && S < 3*K) || !epsFLAG) Kas = Kvs = Kws = K;
		logdxFs(pm, static_cast<int>(Kas), static_cast<int>(Kvs), static_cast<int>(Kws), q, a, v, w, lp, Fa, Fv, Fw);
	}
	// -----------d/da------------
	*da = Fa;

	// -----------d/dv------------
	*dv = Fv;

	// -----------d/dw------------
	*dw = Fw;

}
/*-----------------------------------------------*/
