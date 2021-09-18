

#include "tools.h"
#include "cdf_fncs.h"
#include "pdf_fncs.h"
#include "rwiener.h"
#include "cubature.h"
#include <cmath>
#include <algorithm>    // std::sort

#define t_tilde 2.5


double fun_upper(int k, double x, std::vector<piece> upper) {
	int i = 1;
	//	int k = static_cast<int>(upper.size());
	while ((i != k) && (x >= upper[i].z)) i++;
	i = i - 1;
	double t = upper[i].absc + upper[i].slope * (x - upper[i].center);
	return t;
}

void generate_intervals(int& k, double totallow, std::vector<point> h, std::vector<piece>& lower, std::vector<piece>& upper, std::vector<double>& s) {
	k = static_cast<int>(h.size());

	lower.clear(); upper.clear(); piece low, up;
	for (int j = 0; j != k; j++) {
		double z;
		if (j == 0) z = totallow;
		else z = (h[j].h - h[j - 1].h - h[j].x * h[j].dh + h[j - 1].x * h[j - 1].dh) / (h[j - 1].dh - h[j].dh);

		up.z = z;
		up.slope = h[j].dh; up.absc = h[j].h; up.center = h[j].x;
		upper.push_back(up);
		if (j == 0) low.z = totallow;
		else low.z = h[j - 1].x;
		lower.push_back(low);
	}
	low.z = h[k - 1].x; lower.push_back(low);

	double sum = -INFINITY, t; s.clear();
	for (int i = 0; i != k; i++) {
		if (i == 0) t = fun_upper(k, upper[i + 1].z, upper);
		else
			if (i < k - 1) {
				double sl = upper[i].slope;
				t = upper[i].absc - upper[i].center * sl + logdiff(upper[i + 1].z * sl,
					upper[i].z * sl);
			}
			else {
				t = (fun_upper(k, upper[i].z, upper));
			}
		t -= std::log(fabs(upper[i].slope));

		sum = logsum(sum, t);
		s.push_back(sum);
	}
}

bool update_intervals(int k, double totallow, point new_point, std::vector<point>& h, std::vector<piece>& lower, std::vector<piece>& upper, std::vector<double>& s) {
	double x = new_point.x;
	bool flag = false;
	int i = 0;
	k = static_cast<int>(h.size());
	while ((i != k) && (x > h[i].x))  i++;

	h.insert(h.begin() + i, new_point);
	piece low;
	int j = i + 1;
	low.z = h[i].x;
	lower.insert(lower.begin() + j, low);
	j = i;
	piece up;
	double z;
	if (j == 0) z = totallow;
	else z = (h[j].h - h[j - 1].h - h[j].x * h[j].dh + h[j - 1].x * h[j - 1].dh) / (h[j - 1].dh - h[j].dh);
	up.z = z;

	up.slope = h[j].dh; up.absc = h[j].h; up.center = h[j].x;
	if (i < k) upper[i] = up; else upper.push_back(up);


	if (i < k) {
		j = i + 1;
		z = (h[j].h - h[j - 1].h - h[j].x * h[j].dh + h[j - 1].x * h[j - 1].dh) / (h[j - 1].dh - h[j].dh);
		up.z = z;
		up.slope = h[j].dh; up.absc = h[j].h; up.center = h[j].x;
		upper.insert(upper.begin() + j, up);
	}

	k = k + 1;

	double sum = 0, t;
	std::vector<double> sold = s;
	//	s.clear();
//	if (i + 1 == k) std::cout << "!!!!!!!!!!!!!!!!!";
	{
		if (i > 1) sum = sold[i - 2];
		int iup = (i + 1 == k) ? i + 1 : i + 2;
		int ilow = (i == 0) ? 0 : i - 1;
		for (int j = ilow; j != iup; j++) {
			if (j == 0) t = fun_upper(k, upper[j + 1].z, upper);
			else
				if (j < k - 1) {
					double sl = upper[j].slope;
					t = upper[j].absc - upper[j].center * sl + logdiff(upper[j + 1].z * sl,
						upper[j].z * sl);
				}
				else {
					t = (fun_upper(k, upper[j].z, upper));
				}
			t -= std::log(fabs(upper[j].slope));
			if (j == 0) sum = t;
			else sum = logsum(sum, t);
			if (j != i) s[j] = sum; else s.insert(s.begin() + j, sum);
		}
	}
	if (i + 1 < k) {
		if (sold[i] < s[i + 1]) {
			Rprintf("Denkfehler %20g%20g\n", sold[i], s[i + 1]);
			flag = true;
		}
		double temp = logdiff(sold[i], s[i + 1]);
		for (int j = i + 2; j < k; j++) {
			s[j] = logdiff(sold[j - 1], temp);
		}
	}
	if (flag) {
		generate_intervals(k, totallow, h, lower, upper, s);
		flag = false;
	}
	return flag;
}

double fun_lower(int k, double x, std::vector<point> h, std::vector<piece> lower) {
	int i = 1;
	//	int k = static_cast<int>(lower.size());
	k = k + 1;
	while ((i != k) && (x >= lower[i].z)) i++;
	i = i - 1; double t;
	if ((i == 0) || (i == k - 1)) t = -INFINITY;
	else t = ((h[i].x - x) * h[i - 1].h + (x - h[i - 1].x) * h[i].h) / (h[i].x - h[i - 1].x);

	return t;
}

double inverse_distribution(int k, double xstar, std::vector<piece> upper, std::vector<double> s, double bound, bool& flag) {
	double sum = 0, t;

	if (bound == INFINITY) sum = s[k - 1];
	else {
		if (bound <= upper[k - 1].z) {
			Rprintf("Problem in inverse\n");
			flag = true;
		}
		double sl = upper[k - 1].slope;
		t = upper[k - 1].absc - upper[k - 1].center * sl + logdiff(bound * sl,
			upper[k - 1].z * sl);
		t -= std::log(fabs(sl));
		s[k - 1] = logsum(t, s[k - 2]);
		sum = s[k - 1];
	}
	int j = 0;
	double temp = std::log(xstar) + sum;
	while (temp > s[j]) j++;
	if (j > k - 1) { Rprintf("Wie das?\n"); }


	double sl = upper[j].slope;
	double help = std::log(fabs(sl)); int sign = sl > 0 ? 1 : -1;
	if (std::isnan(sl)) {
		flag = true;
		Rprintf("slope is infinity\n");
	}

	if (j > 0) temp = logdiff(temp, s[j - 1]);
	help = help + temp - upper[j].absc + upper[j].center * sl;
	if (sign == 1) temp = logsum(help, upper[j].z * sl);
	else temp = logdiff(upper[j].z * sl, help);
	t = temp / sl;

	if (t < upper[j].z) {
		Rprintf("\nnanu j=%d; k-1=%d; t=%g; upper[j]=%g; upper[j+1]=%g; s[j-1]=%g; upper slope=%g; upper absc=%g; temp=%g; fun_upper[j]=%g; fun_upper[j+1]=%g\n",
		j, k - 1, t, upper[j].z, upper[j + 1].z, s[j - 1], upper[j].slope, upper[j].absc, temp, fun_upper(k, upper[j].z, upper), fun_upper(k, upper[j + 1].z, upper));
	// else if ((j+1<k) && (t>upper[j+1].z)) std::cout << "nanu2";
//				;	char x; std::cin >> x;
		t = upper[j].z;
  		flag = true;
	}

	//	if (j == k - 1) std::cout << setw(20) << upper[j].z << setw(20) << t << setw(20) << upper[j].center << setw(5) << k << setw(5) << upper.size() << std::endl;
// END:;
	return t;
}




double dwiener_d(double q, double a, double vn, double wn, double leps) {
	double kll, kss, ans, v, w;
	double errziel = leps;
	double err = leps * 1.2;
	int zahl = 0;

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
	double lg1 = (-v * a * w - (pow(v, 2)) * q / 2) - 2 * std::log(a);

NEW:
	ans = 0;
	double es = (err - lg1);
	kss = ks(q_asq, w, es);
	double el = es;
	kll = kl(q_asq, v, w, el);


	if (2 * kss < kll)  // if small t is better
		ans = lg1 + logfs(q_asq, w, static_cast<int>(kss));
	else
	// if large t is better...
	ans = lg1 + logfl(q_asq, v, w, static_cast<int>(kll));

	zahl++; // MONITOR(0, 5)++;
	if (zahl == 10) {
		Rprintf("Zahl = 10 %20g%20g%20g%20g%20g\n", q, a, vn, wn, leps);
		return ans;
	}
	if (err - ans > errziel) {
//		MONITOR(1, 5)++;
		err = errziel*(1 + zahl*0.1) + ans;
		goto NEW;
	}

	return ans;
}

double dtdwiener_d(double q, double a, double v, double w, double &ld, double leps)
{
	double kll, kss, ans;
	//	double err = 1e-12;
	double errziel = leps;
	double err = errziel*1.2;
	int zahl = 0;
	//1.e-12

	double q_asq = q / pow(a, 2);

	double ans0 = -pow(v, 2) / 2.0;
	double la = 2.0 * std::log(a);
	double lg1 = -v * a * w - pow(v, 2) * q / 2.0 - la;
	double factor = lg1 - la;


NEW:
	double es = err - lg1;// + ld;
	es = es + la; //wegen q_asq in dtk1
	kss = dtaks(q_asq, w, es);
	double el = err - lg1;// + ld;
	el = el +la; 	// wegen q_asq in dtkl
	kll = dtakl(q_asq, v, a, el);

	// if small t is better
	if (2*kss < kll)
	{
		double erg; int newsign;
		logdtfs(q_asq, w, static_cast<int>(kss), erg, newsign);
		ans = ans0 - 1.5 / q + newsign * exp(factor - 1.5 * M_LN2 - M_LN_SQRT_PI - 3.5 * std::log(q_asq) + erg - ld);
	}
	// if large t is better...
	else
	{
		double erg; int newsign;
		logdtfl(q_asq, w, static_cast<int>(kll), erg, newsign);
		ans = ans0 - newsign * exp(factor + 3 * M_LNPI - M_LN2 + erg - ld);

	}

	zahl++;
	if(zahl == 10) {
		Rprintf("Zahl dt = 10\n");
		return ans;
	}

	double temp = std::log(fabs(ans)) + ld;
	if (temp < ld) {
		double check = temp - ld;
		if (err - check > errziel) {
			err = errziel*(1 + zahl*0.1) + check;
			goto NEW;
		}
	}
	double check = temp + M_LN2 - ld;

//	MONITOR(0, 6)++;

	if (err  + check > errziel) {
//		std::cout << setw(20) << err << setw(20) << check << std::endl;
//		MONITOR(1, 6)++;
		err = errziel*(1 + zahl*0.1) - check;
		ld = dwiener_d(-q, a, v, w, err);
		goto NEW;
	}
//	else std::cout << "okay" << std::endl;

	return ans;
}

int int_ddiff_d(unsigned dim, const double *x, void *p, unsigned fdim, double *retval) {
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
	// int K = (params->K);
	// int epsFLAG = (params->epsFLAG);
	// double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t - tau), a, nu, omega, errorW, 0, 1);
		// double ldW = dwiener_d(low_or_up * (t - tau), a, nu, omega, std::log(errorW));

		double temp2 = 0;
		if (sv) temp2 = -0.5 * pow(y, 2) - M_LN_SQRT_PI - 0.5 * M_LN2 + log1p(temp) - 2 * log1p(-temp);

		double integrand = exp(ldW + temp2);
// if (integrand == 0) integrand = exp(-700);
		retval[0] = integrand;
	}

	return 0;
}

double ddiff_d(double t, int low_or_up, double a, double v, double t0, double w, double sw, double sv, double st, double myerr) {

	double value;

	double *val_ptr = &value;

	int cnt = 0;

	double newerror = pow(myerr, 1.1);

NEW:

	double errorW = newerror ? newerror*0.1 : 1.e-12*0.1;

	my_params params = {t, low_or_up, a, v, t0, w, sw, sv, st, errorW, 0, 1, val_ptr};

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
	double abstol = newerror ? newerror*0.9 : 1.e-12*0.9;

	double val, err;

	int Meval = 6000;

	hcubature(1, int_ddiff_d, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err);
// Rprintf("val = %g", val);
	double logval = std::log(val);
	cnt++;
	if (cnt == 10) {
		Rprintf("cnt = 10 %20g%20g%20g%20g%20g\n", t, a, v, w, newerror*.1);
		free(xmin); free(xmax);
		return logval;
	}

	if (std::log(newerror) - logval > std::log(myerr)) {
		newerror = exp(std::log(myerr)*(1+cnt*0.1) + logval);
		goto NEW;
	}

	free(xmin); free(xmax);
	return logval;

}

int int_dtddiff_d(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
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
	// int K = (params->K);
	// int epsFLAG = (params->epsFLAG);
	// double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	double temp = sv ? pow(x[0], 2) : 0;
	double y = sv ? x[0] / (1 - temp) : 0;
	double nu = sv ? v + sv * y : v;
	double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
	double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, nu, omega, errorW, 0, 1);

		double temp2 = 0;
		if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		double dt;
		dtdwiener(t-tau, a, -low_or_up*nu, wn, ldW, &dt, errorW, 0, 1);

		double integrand = dt * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

double dtdiff_d(double t, int low_or_up, double a, double v, double t0, double w, double sw, double sv, double st, double myerr, double &d) {

	double value;

	double *val_ptr = &value;
	int cnt = 0;

	double newerror = pow(myerr, 1.1);

NEW:

	double errorW = newerror ? newerror*0.1 : 1.e-12*0.1;

	my_params params = {t, low_or_up, a, v, t0, w, sw, sv, st, errorW, 0, 1, val_ptr};

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
	double abstol = newerror ? newerror*0.9 : 1.e-12*0.9;

	double val, err;

	int Meval = 6000;

	hcubature(1, int_dtddiff_d, &params, dim, xmin, xmax, Meval, abstol, reltol, ERROR_INDIVIDUAL, &val, &err);
	// val = val;

	cnt++;
	if(cnt == 10) {
		Rprintf("cnt dt = 10\n");
		free(xmin); free(xmax);
		return val * exp(-d);
	}

	double temp = std::log(fabs(val));
	if (temp + std::log(myerr) < std::log(newerror) && newerror != 1.e-15) {
			newerror = exp(std::log(myerr)*(1+cnt*0.1) + temp);
			newerror = newerror < 1.e-15 ? 1.e-15 : newerror;
			goto NEW;
	}
	double check = temp + M_LN2 - d;

	if (std::log(newerror) + check > std::log(myerr) && newerror != 1.e-15) {
		double lnewerror = std::log(myerr)*(1+cnt*0.1) - check;
		newerror = exp(lnewerror);
		newerror = newerror < 1.e-15 ? 1.e-15 : newerror;
		d = ddiff_d(t, low_or_up, a, v, 0, w, sw, sv, 0, newerror);
		goto NEW;
	}


	free(xmin); free(xmax);
	return val * exp(-d);

}

void wiener_comp(double start, double scale, double norm, double alpha, double a, double v, double w, double sw, double sv, point& h) {
	h.x = alpha;
	double t = exp(alpha * scale + start);

	if (sv == 0.0 && sw == 0.0) {
		h.h = dwiener_d(-t, a, v, w, -27.63102);
		h.dh = dtdwiener_d(t, a, v, w, h.h, -23.02585);
	} else {
		h.h = ddiff_d(t, -1, a, v, 0.0, w, sw, sv, 0.0, 1.e-9);
		h.dh = dtdiff_d(t, -1, a, v, 0.0, w, sw, sv, 0.0, 1.e-7, h.h);
	}
if(h.h == -INFINITY) Rprintf("t = %g\n", t);
	h.h += start + h.x * scale + std::log(scale) - norm;
	h.dh = (1.0 + t * h.dh) * scale;
}

bool compare(point h1, point h2) {
	return (h1.x < h2.x);
}

void initialize_ars(double a, double v, double w, double sw, double sv, double bound, ars_archiv& ars_store) {

	double norm = 0.0, step = 1.0;
	//			double temp = fmax(0.2, fabs(v)) / a * w;
	//			double start = sqrt(temp * temp + 2.25) - 1.5;
	//			start = std::log(exp_mean(0, a, v, w) * start / temp);
				// based on mode of inverse-Gaussian
	double t0 = exp_mean(0, a, v, w);
	double start = std::log(t0);

	double scale;
	if (fabs(v) > 0.1)
		scale = sqrt(lower_bound_var(a, v, w)) / 1.0;
	else {
		int sign = (v > 0) ? 1 : -1;
		scale = sqrt(lower_bound_var(a, sign * 0.1, w)) / 1.0;
	}

	scale = fmax(scale, 0.05);
	scale = std::log(t0 + scale) - start;

	point begun, one;
	double ub, lb;
	point high, low;
	begun.x = 0.0;
	wiener_comp(start, scale, norm, begun.x, a, v, w, sw, sv, begun);
	if (fabs(begun.dh) > 10) { Rprintf("begun %5g%20g%20g%20g\n", begun.dh, a, v, w); }
	norm = begun.h; begun.h = 0.0;
	int sign = begun.dh > 0 ? 1 : -1;
	// make ldh >2.0 <5.0
	std::vector<point> h;
	for (int i = 0; i != 2; i++) {
		one = begun;
		double dh = sign * one.dh;
		if ((dh <= 2.0) || (dh >= 5.0))
		{
			if (dh <= 2.0) {
				while (dh <= 2.0) {
					double dold = dh;
					lb = one.x;
					one.x -= sign * step;

					wiener_comp(start, scale, norm, one.x, a, v, w, sw, sv, one);
					if ((abs(one.dh) > 2.0) && (abs(one.dh) < 5)) h.push_back(one);
					dh = sign * one.dh;
					if (dold > dh) {
						Rprintf("convex2?\n");
						//char xx; std::cin >> xx;
					}
				}
				ub = one.x;
			}
			else {
				while (dh >= 5.0) {
					double dold = dh;
					ub = one.x;
					one.x += sign * step;
					wiener_comp(start, scale, norm, one.x, a, v, w, sw, sv, one);
					if ((abs(one.dh) > 2.0) && (abs(one.dh) < 5)) h.push_back(one);
					dh = sign * one.dh;
					if (dh > dold) {
						Rprintf("convex2?\n");
						//						std::cin >> xx;
					}
				}
				//							}
				lb = one.x;
			}
		}
		while ((dh <= 2.0) || (dh >= 5.0)) {
			one.x = (lb + ub) / 2.0;
			wiener_comp(start, scale, norm, one.x, a, v, w, sw, sv, one);
			if (abs(one.dh) < 5) h.push_back(one);
			dh = sign * one.dh;
			if (dh <= 2.0) { lb = one.x; }
			if (dh >= 5.0) { ub = one.x; }
		}
		if (sign == 1) low = one; else high = one;
		sign = -sign;
	}

	wiener_comp(start, scale, norm, (low.x + high.x) / 2, a, v, w, sw, sv, one);


	h.push_back(one); h.push_back(low); h.push_back(high);

	bound = (std::log(bound) - start) / scale;
	if (high.x > bound) {
		one.x = bound - step;
		wiener_comp(start, scale, norm, one.x, a, v, w, sw, sv, one);
		h.push_back(one);
	}
	if (low.x > bound) {
		one.x = bound - 2 * step;
		wiener_comp(start, scale, norm, one.x, a, v, w, sw, sv, one);
		h.push_back(one);
	}

	std::sort(h.begin(), h.end(), compare);


	std::vector<point> xemp; xemp.clear();
	xemp.push_back(h[0]);
	for (long long unsigned int i = 1; i != h.size(); i++) if (h[i].x > h[i - 1].x) xemp.push_back(h[i]);
	h.assign(xemp.begin(), xemp.end());


	int k = 0;
	std::vector<piece> lower, upper; std::vector<double> s;
	generate_intervals(k, -INFINITY, h, lower, upper, s);

	ars_store.startstore = start;
	ars_store.scalestore = scale;
	ars_store.normstore = norm;
	ars_store.hstore.assign(h.begin(), h.end());
	ars_store.lowerstore.assign(lower.begin(), lower.end());
	ars_store.upperstore.assign(upper.begin(), upper.end());
	ars_store.sstore.assign(s.begin(), s.end());

	//			if (h.size() < 10) {
	//  		std::cout << setw(20) << a << setw(20) << v << setw(20) << w << setw(20) << scale << setw(20) << h.size() << std::endl;
	// 				char zz; std::cin >> zz;
	//			}
}




bool accept(double logf_t_prop,double c2) {
	if (c2 <= 0.06385320297074884) Rprintf("hm\n"); // std::log(5/3) / 16, req. for convergence
	double z=std::log(oneuni())+logf_t_prop;
	double b=-c2;
	int k=3; double lk=std::log(static_cast<double>(k));
	while (true)  {
		if (z>b) return false;
		b=logdiff(b,(lk-c2*k*k));
		if (z<b) return true;
		k=k+2; lk=std::log(static_cast<double>(k));
		b=logsum(b,(lk-c2*k*k));
		k=k+2; lk=std::log(static_cast<double>(k));
	}
}

double norm_exp_proposal(double drift) {
	double drift2 = pow(drift, 2);
	double rate = (4 * drift2 + M_PISQ) / 8.0;
	double a = (3.0 + sqrt(9.0 + 4 * drift2)) / 6.0;
	double p = exp(logsum(drift,-drift));
	double temp = -sqrt(fmax((a - 1) / a,0.0))*drift;
	double  cf1s = sqrt(a) * p * exp(temp), cf1l = 2 * M_PI * p / (4 * drift2 + M_PISQ);
	//	double t_prop;
	double prob_that_large = cf1l * exp(-rate * t_tilde);
	double prob_that_small = cf1s * gsl_sf_erfc(1.0 / sqrt(2 * a * t_tilde));

	while (true) {
		double z = oneuni() * (prob_that_small + prob_that_large);
		if (z <= prob_that_small) {
			double up = 1.0 / sqrt(a * t_tilde), t_prop = 1.0 / (a * pow(gsl_ran_ugaussian_tail(up), 2));
			double f_t_prop = (-1.0 / (2 * a * t_prop) + temp + drift2 * t_prop / 2.0);
			double c2 = 1.0 / (2 * t_prop);
			if (accept(f_t_prop, c2)) return t_prop;
		}
		else {
			double z = oneuni(); //zc=log1p(-z);
			double t_prop = t_tilde - std::log(z) / rate;
			double c2 = M_PISQ * t_prop / 8.0;  double f_t_prop = (-c2);
			if (accept(f_t_prop, c2)) return t_prop;
		}
	}
}

double invnorm(double drift) {
  double mu=1.0/drift;
  double result=0.0;
  double y;
    y=onenorm();
    y=pow(y,2);
	double w = mu * y;
    double x=mu + mu/2.0*(w-sqrt(w*(4.0+w)));
    double z=oneuni();
    if (z <= mu/(mu+x)) result=x;
	else result=mu*mu/x;
    return result ;
}

double invgauss_proposal(double drift) {
    double t_prop, cs, cl;
    if (t_tilde >= 0.63662) { cs = 0.0; cl = -M_LNPI + 2 * M_LN2 - 0.5*(M_LNPI + M_LN2); }
    else { cl = -M_PISQ / 8.0*t_tilde + 1.5*std::log(t_tilde) + 1.0 / (2 * t_tilde); cs = cl + 0.5*(M_LNPI + M_LN2) + M_LNPI - 2 * M_LN2; }

    while (true) {
        t_prop = invnorm(drift);


        if (t_prop <= t_tilde) { if (accept(cs - 1.0 / (2 * t_prop), 1.0 / (2 * t_prop))) return t_prop; }
        else if (accept(cl - 1 / (2 * t_prop) - 1.5 * std::log(t_prop), M_PISQ * t_prop / 8)) return t_prop;
    }

}

double rdiffusion(double drift, double a) {
	// a = threshold separation -> distance of upper/lower to zero = a/2
	a = a / 2.0; drift = drift * a; double a2 = pow(a, 2);
	if (fabs(drift) <= 1) return(a2*norm_exp_proposal(fabs(drift)));
	else return(a2*invgauss_proposal(fabs(drift)));
}

double rdiffusion_UPbound(double bound, double a, double drift, double w)
{
	double b_lo = -w * a; double b_up = (1 - w)*a;
START:
	double x = 0.0, t = 0.0;
	while (true) {
		const double xlo = fabs(x - b_lo);
		const double xup = fabs(x - b_up);
		if (fabs(xlo - xup) < 1e-5) {
			// symmetric bounds, diffusion model in [x - xup, x + xup]
			t += rdiffusion(drift, (2 * xup));
			if (t > bound) goto START;
			if (oneuni() < 1 / (1 + exp(-2 * drift * xup)))
				return t;
			else
				return -t;
// symmetric; need lower bound
		}
		else if (xlo > xup) {
			// x closer to upper bound, diffusion model in [x - xup, x + xup]
			t += rdiffusion(drift, (2 * xup));
			if (t > bound) goto START;
			if (oneuni() < 1 / (1 + exp(-2 * drift*xup)))
				return t;
			x -= xup;
			// 			bound-=step;
		}
		else {
			// x closer to lower bound, diffusion model in [x-xlo, x+xlo]
			t += rdiffusion(drift, (2 * xlo));
			if (t > bound) goto START;
			if (oneuni() > 1 / (1 + exp(-2 * drift*xlo)))
				return -t;
			x += xlo;
		}
	}
}



double rdiffusion_lower_trunc(double bound, double a, double drift, double w)
{
	double b_lo = -w * a; double b_up = (1 - w) * a;
START:
	double x = 0.0, t = 0.0;
	while (true) {
		const double xlo = fabs(x - b_lo);
		const double xup = fabs(x - b_up);
		if (fabs(xlo - xup) < 1e-5) {
			// symmetric bounds, diffusion model in [x - xup, x + xup]
			t += rdiffusion(drift, (2 * xup));
			if (t > bound) goto START;
//			if (oneuni(rst) < 1 / (1 + exp(-2 * drift * xup)))
//				goto START;
			//				return t;
			return -t;
			// symmetric; need lower bound
		}	else if (xlo > xup) {
			// x closer to upper bound, diffusion model in [x - xup, x + xup]
			t+= rdiffusion(drift, (2 * xup));
			if (t > bound) goto START;
			if (oneuni() < 1 / (1 + exp(-2 * drift * xup))) goto START;
//				return t;
			x -= xup;
			// 			bound-=step;
		}	else {
			// x closer to lower bound, diffusion model in [x-xlo, x+xlo]
			t+=rdiffusion(drift, (2 * xlo));
			if (t > bound) goto START;
			if (oneuni() > 1 / (1 + exp(-2 * drift * xlo)))
				return -t;
			x += xlo;
		}
	}
}

double rwiener_diag2(int pm, double bound, double a, double v, double w, double err, int K, int epsFLAG) {
	double eps = 1e-5;
	double pmid = 0;
	double qmin = 0;

	double qmax = bound;
	double q = std::isinf(bound) ? 1.0 : bound / 2;
	double p = std::log(oneuni());

	double qold;

	if (pm == 1) {
		v = -v;
		w = 1.0 - w;
	}
  double total = pwiener(bound, a, v, w, err, K, epsFLAG);
//	double total = logFjoint_lower(bound, a, v, w);

	do {
		qold = q;
    pmid = pwiener(q, a, v, w, err, K, epsFLAG) - total;
		//pmid = logFjoint_lower(q, a, v, w) - total;
		if (p <= pmid) {
			qmax = q;
			q = qmin + (qmax - qmin) / 2.0;
		}
		else {
			qmin = q;
			q = std::isinf(qmax) ? q * 2 : qmin + (qmax - qmin) / 2.0;
		}

	} while (fabs(q - qold) > eps);
	return(q);
}
