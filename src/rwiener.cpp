

#include "rwiener.h"
#include <cmath>

#define t_tilde 2.5


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
	double lg1 = (-v * a * w - (pow(v, 2)) * q / 2) - 2 * log(a);

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
	double la = 2.0 * log(a);
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
		ans = ans0 - 1.5 / q + newsign * exp(factor - 1.5 * M_LN2 - M_LN_SQRT_PI - 3.5 * log(q_asq) + erg - ld);
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

	double temp = log(fabs(ans)) + ld;
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
		// double ldW = dwiener_d(low_or_up * (t - tau), a, nu, omega, log(errorW));

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
	double logval = log(val);
	cnt++;
	if (cnt == 10) {
		Rprintf("cnt = 10 %20g%20g%20g%20g%20g\n", t, a, v, w, newerror*.1);
		free(xmin); free(xmax);
		return logval;
	}

	if (log(newerror) - logval > log(myerr)) {
		newerror = exp(log(myerr)*(1+cnt*0.1) + logval);
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

	double temp = log(fabs(val));
	if (temp + log(myerr) < log(newerror) && newerror != 1.e-15) {
			newerror = exp(log(myerr)*(1+cnt*0.1) + temp);
			newerror = newerror < 1.e-15 ? 1.e-15 : newerror;
			goto NEW;
	}
	double check = temp + M_LN2 - d;

	if (log(newerror) + check > log(myerr) && newerror != 1.e-15) {
		double lnewerror = log(myerr)*(1+cnt*0.1) - check;
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
	h.h += start + h.x * scale + log(scale) - norm;
	h.dh = (1.0 + t * h.dh) * scale;
}

bool compare(point h1, point h2) {
	return (h1.x < h2.x);
}

void initialize_ars(double a, double v, double w, double sw, double sv, double bound, ars_archiv& ars_store) {

	double norm = 0.0, step = 1.0;
	//			double temp = fmax(0.2, fabs(v)) / a * w;
	//			double start = sqrt(temp * temp + 2.25) - 1.5;
	//			start = log(exp_mean(0, a, v, w) * start / temp);
				// based on mode of inverse-Gaussian
	double t0 = exp_mean(0, a, v, w);
	double start = log(t0);

	double scale;
	if (fabs(v) > 0.1)
		scale = sqrt(lower_bound_var(a, v, w)) / 1.0;
	else {
		int sign = (v > 0) ? 1 : -1;
		scale = sqrt(lower_bound_var(a, sign * 0.1, w)) / 1.0;
	}

	scale = fmax(scale, 0.05);
	scale = log(t0 + scale) - start;

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

	bound = (log(bound) - start) / scale;
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
	h = xemp;


	int k = 0;
	std::vector<piece> lower, upper; std::vector<double> s;
	generate_intervals(k, -INFINITY, h, lower, upper, s);

	ars_store.startstore = start;
	ars_store.scalestore = scale;
	ars_store.normstore = norm;
	ars_store.hstore = h;
	ars_store.lowerstore = lower;
	ars_store.upperstore = upper;
	ars_store.sstore = s;

	//			if (h.size() < 10) {
	//  		std::cout << setw(20) << a << setw(20) << v << setw(20) << w << setw(20) << scale << setw(20) << h.size() << std::endl;
	// 				char zz; std::cin >> zz;
	//			}
}




bool accept(double logf_t_prop,double c2) {
	if (c2 <= 0.06385320297074884) Rprintf("hm\n"); // log(5/3) / 16, req. for convergence
	double z=log(oneuni())+logf_t_prop;
	double b=-c2;
	int k=3; double lk=log(1.0*k);
	while (true)  {
		if (z>b) return false;
		b=logdiff(b,(lk-c2*k*k));
		if (z<b) return true;
		k=k+2; lk=log(1.0*k);
		b=logsum(b,(lk-c2*k*k));
		k=k+2; lk=log(1.0*k);
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
			double t_prop = t_tilde - log(z) / rate;
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
    else { cl = -M_PISQ / 8.0*t_tilde + 1.5*log(t_tilde) + 1.0 / (2 * t_tilde); cs = cl + 0.5*(M_LNPI + M_LN2) + M_LNPI - 2 * M_LN2; }

    while (true) {
        t_prop = invnorm(drift);


        if (t_prop <= t_tilde) { if (accept(cs - 1.0 / (2 * t_prop), 1.0 / (2 * t_prop))) return t_prop; }
        else if (accept(cl - 1 / (2 * t_prop) - 1.5 * log(t_prop), M_PISQ * t_prop / 8)) return t_prop;
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
	double p = log(oneuni());

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

double arst(ars_archiv& ars_store, double scale, double totallow, double start, double bound, double a, double v, double w, double sw, double sv,
	void generic2(double start, double scale, double norm, double alpha, double a, double v, double w, double sw, double sv, point& h)) {
	// intialize h's
	double norm = ars_store.normstore; bool flag;

// NEW:
	flag = false;
	std::vector<point> h; std::vector<piece> lower, upper; std::vector<double> s;
	double ww, tt, xstar;
	point one;


	h = ars_store.hstore;
	lower = ars_store.lowerstore;
	upper = ars_store.upperstore;
	s = ars_store.sstore;

	if (s.size() != h.size())
		Rprintf("Problem in ars\n");
	int k = static_cast<int>(h.size());
	bool update = false;

	if (bound < INFINITY)
	{
		int l = 0;
		while ((l != k) && (h[l].x < bound)) l++;
		k = l;
	}
	if (k < 2) {
		Rprintf("k Probel %g %d\n", exp(start + bound / scale), k);
		xstar = -INFINITY;
		goto END;
	}


	double ss;

WEITER:
//	MONITOR(0, 2)++;
	xstar = oneuni();

	xstar = inverse_distribution(k, xstar, upper, s, bound, flag);
	if (flag) {
//		std::cout << "NEW0 in ars";
//		std::cout << "arst "  << setw(20) << a << setw(20) << v << setw(20) << w << std::endl;
		xstar = -INFINITY;
		goto END;
	}

	ww = log(oneuni()); tt = fun_upper(k, xstar, upper);  ss = fun_lower(k, xstar, h, lower);

	if(xstar > 12) Rprintf("ww = %g   tt = %g   ss = %g\n", ww, tt, ss);
	if (ww <= (ss - tt))  goto STOP;
	one.x = xstar; generic2(start, scale, norm, xstar, a, v, w, sw, sv, one);
	if(xstar > 12) Rprintf("bin hier\n");
	if (ww <= (one.h - tt))  goto STOP;
	flag = update_intervals(k, totallow, one, h, lower, upper, s);
if(xstar > 12) Rprintf("flag = %d\n", flag);
	if (flag) {
		xstar = -INFINITY;
		goto END;
	}
	k = k + 1;
	update = true;

//	MONITOR(1, 2)++;
	goto WEITER;
STOP:
	if (update) {
		ars_store.hstore = h;
		ars_store.lowerstore = lower;
		ars_store.upperstore = upper;
		ars_store.sstore = s;
	}
END:	return xstar;

}

double make_rwiener2(ars_archiv& ars_store, double bound, double a, double v, double w, double sw, double sv, double err, int K, int epsFLAG) {
	double temp;
NEW:
	double start = ars_store.startstore;
// Rprintf("ars hstore length = %d", static_cast<int>(ars_store.hstore.size()));
	double  scale = ars_store.scalestore;

	double bound2 = (bound == INFINITY) ? bound : (log(bound) - start) / scale;

	//	start = start * scale;
	temp = arst(ars_store, scale, -INFINITY, start, bound2, a, v, w, sw, sv, wiener_comp);

	if (temp != -INFINITY) temp = exp(start + temp * scale);
	//	else temp = -rdiffusion_lower_trunc(bound, v, w, a, rst);
	else
	{
		Rprintf("ars hat nicht geklappt\n");
		initialize_ars(a, v, w, sw, sv, bound, ars_store);
		goto NEW;
		// double vs = v, ws = w;
		// bool REPEAT;
		// if(sv || sw) {
		// 	REPEAT = true;
		// 	while (REPEAT) {
		// 		vs = v, ws = w;
		// 		if (sv) vs += sv * onenorm();
		// 		if (sw) ws += sw * (oneuni()-0.5);
		// 		double P = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
		// 		REPEAT = P < oneuni();
		// 	}
		// }
		// temp = rwiener_diag2(0, bound, a, vs, ws, err, K, epsFLAG);
	}

	return temp;
}

// R=0 is lower bound and R=1 is upper bound
void run_make_rwiener(int choice, int N, double a, double v, double w, double sv, double sw, int R, double bound, double err, int K, int epsFLAG, double *q, int *resp, ars_archiv *ars_store1, ars_archiv *ars_store2, int use_store) {
	double vs, ws;
	// printf("h1 size = %d\n", static_cast<int>(ars_store1->hstore.size()));
	// printf("h2 size = %d\n", static_cast<int>(ars_store2->hstore.size()));
	switch (choice) {
		case 1:
			if (R != 0) { // one-sided

				// ars_archiv ars_store;
				if(R == 2) {
					v = -v;
					w = 1-w;
				}
				// initialize_ars(a, v, w, sw, sv, bound, ars_store);
				if (!use_store) initialize_ars(a, v, w, sw, sv, bound, *ars_store1);
				for (int i = 0; i != N; i++) {
					// q[i] = make_rwiener2(ars_store, bound, a, v, w, sw, sv, err, K, epsFLAG);
					q[i] = make_rwiener2(*ars_store1, bound, a, v, w, sw, sv, err, K, epsFLAG);
					resp[i] = R;
				}

			} else { // R = 0 -- both sides

				// ars_archiv ars_store1;
				// ars_archiv ars_store2;
				double p_up, Rerr = 99.9; //, p_lo;
				if (std::isfinite(bound)) { // truncated
					if (sv || sw) {
						double temp_u, temp_l;
	          pdiff(0, bound, 1.0, a, v, 0, w, sw, sv, 0, err, K, epsFLAG, &temp_u, &Rerr);
						pdiff(0, bound, -1.0, a, v, 0, w, sw, sv, 0, err, K, epsFLAG, &temp_l, &Rerr);
						p_up = temp_u / (temp_u + temp_l);
						// p_lo = 1-p_up;
					} else {
	          double temp_u = exp(pwiener(bound, a, -v, 1-w, err, K, epsFLAG));
						double temp_l = exp(pwiener(bound, a, v, w, err, K, epsFLAG));
						p_up = temp_u / (temp_u + temp_l);
						// p_lo = 1-p_up;
					}
				} else { // not truncated
					if (sv || sw) {
	          pdiff(0, bound, 1.0, a, v, 0, w, sw, sv, 0, err, K, epsFLAG, &p_up, &Rerr);
						// p_lo = 1-p_up;
					} else {
	          // p_up = exp(pwiener(bound, a, -v, 1-w, err, K, epsFLAG));
						p_up = (1-exp(2*v*a*w))/(exp(-2*v*a*(1-w))-exp(2*v*a*w));
						// p_lo = 1-p_up;
					}
				}
				int cnt_up = 0;
				for (int i = 0; i !=N; i++) {
					cnt_up += oneuni()<=p_up ? 1 : 0;
				}

				// printf("ars h1 size = %d\n", static_cast<int>(ars_store1->hstore.size()));
				if (!use_store) initialize_ars(a, -v, 1-w, sw, sv, bound, *ars_store1);
				for (int i = 0; i != cnt_up; i++) {
					q[i] = make_rwiener2(*ars_store1, bound, a, -v, 1-w, sw, sv, err, K, epsFLAG);
					resp[i] = 2;
				}

				// printf("ars h2 size = %d\n", static_cast<int>(ars_store2->hstore.size()));
				if (!use_store) initialize_ars(a, v, w, sw, sv, bound, *ars_store2);
				for (int i = cnt_up; i != N; i++) {
					q[i] = make_rwiener2(*ars_store2, bound, a, v, w, sw, sv, err, K, epsFLAG);
					resp[i] = 1;
				}

			}
			break;

		case 2:
			if(R != 0) { // one-sided
				double P;
				bool REPEAT;
				for (int i = 0; i != N; i++) {
					vs = v; ws = w;
					if(sv || sw) {
						REPEAT = true;
						while (REPEAT) {
							vs = v, ws = w;
							if (sv) vs += sv * onenorm();
							if (sw) ws += sw * (oneuni()-0.5);
							if (R == 2) {vs = -vs; ws = 1 - ws;}
							if (std::isfinite(bound)) P = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
							else P = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
							REPEAT = oneuni() > P;
						}
					} else {
						if (R == 2) {
							vs = -vs;
							ws = 1 - ws;
						}
					}
					q[i] = -rdiffusion_lower_trunc(bound, a, vs, ws);
					resp[i] = R;
				}

			} else { // R = 0 -- both sides
				if (std::isfinite(bound)) { // truncated
					double p_up, p_lo;
					bool REPEAT;
					for (int i = 0; i != N; i++) {
						vs = v; ws = w;
						if(sv || sw) {
							REPEAT = true;
							while (REPEAT) {
								vs = v; ws = w;
								if (sv) vs += sv*onenorm();
								if (sw) ws += sw*(oneuni()-0.5);
								p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
								p_up = exp(pwiener(bound, a, -vs, 1-ws, err, K, epsFLAG));
								REPEAT = oneuni() > (p_up + p_lo);
							}
						}
						q[i] = rdiffusion_UPbound(bound, a, vs, ws);
						resp[i] = q[i] > 0 ? 2 : 1;
						if (resp[i] == 1) q[i] = fabs(q[i]);
					}
				} else { // not truncated
					for (int i = 0; i != N; i++) {
						vs = v; ws = w;
						if (sv) vs += sv*onenorm();
						if (sw) ws += sw*(oneuni()-0.5);
						q[i] = rdiffusion_UPbound(bound, a, vs, ws);
						resp[i] = q[i] > 0 ? 2 : 1;
						if (resp[i] == 1) q[i] = fabs(q[i]);
					}
				}
			}
			break;

		case 3:
			double p_lo, p_up;
			if (R != 0) { // one-sided
				vs = v; ws = w;
				bool REPEAT;
				for (int i = 0; i != N; i++) {
					if(sv || sw) {
						REPEAT = true;
						while (REPEAT) {
							vs = v; ws = w;
							if (sv) vs += sv * onenorm();
							if (sw) ws += sw * (oneuni()-0.5);
							if (std::isfinite(bound)) { // truncated
								if (R == 1) {
									p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
									REPEAT = oneuni() > p_lo;
								}
								if (R == 2) {
									p_up = exp(pwiener(bound, a, -vs, 1-ws, err, K, epsFLAG));
									REPEAT = oneuni() > p_up;
								}
							} else { // not truncated
								// p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
								p_lo = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
								if (R == 1) REPEAT = oneuni() > p_lo;
								if (R == 2) REPEAT = oneuni() > 1-p_lo;
							}
						}
					}
					q[i] = rwiener_diag2(R-1, bound, a, vs, ws, err, K, epsFLAG);
					resp[i] = R;
				}

			} else { // R = 0 -- both sides

				int up_or_down;
				if (std::isfinite(bound)) { // truncated
					bool REPEAT;
					double u;
					for (int i = 0; i != N; i++) {
						vs = v; ws = w;
						if(sv || sw) {
							REPEAT = true;
							while (REPEAT) {
								vs = v; ws = w;
								if (sv) vs += sv * onenorm();
								if (sw) ws += sw * (oneuni()-0.5);
								p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
								p_up = exp(pwiener(bound, a, -vs, 1-ws, err, K, epsFLAG));
								u = oneuni();
								if (u <= p_lo) {
									REPEAT = false;
									up_or_down = 0;
								} else if (u >= 1-p_up) {
									REPEAT = false;
									up_or_down = 1;
								} else {
									REPEAT = true;
								}
							}
						} else {
							p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
							p_up = exp(pwiener(bound, a, -vs, 1-ws, err, K, epsFLAG));
							up_or_down = oneuni() <= p_up / (p_up + p_lo) ? 1 : 0;
						}
						q[i] = rwiener_diag2(up_or_down, bound, a, vs, ws, err, K, epsFLAG);
						resp[i] = up_or_down + 1;
					}
				} else { // not truncated
					for (int i = 0; i != N; i++) {
						vs = v; ws = w;
						if (sv) vs += sv * onenorm();
						if (sw) ws += sw * (oneuni()-0.5);
						// p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
						p_lo = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
						up_or_down = oneuni() < p_lo ? 0 : 1;
						q[i] = rwiener_diag2(up_or_down, bound, a, vs, ws, err, K, epsFLAG);
						resp[i] = up_or_down + 1;
					}
				}

			}
			break;

	}
}
