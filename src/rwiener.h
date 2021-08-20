


#ifndef RWIENER_H
#define RWIENER_H

#include "tools.h"

double bla(double a, double b, double c);
double fun_upper(int k, double x, std::vector<piece> upper);
void generate_intervals(int& k, double totallow, std::vector<point> h, std::vector<piece>& lower, std::vector<piece>& upper, std::vector<double>& s);
bool update_intervals(int k, double totallow, point new_point, std::vector<point>& h, std::vector<piece>& lower, std::vector<piece>& upper, std::vector<double>& s);
double fun_lower(int k, double x, std::vector<point> h, std::vector<piece> lower);
double inverse_distribution(int k, double xstar, std::vector<piece> upper, std::vector<double> s, double bound, bool& flag);

double dwiener_d(double q, double a, double vn, double wn, double leps);
double dtdwiener_d(double q, double a, double v, double w, double &ld, double leps);
int int_ddiff_d(unsigned dim, const double *x, void *p, unsigned fdim, double *retval);
double ddiff_d(double t, int low_or_up, double a, double v, double t0, double w, double sw, double sv, double st, double myerr);
int int_dtddiff_d(unsigned dim, const double* x, void* p, unsigned fdim, double* retval);
double dtdiff_d(double t, int low_or_up, double a, double v, double t0, double w, double sw, double sv, double st, double myerr, double &d);

void wiener_comp(double start, double scale, double norm, double alpha, double a, double v, double w, double sw, double sv, point& h);
bool compare(point h1, point h2);
void initialize_ars(double a, double v, double w, double sw, double sv, double bound, ars_archiv& ars_store);

bool accept(double logf_t_prop,double c2);

double norm_exp_proposal(double drift);
double invnorm(double drift);
double invgauss_proposal(double drift);

double rdiffusion(double drift, double a);
double rdiffusion_UPbound(double bound, double a, double drift, double w);
double rdiffusion_lower_trunc(double bound, double a, double drift, double w);
double rwiener_diag2(int pm, double bound, double a, double v, double w, double err, int K, int epsFLAG);

#endif
