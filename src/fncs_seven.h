
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#ifndef FNCS_SEVEN_H
#define FNCS_SEVEN_H



/* dependencies for Wrapper */
int int_ddiff(unsigned, const double *, void *, unsigned, double *);
int int_daddiff(unsigned, const double *, void *, unsigned, double *);
int int_dvddiff(unsigned, const double *, void *, unsigned, double *);
int int_dt0ddiff(unsigned, const double *, void *, unsigned, double *);
int int_dwddiff(unsigned, const double *, void *, unsigned, double *);
int int_dswddiff(unsigned, const double *, void *, unsigned, double *);
int int_dsvddiff(unsigned, const double *, void *, unsigned, double *);
int int_dst0ddiff(unsigned, const double *, void *, unsigned, double *);
int int_pdiff(unsigned, const double *, void *, unsigned, double *);
int int_dapdiff(unsigned, const double *, void *, unsigned, double *);
int int_dvpdiff(unsigned, const double *, void *, unsigned, double *);
int int_dt0pdiff(unsigned, const double *, void *, unsigned, double *);
int int_dwpdiff(unsigned, const double *, void *, unsigned, double *);
int int_dswpdiff(unsigned, const double *, void *, unsigned, double *);
int int_dsvpdiff(unsigned, const double *, void *, unsigned, double *);
int int_dst0pdiff(unsigned, const double *, void *, unsigned, double *);
/*-----------------------------------------------*/


/* Wrapper */
void ddiff(int, double, int, double, double, double, double, double, double, double, double, int, int, double *, double *);
void pdiff(int, double, int, double, double, double, double, double, double, double, double, int, int, double *, double *);
/*-----------------------------------------------*/



#endif
