
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#ifndef CDF_FNCS_H
#define CDF_FNCS_H



/* dependencies for DISTRIBUTION */
double logP(int, double, double, double);
double Ks(double, double, double, double, double);
double Kl(double, double, double, double, double);
double logFs(double, double, double, double, int);
double logFl(double, double, double, double, int);
/*-----------------------------------------------*/

/* dependencies for d/da DISTRIBUTION */
double davlogP(int, double, double, double);
double dalogP(int, double, double, double, double);
void dakS(double, double, double, double, double, double &);
void dakL(double, double, double, double, double, double &);
void logdaFs(int, int, double, double, double, double, double &, double);
void logdaFl(int, int, double, double, double, double, double &, double);
/*-----------------------------------------------*/

/* dependencies for d/dv DISTRIBUTION */
double dvlogP(int, double, double, double, double);
void dvkS(double, double, double, double, double, double &);
void dvkL(double, double, double, double, double, double &);
void logdvFs(int, int, double, double, double, double, double &, double);
void logdvFl(int, int, double, double, double, double, double &, double);
/*-----------------------------------------------*/

/* dependencies for d/dw DISTRIBUTION */
double dwlogP(int, double, double, double);
void dwkS(double, double, double, double, double, double &);
void dwkL(double, double, double, double, double, double &);
void logdwFs(int, int, double, double, double, double, double &, double);
void logdwFl(int, int, double, double, double, double, double &, double);
/*-----------------------------------------------*/

/* dependencies for d/dt d/da d/dv d/dw DISTRIBUTION */
void dxkS(double, double, double, double, double, double &);
void dxkL(double, double, double, double, double, double &);
void logdxFs(int, int, int, int, double, double, double, double, double, double &, double &, double &);
void logdxFl(int, int, int, int, double, double, double, double, double, double &, double &, double &);
/*-----------------------------------------------*/

/* DISTRIBUTION */
double pwiener(double, double, double, double, double, int, int);
/*-----------------------------------------------*/

/* d/da DISTRIBUTION */
void dapwiener(int, double, double, double, double, double, double *, double, int, int);
/*-----------------------------------------------*/

/* d/dv DISTRIBUTION */
void dvpwiener(int, double, double, double, double, double, double *, double, int, int);
/*-----------------------------------------------*/

/* d/dw DISTRIBUTION */
void dwpwiener(int, double, double, double, double, double, double *, double, int, int);
/*-----------------------------------------------*/

/* d/dt d/da d/dv d/dw DISTRIBUTION */
void dxpwiener(int, double, double, double, double, double, double, int, int, double *, double *, double *);
/*-----------------------------------------------*/


#endif
