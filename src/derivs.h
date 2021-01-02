
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#ifndef DERIVS_H
#define DERIVS_H

/* PDF and CDF of Wiener diffusion */
  /* PDF */
  void PDF(double *, double *, double *, double *, double, int *, int, int, int, double *, double *, int);
  /* CDF */
  void CDF(double *, double *, double *, double *, double, int *, int, int, int, double *, double *, int);
/* ------------------------------------------------ */

/* partial derivs with respect to each param of PDF */
  /* derivative of PDF with respect to t */
  void dtPDF(double *, double *, double *, double *, double, int *, int, int, int, double *, int);
  /* derivative of PDF with respect to a */
  void daPDF(double *, double *, double *, double *, double, int *, int, int, int, double *, int);
  /* derivative of PDF with respect to v */
  void dvPDF(double *, double *, double *, double *, double, int *, int, int, int, double *, int);
  /* derivative of PDF with respect to w */
  void dwPDF(double *, double *, double *, double *, double, int *, int, int, int, double *, int);
/* ------------------------------------------------ */

/* partial derivs with respect to each param of CDF */
  /* derivative of CDF with respect to a */
  void daCDF(double *, double *, double *, double *, double, int *, int, int, int, double *, int);
  /* derivative of CDF with respect to v */
  void dvCDF(double *, double *, double *, double *, double, int *, int, int, int, double *, int);
  /* derivative of CDF with respect to w */
  void dwCDF(double *, double *, double *, double *, double, int *, int, int, int, double *, int);
/* ------------------------------------------------ */

/* gradient of PDF and CDF */
  /* derivative of PDF with respect to all params */
  void dxPDF_old(double *, double *, double *, double *, double, int *, int, int, int, double *, double *, double *, int);
  void dxPDF(double *, double *, double *, double *, double, int *, int, int, int, double *, double *, double *, int);
  /* derivative of CDF with respect to all params */
  void dxCDF_old(double *, double *, double *, double *, double, int *, int, int, int, double *, double *, double *, int);
  void dxCDF(double *, double *, double *, double *, double, int *, int, int, int, double *, double *, double *, int);
/* ------------------------------------------------ */


/* PDF and CDF of 7-param diffusion */
  /* PDF */
  void PDF7(int, double *, int *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, double *, double *, int);

  /* CDF */
  void CDF7(int, double *, int *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, double *, double *, int);
/* ------------------------------------------------ */

/* gradient of PDF and CDF */
  /* derivative of PDF with respect to all params */
  void dxPDF7(double *, int *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, double *, double *, double *, double *, double *, double *, double *, int);
  /* derivative of CDF with respect to all params */
  void dxCDF7(double *, int *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, double *, double *, double *, double *, double *, double *, double *, int);
/* ------------------------------------------------ */






#endif
