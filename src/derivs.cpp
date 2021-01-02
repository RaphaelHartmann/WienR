
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#include "cstdio"
#include "pdf_fncs.h"
#include "cdf_fncs.h"
#include "fncs_seven.h"
#include "tools.h"
#include "Rinternals.h"
#include <thread>
#include <vector>

/* PDF and CDF of Wiener diffusion */
  /* PDF of Wiener diffusion */
void PDF(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *Rpdf, double *Rlogpdf, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          double pm = (resp[i]==1) ? 1.0 : -1.0;
          // double mp = -1.0*pm;
          // double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
          double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
          Rlogpdf[i] = ld;
          Rpdf[i] = exp(ld);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      // double mp = -1.0*pm;
      // double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      Rlogpdf[i] = ld;
      Rpdf[i] = exp(ld);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      // double mp = -1.0*pm;
      // double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      Rlogpdf[i] = ld;
      Rpdf[i] = exp(ld);
    }
  }

}

  /* CDF of Wiener diffusion */
void CDF(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *Rcdf, double *Rlogcdf, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          int pm = (resp[i]==1) ? 1 : -1;
          int mp = -1*pm;
          double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
          Rlogcdf[i] = lp;
          Rcdf[i] = exp(lp);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      int pm = (resp[i]==1) ? 1 : -1;
      int mp = -1*pm;
      double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      Rlogcdf[i] = lp;
      Rcdf[i] = exp(lp);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      int pm = (resp[i]==1) ? 1 : -1;
      int mp = -1*pm;
      double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      Rlogcdf[i] = lp;
      Rcdf[i] = exp(lp);
    }
  }

}
/* ------------------------------------------------ */


/* partial derivs with respect to each param of PDF */
  /* derivative of PDF with respect to t */
void dtPDF(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *Rderiv, double *Rderiv_ln, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          double pm = (resp[i]==1) ? 1.0 : -1.0;
          double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
          double mp = -1.0*pm;
          dtdwiener(t[i], a[i], v[i]*mp, (resp[i]-w[i])*pm, ld, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      double mp = -1.0*pm;
      dtdwiener(t[i], a[i], v[i]*mp, (resp[i]-w[i])*pm, ld, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      double mp = -1.0*pm;
      dtdwiener(t[i], a[i], v[i]*mp, (resp[i]-w[i])*pm, ld, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
    }
  }

}

  /* derivative of PDF with respect to a */
void daPDF(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *Rderiv, double *Rderiv_ln, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          double pm = (resp[i]==1) ? 1.0 : -1.0;
          double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
          dadwiener(t[i]*pm, a[i], v[i], w[i], ld, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      dadwiener(t[i]*pm, a[i], v[i], w[i], ld, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      dadwiener(t[i]*pm, a[i], v[i], w[i], ld, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
    }
  }

}

  /* derivative of PDF with respect to v */
void dvPDF(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *Rderiv, double *Rderiv_ln, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          double pm = (resp[i]==1) ? 1.0 : -1.0;
          double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
          dvdwiener(t[i]*pm, a[i], v[i], w[i], ld, &Rderiv[i], &Rderiv_ln[i]);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      dvdwiener(t[i]*pm, a[i], v[i], w[i], ld, &Rderiv[i], &Rderiv_ln[i]);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      dvdwiener(t[i]*pm, a[i], v[i], w[i], ld, &Rderiv[i], &Rderiv_ln[i]);
    }
  }

}

  /* derivative of PDF with respect to w */
void dwPDF(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *Rderiv, double *Rderiv_ln, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          double pm = (resp[i]==1) ? 1.0 : -1.0;
          double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
          dwdwiener(t[i]*pm, a[i], v[i], w[i], ld, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      dwdwiener(t[i]*pm, a[i], v[i], w[i], ld, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      dwdwiener(t[i]*pm, a[i], v[i], w[i], ld, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
    }
  }

}
/* ------------------------------------------------ */


/* partial derivs with respect to each param of CDF */
  /* derivative of CDF with respect to a */
void daCDF(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *Rderiv, double *Rderiv_ln, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          int pm = (resp[i]==1) ? 1 : -1;
          int mp = -1*pm;
          double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
          dapwiener(pm, t[i], a[i], v[i], w[i], lp, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      int pm = (resp[i]==1) ? 1 : -1;
      int mp = -1*pm;
      double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      dapwiener(pm, t[i], a[i], v[i], w[i], lp, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      int pm = (resp[i]==1) ? 1 : -1;
      int mp = -1*pm;
      double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      dapwiener(pm, t[i], a[i], v[i], w[i], lp, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
    }
  }

}

  /* derivative of CDF with respect to v */
void dvCDF(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *Rderiv, double *Rderiv_ln, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          int pm = (resp[i]==1) ? 1 : -1;
          int mp = -1*pm;
          double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
          dvpwiener(pm, t[i], a[i], v[i], w[i], lp, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      int pm = (resp[i]==1) ? 1 : -1;
      int mp = -1*pm;
      double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      dvpwiener(pm, t[i], a[i], v[i], w[i], lp, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      int pm = (resp[i]==1) ? 1 : -1;
      int mp = -1*pm;
      double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      dvpwiener(pm, t[i], a[i], v[i], w[i], lp, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
    }
  }

}

  /* derivative of CDF with respect to w */
void dwCDF(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *Rderiv, double *Rderiv_ln, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          int pm = (resp[i]==1) ? 1 : -1;
          int mp = -1*pm;
          double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
          dwpwiener(pm, t[i], a[i], v[i], w[i], lp, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      int pm = (resp[i]==1) ? 1 : -1;
      int mp = -1*pm;
      double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      dwpwiener(pm, t[i], a[i], v[i], w[i], lp, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      int pm = (resp[i]==1) ? 1 : -1;
      int mp = -1*pm;
      double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      dwpwiener(pm, t[i], a[i], v[i], w[i], lp, &Rderiv[i], &Rderiv_ln[i], eps, K, epsFLAG);
    }
  }

}
/* ------------------------------------------------ */


/* gradient of PDF and CDF */
  /* gradient of PDF */
void dxPDF_old(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *da, double *da_ln, double *dv, double *dv_ln, double *dw, double *dw_ln, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          double pm = (resp[i]==1) ? 1.0 : -1.0;
          double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
          dxdwiener(t[i]*pm, a[i], v[i], w[i], ld, eps, K, epsFLAG, &da[i], &da_ln[i], &dv[i], &dv_ln[i], &dw[i], &dw_ln[i]);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      dxdwiener(t[i]*pm, a[i], v[i], w[i], ld, eps, K, epsFLAG, &da[i], &da_ln[i], &dv[i], &dv_ln[i], &dw[i], &dw_ln[i]);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      dxdwiener(t[i]*pm, a[i], v[i], w[i], ld, eps, K, epsFLAG, &da[i], &da_ln[i], &dv[i], &dv_ln[i], &dw[i], &dw_ln[i]);
    }
  }

}

void dxPDF(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *da, double *da_ln, double *dv, double *dv_ln, double *dw, double *dw_ln, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          double pm = (resp[i]==1) ? 1.0 : -1.0;
          double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
          dadwiener(t[i]*pm, a[i], v[i], w[i], ld, &da[i], &da_ln[i], eps, K, epsFLAG);
          dvdwiener(t[i]*pm, a[i], v[i], w[i], ld, &dv[i], &dv_ln[i]);
          dwdwiener(t[i]*pm, a[i], v[i], w[i], ld, &dw[i], &dw_ln[i], eps, K, epsFLAG);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      dadwiener(t[i]*pm, a[i], v[i], w[i], ld, &da[i], &da_ln[i], eps, K, epsFLAG);
      dvdwiener(t[i]*pm, a[i], v[i], w[i], ld, &dv[i], &dv_ln[i]);
      dwdwiener(t[i]*pm, a[i], v[i], w[i], ld, &dw[i], &dw_ln[i], eps, K, epsFLAG);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], eps, K, epsFLAG);
      dadwiener(t[i]*pm, a[i], v[i], w[i], ld, &da[i], &da_ln[i], eps, K, epsFLAG);
      dvdwiener(t[i]*pm, a[i], v[i], w[i], ld, &dv[i], &dv_ln[i]);
      dwdwiener(t[i]*pm, a[i], v[i], w[i], ld, &dw[i], &dw_ln[i], eps, K, epsFLAG);
    }
  }

}

  /* gradient of CDF */
void dxCDF_old(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *da, double *da_ln, double *dv, double *dv_ln, double *dw, double *dw_ln, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          int pm = (resp[i]==1) ? 1 : -1;
          int mp = -1*pm;
          double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
          dxpwiener(pm, t[i], a[i], v[i], w[i], lp, eps, K, epsFLAG, &da[i], &da_ln[i], &dv[i], &dv_ln[i], &dw[i], &dw_ln[i]);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      int pm = (resp[i]==1) ? 1 : -1;
      int mp = -1*pm;
      double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      dxpwiener(pm, t[i], a[i], v[i], w[i], lp, eps, K, epsFLAG, &da[i], &da_ln[i], &dv[i], &dv_ln[i], &dw[i], &dw_ln[i]);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      int pm = (resp[i]==1) ? 1 : -1;
      int mp = -1*pm;
      double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      dxpwiener(pm, t[i], a[i], v[i], w[i], lp, eps, K, epsFLAG, &da[i], &da_ln[i], &dv[i], &dv_ln[i], &dw[i], &dw_ln[i]);
    }
  }

}

void dxCDF(double *t, double *a, double *v, double *w, double eps, int *resp, int K, int N, int epsFLAG, double *da, double *da_ln, double *dv, double *dv_ln, double *dw, double *dw_ln, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          int pm = (resp[i]==1) ? 1 : -1;
          int mp = -1*pm;
          double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
          dapwiener(pm, t[i], a[i], v[i], w[i], lp, &da[i], &da_ln[i], eps, K, epsFLAG);
          dvpwiener(pm, t[i], a[i], v[i], w[i], lp, &dv[i], &dv_ln[i], eps, K, epsFLAG);
          dwpwiener(pm, t[i], a[i], v[i], w[i], lp, &dw[i], &dw_ln[i], eps, K, epsFLAG);
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      int pm = (resp[i]==1) ? 1 : -1;
      int mp = -1*pm;
      double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      dapwiener(pm, t[i], a[i], v[i], w[i], lp, &da[i], &da_ln[i], eps, K, epsFLAG);
      dvpwiener(pm, t[i], a[i], v[i], w[i], lp, &dv[i], &dv_ln[i], eps, K, epsFLAG);
      dwpwiener(pm, t[i], a[i], v[i], w[i], lp, &dw[i], &dw_ln[i], eps, K, epsFLAG);
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      int pm = (resp[i]==1) ? 1 : -1;
      int mp = -1*pm;
      double lp = pwiener(t[i], a[i], v[i]*mp, pm*(resp[i]-w[i]), eps, K, epsFLAG);
      dapwiener(pm, t[i], a[i], v[i], w[i], lp, &da[i], &da_ln[i], eps, K, epsFLAG);
      dvpwiener(pm, t[i], a[i], v[i], w[i], lp, &dv[i], &dv_ln[i], eps, K, epsFLAG);
      dwpwiener(pm, t[i], a[i], v[i], w[i], lp, &dw[i], &dw_ln[i], eps, K, epsFLAG);
    }
  }

}
/* ------------------------------------------------ */



/* PDF and CDF of 7-param diffusion */
  /* PDF of 7-param diffusion */
void PDF7(int choice, double *t, int *resp, double *a, double *v, double *t0, double *w, double *sw, double *sv, double *st, double err, int K, int N, int epsFLAG, double *Rval, double *Rlogval, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
          ddiff(choice, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rval[i], &Rlogval[i]);
          if (choice == 0) {
            Rlogval[i] = log(Rval[i]);
          }
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
      ddiff(choice, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rval[i], &Rlogval[i]);
      if (choice == 0) {
        Rlogval[i] = log(Rval[i]);
      }
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
      ddiff(choice, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rval[i], &Rlogval[i]);
      if (choice == 0) {
        Rlogval[i] = log(Rval[i]);
      }
    }
  }

}



  /* CDF of 7-param diffusion */
void CDF7(int choice, double *t, int *resp, double *a, double *v, double *t0, double *w, double *sw, double *sv, double *st, double err, int K, int N, int epsFLAG, double *Rval, double *Rlogval, int NThreads) {

  if (NThreads) {
    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* calculate derivative with parallelization */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
          pdiff(choice, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rval[i], &Rlogval[i]);
          if (choice == 0) {
            Rlogval[i] = log(Rval[i]);
          }
        }
      });
    }

    int last = NperThread * (AmntOfThreads-1);
    for (int i = last; i < N; i++) {
      double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
      pdiff(choice, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rval[i], &Rlogval[i]);
      if (choice == 0) {
        Rlogval[i] = log(Rval[i]);
      }
    }

    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } else {
    /* calculate derivative without parallelization */
    for(int i = 0; i < N; i++) {
      double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
      pdiff(choice, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rval[i], &Rlogval[i]);
      if (choice == 0) {
        Rlogval[i] = log(Rval[i]);
      }
    }
  }
}

  /* gradient of CDF */
void dxPDF7(double *t, int *resp, double *a, double *v, double *t0, double *w, double *sw, double *sv, double *st, double err, int K, int N, int epsFLAG, double *Rda, double *Rda_ln, double *Rdv, double *Rdv_ln, double *Rdt0, double *Rdt0_ln, double *Rdw, double *Rdw_ln, double *Rdsw, double *Rdsw_ln, double *Rdsv, double *Rdsv_ln, double *Rdst, double *Rdst_ln, int NThreads) {

    if (NThreads) {
      /* prepare threads */
      int maxThreads = std::thread::hardware_concurrency();
      if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
      int suppThreads = maxThreads == 0 ? 2 : maxThreads;
      int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
      int NperThread = N / AmntOfThreads;
      std::vector<std::thread> threads(AmntOfThreads-1);

      /* calculate derivative with parallelization */
      for (int j = 0; j < AmntOfThreads-1; j++) {
        threads[j] = std::thread([=]() {
          for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
            double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
            ddiff(1, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rda[i], &Rda_ln[i]);
            ddiff(2, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdv[i], &Rdv_ln[i]);
            ddiff(3, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdt0[i], &Rdt0_ln[i]);
            ddiff(4, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdw[i], &Rdw_ln[i]);
            if (sw[0]) ddiff(5, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdsw[i], &Rdsw_ln[i]);
            else Rdsw[i] = Rdsw_ln[i] = NAN;
            if (sv[0]) ddiff(6, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdsv[i], &Rdsv_ln[i]);
            else Rdsv[i] = Rdsv_ln[i] = NAN;
            if (st[0]) ddiff(7, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdst[i], &Rdst_ln[i]);
            else Rdst[i] = Rdst_ln[i] = NAN;
          }
        });
      }

      int last = NperThread * (AmntOfThreads-1);
      for (int i = last; i < N; i++) {
        double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
        ddiff(1, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rda[i], &Rda_ln[i]);
        ddiff(2, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdv[i], &Rdv_ln[i]);
        ddiff(3, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdt0[i], &Rdt0_ln[i]);
        ddiff(4, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdw[i], &Rdw_ln[i]);
        if (sw[0]) ddiff(5, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdsw[i], &Rdsw_ln[i]);
        else Rdsw[i] = Rdsw_ln[i] = NAN;
        if (sv[0]) ddiff(6, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdsv[i], &Rdsv_ln[i]);
        else Rdsv[i] = Rdsv_ln[i] = NAN;
        if (st[0]) ddiff(7, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdst[i], &Rdst_ln[i]);
        else Rdst[i] = Rdst_ln[i] = NAN;
      }

      for (int j = 0; j < AmntOfThreads-1; j++) {
        threads[j].join();
      }

    } else {
      /* calculate derivative without parallelization */
      for(int i = 0; i < N; i++) {
        double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
        ddiff(1, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rda[i], &Rda_ln[i]);
        ddiff(2, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdv[i], &Rdv_ln[i]);
        ddiff(3, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdt0[i], &Rdt0_ln[i]);
        ddiff(4, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdw[i], &Rdw_ln[i]);
        if (sw[0]) ddiff(5, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdsw[i], &Rdsw_ln[i]);
        else Rdsw[i] = Rdsw_ln[i] = NAN;
        if (sv[0]) ddiff(6, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdsv[i], &Rdsv_ln[i]);
        else Rdsv[i] = Rdsv_ln[i] = NAN;
        if (st[0]) ddiff(7, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdst[i], &Rdst_ln[i]);
        else Rdst[i] = Rdst_ln[i] = NAN;
      }
    }

  }

void dxCDF7(double *t, int *resp, double *a, double *v, double *t0, double *w, double *sw, double *sv, double *st, double err, int K, int N, int epsFLAG, double *Rda, double *Rda_ln, double *Rdv, double *Rdv_ln, double *Rdt0, double *Rdt0_ln, double *Rdw, double *Rdw_ln, double *Rdsw, double *Rdsw_ln, double *Rdsv, double *Rdsv_ln, double *Rdst, double *Rdst_ln, int NThreads) {

    if (NThreads) {
      /* prepare threads */
      int maxThreads = std::thread::hardware_concurrency();
      if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
      int suppThreads = maxThreads == 0 ? 2 : maxThreads;
      int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
      int NperThread = N / AmntOfThreads;
      std::vector<std::thread> threads(AmntOfThreads-1);

      /* calculate derivative with parallelization */
      for (int j = 0; j < AmntOfThreads-1; j++) {
        threads[j] = std::thread([=]() {
          for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
            double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
            pdiff(1, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rda[i], &Rda_ln[i]);
            pdiff(2, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdv[i], &Rdv_ln[i]);
            pdiff(3, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdt0[i], &Rdt0_ln[i]);
            pdiff(4, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdw[i], &Rdw_ln[i]);
            if (sw[0]) pdiff(5, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdsw[i], &Rdsw_ln[i]);
            else Rdsw[i] = Rdsw_ln[i] = NAN;
            if (sv[0]) pdiff(6, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdsv[i], &Rdsv_ln[i]);
            else Rdsv[i] = Rdsv_ln[i] = NAN;
            if (st[0]) pdiff(7, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdst[i], &Rdst_ln[i]);
            else Rdst[i] = Rdst_ln[i] = NAN;
          }
        });
      }

      int last = NperThread * (AmntOfThreads-1);
      for (int i = last; i < N; i++) {
        double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
        pdiff(1, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rda[i], &Rda_ln[i]);
        pdiff(2, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdv[i], &Rdv_ln[i]);
        pdiff(3, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdt0[i], &Rdt0_ln[i]);
        pdiff(4, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdw[i], &Rdw_ln[i]);
        if (sw[0]) pdiff(5, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdsw[i], &Rdsw_ln[i]);
        else Rdsw[i] = Rdsw_ln[i] = NAN;
        if (sv[0]) pdiff(6, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdsv[i], &Rdsv_ln[i]);
        else Rdsv[i] = Rdsv_ln[i] = NAN;
        if (st[0]) pdiff(7, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdst[i], &Rdst_ln[i]);
        else Rdst[i] = Rdst_ln[i] = NAN;
      }

      for (int j = 0; j < AmntOfThreads-1; j++) {
        threads[j].join();
      }

    } else {
      /* calculate derivative without parallelization */
      for(int i = 0; i < N; i++) {
        double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
        pdiff(1, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rda[i], &Rda_ln[i]);
        pdiff(2, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdv[i], &Rdv_ln[i]);
        pdiff(3, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdt0[i], &Rdt0_ln[i]);
        pdiff(4, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdw[i], &Rdw_ln[i]);
        if (sw[0]) pdiff(5, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdsw[i], &Rdsw_ln[i]);
        else Rdsw[i] = Rdsw_ln[i] = NAN;
        if (sv[0]) pdiff(6, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdsv[i], &Rdsv_ln[i]);
        else Rdsv[i] = Rdsv_ln[i] = NAN;
        if (st[0]) pdiff(7, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, &Rdst[i], &Rdst_ln[i]);
        else Rdst[i] = Rdst_ln[i] = NAN;
      }
    }

  }
