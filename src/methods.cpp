
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#include "tools.h"
#include "cdf_fncs.h"
#include "fncs_seven.h"
#include "rwiener.h"
#include <cmath>
#include <thread>
#include <mutex>
#include <atomic>


std::mutex mtx_ars1, mtx_ars2; // mutexes for global ars_archiv's
std::atomic<int> atm_ars1 (3), atm_ars2 (3); // atomic counters



double arst(ars_archiv& ars_store, ars_archiv* ars_str_glob, double scale, double totallow, double start, double bound, double a, double v, double w, double sw, double sv, int Nars_parallel,
	void generic2(double start, double scale, double norm, double alpha, double a, double v, double w, double sw, double sv, point& h)) {
	// intialize h's
	double norm = ars_store.normstore; bool flag;

// NEW:
	flag = false;
	std::vector<point> h; std::vector<piece> lower, upper; std::vector<double> s;
	double ww, tt, ss, xstar;
	point one;

	if (Nars_parallel == 1) {
		if (static_cast<int>(ars_store.hstore.size()) < atm_ars1.load()) {
			std::lock_guard<std::mutex> guard(mtx_ars1);
			ars_store.hstore.assign(ars_str_glob->hstore.begin(), ars_str_glob->hstore.end());
			ars_store.lowerstore.assign(ars_str_glob->lowerstore.begin(), ars_str_glob->lowerstore.end());
			ars_store.upperstore.assign(ars_str_glob->upperstore.begin(), ars_str_glob->upperstore.end());
			ars_store.sstore.assign(ars_str_glob->sstore.begin(), ars_str_glob->sstore.end());
		}
	}
	if (Nars_parallel == 2) {
		if (static_cast<int>(ars_store.hstore.size()) < atm_ars2.load()) {
			std::lock_guard<std::mutex> guard(mtx_ars2);
			ars_store.hstore.assign(ars_str_glob->hstore.begin(), ars_str_glob->hstore.end());
			ars_store.lowerstore.assign(ars_str_glob->lowerstore.begin(), ars_str_glob->lowerstore.end());
			ars_store.upperstore.assign(ars_str_glob->upperstore.begin(), ars_str_glob->upperstore.end());
			ars_store.sstore.assign(ars_str_glob->sstore.begin(), ars_str_glob->sstore.end());
		}
	}

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

WEITER:

	xstar = oneuni();

	xstar = inverse_distribution(k, xstar, upper, s, bound, flag);
	if (flag) {
		xstar = -INFINITY;
		goto END;
	}

	ww = std::log(oneuni()); tt = fun_upper(k, xstar, upper);  ss = fun_lower(k, xstar, h, lower);

	if (ww <= (ss - tt))  goto STOP;
	one.x = xstar; generic2(start, scale, norm, xstar, a, v, w, sw, sv, one);
	if (ww <= (one.h - tt))  goto STOP;


	if (Nars_parallel) { // mutex needed for parallelization

		std::vector<point> h_g; std::vector<piece> lower_g, upper_g; std::vector<double> s_g;

		if (Nars_parallel == 1) { // lockguard to guard global ars_archiv 1

			std::lock_guard<std::mutex> guard(mtx_ars1);
			// Rprintf("in guard\n");
			h_g.assign(ars_str_glob->hstore.begin(), ars_str_glob->hstore.end());
			lower_g.assign(ars_str_glob->lowerstore.begin(), ars_str_glob->lowerstore.end());
			upper_g.assign(ars_str_glob->upperstore.begin(), ars_str_glob->upperstore.end());
			s_g.assign(ars_str_glob->sstore.begin(), ars_str_glob->sstore.end());

			flag = update_intervals(k, totallow, one, h_g, lower_g, upper_g, s_g);
			if (flag) {
				xstar = -INFINITY;
				return xstar;
			}

			ars_str_glob->hstore.assign(h_g.begin(), h_g.end());
			ars_str_glob->lowerstore.assign(lower_g.begin(), lower_g.end());
			ars_str_glob->upperstore.assign(upper_g.begin(), upper_g.end());
			ars_str_glob->sstore.assign(s_g.begin(), s_g.end());
		}

		if (Nars_parallel == 2) { // lockguard to guard global ars_archiv 2

			std::lock_guard<std::mutex> guard(mtx_ars2);
			// Rprintf("in guard\n");
			h_g = ars_str_glob->hstore;
			lower_g = ars_str_glob->lowerstore;
			upper_g = ars_str_glob->upperstore;
			s_g = ars_str_glob->sstore;

			flag = update_intervals(k, totallow, one, h_g, lower_g, upper_g, s_g);
			if (flag) {
				xstar = -INFINITY;
				return xstar;
			}

			ars_str_glob->hstore.assign(h_g.begin(), h_g.end());
			ars_str_glob->lowerstore.assign(lower_g.begin(), lower_g.end());
			ars_str_glob->upperstore.assign(upper_g.begin(), upper_g.end());
			ars_str_glob->sstore.assign(s_g.begin(), s_g.end());
		}

		if (Nars_parallel == 1) {
			atm_ars1.store(static_cast<int>(h_g.size()));
		}
		if (Nars_parallel == 2) {
			atm_ars2.store(static_cast<int>(h_g.size()));
		}

		k = static_cast<int>(h_g.size());
		h.assign(h_g.begin(), h_g.end());
		lower.assign(lower_g.begin(), lower_g.end());
		upper.assign(upper_g.begin(), upper_g.end());
		s.assign(s_g.begin(), s_g.end());

		if (bound < INFINITY)
		{
			int l = 0;
			while ((l != k) && (h_g[l].x < bound)) l++;
			k = l;
		}

		update = true;

	} else { // no mutex needed since no parallelization

		flag = update_intervals(k, totallow, one, h, lower, upper, s);
		if (flag) {
			xstar = -INFINITY;
			goto END;
		}
		k = k + 1;
		update = true;

	}

	goto WEITER;

STOP:
	if (update) {
		ars_store.hstore.assign(h.begin(), h.end());
		ars_store.lowerstore.assign(lower.begin(), lower.end());
		ars_store.upperstore.assign(upper.begin(), upper.end());
		ars_store.sstore.assign(s.begin(), s.end());
	}

END:
	return xstar;

}


double make_rwiener2(ars_archiv& ars_store, ars_archiv* ars_str_glob, double bound, double a, double v, double w, double sw, double sv, double err, int K, int epsFLAG, int Nars_parallel) {
	double temp;
NEW:
	double start = ars_store.startstore;
// Rprintf("ars hstore length = %d", static_cast<int>(ars_store.hstore.size()));
	double  scale = ars_store.scalestore;

	double bound2 = (bound == INFINITY) ? bound : (std::log(bound) - start) / scale;

	//	start = start * scale;
	temp = arst(ars_store, ars_str_glob, scale, -INFINITY, start, bound2, a, v, w, sw, sv, Nars_parallel, wiener_comp);

	if (temp != -INFINITY) temp = exp(start + temp * scale);
	//	else temp = -rdiffusion_lower_trunc(bound, v, w, a, rst);
	else
	{
		Rprintf("ars hat nicht geklappt\n");
		initialize_ars(a, v, w, sw, sv, bound, ars_store);
		goto NEW;
	}

	return temp;
}





void method1_one(int N, double a, double v, double w, double t0, double sv, double sw, double st0, int R, double bound, double err, int K, int epsFLAG, int NThreads, double *q, int *resp, ars_archiv *ars_store1, ars_archiv *ars_store2, int use_store) {

	if(R == 2) {
		v = -v;
		w = 1-w;
	}

  /* END "WITH PARALLELIZATION" */
  if (NThreads) {

		/* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

		/* prepare global ars_archiv */
		if (!use_store) initialize_ars(a, v, w, sw, sv, bound-t0, *ars_store1);
		ars_archiv ars_str_tmp = *ars_store1;

    /* starting threads while ... */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
				ars_archiv ars_str_local = ars_str_tmp;
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
          double tau = 0;
          if (t0) {
            tau = t0;
            if (st0) tau += st0*oneuni();
          }
					q[i] = tau + make_rwiener2(ars_str_local, ars_store1, bound-t0, a, v, w, sw, sv, err, K, epsFLAG, 1);
					resp[i] = R;
        }
      });
    }

    /* the main thread also runs */
    int last = NperThread * (AmntOfThreads-1);
		ars_archiv ars_str_local = ars_str_tmp;
    for (int i = last; i < N; i++) {
      // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
      double tau = 0;
      if (t0) {
        tau = t0;
        if (st0) tau += st0*oneuni();
      }
			q[i] = tau + make_rwiener2(ars_str_local, ars_store1, bound-t0, a, v, w, sw, sv, err, K, epsFLAG, 1);
			resp[i] = R;
    }

    /* join threads */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } /* END "WITH PARALLELIZATION" */

  /* START "WITH _NO_ PARALLELIZATION" */
  else {

		if (!use_store) initialize_ars(a, v, w, sw, sv, bound-t0, *ars_store1);
		for (int i = 0; i != N; i++) {
			if (i % 1024 == 0) R_CheckUserInterrupt();
			double tau = 0;
			if (t0) {
			  tau = t0;
			  if (st0) tau += st0*oneuni();
			}
			q[i] = tau + make_rwiener2(*ars_store1, nullptr, bound-t0, a, v, w, sw, sv, err, K, epsFLAG, 0);
			resp[i] = R;
		}

  } /* END "WITH _NO_ PARALLELIZATION" */

}

void method1_both(int N, double a, double v, double w, double t0, double sv, double sw, double st0, int R, double bound, double err, int K, int epsFLAG, int NThreads, double *q, int *resp, ars_archiv *ars_store1, ars_archiv *ars_store2, int use_store) {

  bool truncated = R_FINITE(bound), sv_or_sw = (sv > 0 || sw > 0);
	double Rerr = 99.9;
	int Neval = 6000;

  /* END "WITH PARALLELIZATION" */
  if (NThreads) {

    double p_up;
		if (truncated) { // truncated
			if (sv_or_sw) {
				double temp_u, temp_l;
				pdiff(0, bound, 1.0, a, v, 0, w, sw, sv, 0, err, K, epsFLAG, Neval, &temp_u, &Rerr);
				pdiff(0, bound, -1.0, a, v, 0, w, sw, sv, 0, err, K, epsFLAG, Neval, &temp_l, &Rerr);
				p_up = temp_u / (temp_u + temp_l);
			} else {
				double temp_u = exp(pwiener(bound, a, -v, 1-w, err, K, epsFLAG));
				double temp_l = exp(pwiener(bound, a, v, w, err, K, epsFLAG));
				p_up = temp_u / (temp_u + temp_l);
			}
		} else { // not truncated
			if (sv_or_sw) {
				pdiff(0, bound, 1.0, a, v, 0, w, sw, sv, 0, err, K, epsFLAG, Neval, &p_up, &Rerr);
			} else {
				p_up = (1-exp(2*v*a*w))/(exp(-2*v*a*(1-w))-exp(2*v*a*w));
			}
		}
		int cnt_up = 0;
		for (int i = 0; i !=N; i++) {
			cnt_up += oneuni()<=p_up ? 1 : 0;
		}

		/* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread1 = cnt_up / AmntOfThreads;
		int NperThread2 = (N-cnt_up) / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

		/* prepare global ars_archives */
		if (!use_store) {
			initialize_ars(a, -v, 1-w, sw, sv, bound-t0, *ars_store1);
			initialize_ars(a, v, w, sw, sv, bound-t0, *ars_store2);
		}
		ars_archiv ars_str_tmp1 = *ars_store1;
		ars_archiv ars_str_tmp2 = *ars_store2;

    /* starting threads while ... */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
				ars_archiv ars_str_local1 = ars_str_tmp1;
        for (int i = j*NperThread1; i < (j+1)*NperThread1; i++) {
          // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
          double tau = 0;
          if (t0) {
            tau = t0;
            if (st0) tau += st0*oneuni();
          }
					q[i] = tau + make_rwiener2(ars_str_local1, ars_store1, bound-t0, a, -v, 1-w, sw, sv, err, K, epsFLAG, 1);
					resp[i] = 2;
        }
				ars_archiv ars_str_local2 = ars_str_tmp2;
				for (int i = j*NperThread2+cnt_up; i < (j+1)*NperThread2+cnt_up; i++) {
				  // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
				  double tau = 0;
				  if (t0) {
				    tau = t0;
				    if (st0) tau += st0*oneuni();
				  }
					q[i] = tau + make_rwiener2(ars_str_local2, ars_store2, bound-t0, a, v, w, sw, sv, err, K, epsFLAG, 2);
					resp[i] = 1;
        }
      });
    }

    /* the main thread also runs */
		int last1 = NperThread1 * (AmntOfThreads-1);
		int last2 = NperThread2 * (AmntOfThreads-1)+cnt_up;
		ars_archiv ars_str_local1 = ars_str_tmp1;
		for (int i = last1; i < cnt_up; i++) {
		  // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
		  double tau = 0;
		  if (t0) {
		    tau = t0;
		    if (st0) tau += st0*oneuni();
		  }
			q[i] = tau + make_rwiener2(ars_str_local1, ars_store1, bound-t0, a, -v, 1-w, sw, sv, err, K, epsFLAG, 1);
			resp[i] = 2;
    }
		ars_archiv ars_str_local2 = ars_str_tmp2;
		for (int i = last2; i < N; i++) {
		  // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
		  double tau = 0;
		  if (t0) {
		    tau = t0;
		    if (st0) tau += st0*oneuni();
		  }
			q[i] = tau + make_rwiener2(ars_str_local2, ars_store2, bound-t0, a, v, w, sw, sv, err, K, epsFLAG, 2);
			resp[i] = 1;
    }

    /* join threads */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } /* END "WITH PARALLELIZATION" */

  /* START "WITH _NO_ PARALLELIZATION" */
  else {

		double p_up, Rerr = 99.9;//, p_lo;
		if (truncated) { // truncated
			if (sv_or_sw) {
				double temp_u, temp_l;
				pdiff(0, bound, 1.0, a, v, 0, w, sw, sv, 0, err, K, epsFLAG, Neval, &temp_u, &Rerr);
				pdiff(0, bound, -1.0, a, v, 0, w, sw, sv, 0, err, K, epsFLAG, Neval, &temp_l, &Rerr);
				p_up = temp_u / (temp_u + temp_l);
				// p_lo = 1-p_up;
			} else {
				double temp_u = exp(pwiener(bound, a, -v, 1-w, err, K, epsFLAG));
				double temp_l = exp(pwiener(bound, a, v, w, err, K, epsFLAG));
				p_up = temp_u / (temp_u + temp_l);
				// p_lo = 1-p_up;
			}
		} else { // not truncated
			if (sv_or_sw) {
				pdiff(0, bound, 1.0, a, v, 0, w, sw, sv, 0, err, K, epsFLAG, Neval, &p_up, &Rerr);
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

		if (!use_store) initialize_ars(a, -v, 1-w, sw, sv, bound-t0, *ars_store1);
		for (int i = 0; i != cnt_up; i++) {
			if (i % 1024 == 0) R_CheckUserInterrupt();
			double tau = 0;
			if (t0) {
			  tau = t0;
			  if (st0) tau += st0*oneuni();
			}
			q[i] = tau + make_rwiener2(*ars_store1, nullptr, bound-t0, a, -v, 1-w, sw, sv, err, K, epsFLAG, 0);
			resp[i] = 2;
		}

		if (!use_store) initialize_ars(a, v, w, sw, sv, bound-t0, *ars_store2);
		for (int i = cnt_up; i != N; i++) {
			if (i % 1024 == 0) R_CheckUserInterrupt();
			double tau = 0;
			if (t0) {
			  tau = t0;
			  if (st0) tau += st0*oneuni();
			}
			q[i] = tau + make_rwiener2(*ars_store2, nullptr, bound-t0, a, v, w, sw, sv, err, K, epsFLAG, 0);
			resp[i] = 1;
		}

  } /* END "WITH _NO_ PARALLELIZATION" */

}


void method2_one(int N, double a, double v, double w, double t0, double sv, double sw, double st0, int R, double bound, double err, int K, int epsFLAG, int NThreads, double *q, int *resp) {

  bool truncated = R_FINITE(bound), sv_or_sw = (sv > 0 || sw > 0);

  /* END "WITH PARALLELIZATION" */
  if (NThreads) {

    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* starting threads while ... */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        double P, vs = v, ws = w;
        bool REPEAT;
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
          vs = v; ws = w;
          if(sv_or_sw) {
            REPEAT = true;
            while (REPEAT) {
              vs = v, ws = w;
              if (sv) vs += sv * onenorm();
              if (sw) ws += sw * (oneuni()-0.5);
              if (R == 2) {vs = -vs; ws = 1 - ws;}
              if (truncated) { // truncated
                P = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
              } else { // not truncated
                P = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
              }
              REPEAT = oneuni() > P;
            }
          } else {
            if (R == 2) {
              vs = -vs;
              ws = 1 - ws;
            }
          }
          q[i] = -rdiffusion_lower_trunc(bound, a, vs, ws, t0, st0);
          resp[i] = R;
        }
      });
    }

    /* the main thread also runs */
    int last = NperThread * (AmntOfThreads-1);
    double P, vs = v, ws = w;
    bool REPEAT;
    for (int i = last; i < N; i++) {
      // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
      vs = v; ws = w;
      if(sv_or_sw) {
        REPEAT = true;
        while (REPEAT) {
          vs = v, ws = w;
          if (sv) vs += sv * onenorm();
          if (sw) ws += sw * (oneuni()-0.5);
          if (R == 2) {vs = -vs; ws = 1 - ws;}
          if (truncated) { // truncated
            P = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
          } else { // not truncated
            P = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
          }
          REPEAT = oneuni() > P;
        }
      } else {
        if (R == 2) {
          vs = -vs;
          ws = 1 - ws;
        }
      }
      q[i] = -rdiffusion_lower_trunc(bound, a, vs, ws, t0, st0);
      resp[i] = R;
    }

    /* join threads */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } /* END "WITH PARALLELIZATION" */

  /* START "WITH _NO_ PARALLELIZATION" */
  else {

    double P, vs = v, ws = w;
    bool REPEAT;
    for (int i = 0; i != N; i++) {
			if (i % 1024 == 0) R_CheckUserInterrupt();
      vs = v; ws = w;
      if(sv_or_sw) {
        REPEAT = true;
        while (REPEAT) {
          vs = v, ws = w;
          if (sv) vs += sv * onenorm();
          if (sw) ws += sw * (oneuni()-0.5);
          if (R == 2) {vs = -vs; ws = 1 - ws;}
          if (truncated) { // truncated
            P = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
          } else { // not truncated
            P = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
          }
          REPEAT = oneuni() > P;
        }
      } else {
        if (R == 2) {
          vs = -vs;
          ws = 1 - ws;
        }
      }
      q[i] = -rdiffusion_lower_trunc(bound, a, vs, ws, t0, st0);
      resp[i] = R;
    }

  } /* END "WITH _NO_ PARALLELIZATION" */

}

void method2_both(int N, double a, double v, double w, double t0, double sv, double sw, double st0, int R, double bound, double err, int K, int epsFLAG, int NThreads, double *q, int *resp) {

  bool truncated = R_FINITE(bound), sv_or_sw = (sv > 0 || sw > 0);

  /* END "WITH PARALLELIZATION" */
  if (NThreads) {

    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* starting threads while ... */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        double p_up, p_lo, vs = v, ws = w;
        bool REPEAT;
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
          vs = v; ws = w;
          if (truncated) {
            if(sv_or_sw) {
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
          } else {
            if (sv) vs += sv*onenorm();
            if (sw) ws += sw*(oneuni()-0.5);
          }
          q[i] = rdiffusion_UPbound(bound, a, vs, ws, t0, st0);
          resp[i] = q[i] > 0 ? 2 : 1;
          if (resp[i] == 1) q[i] = fabs(q[i]);
        }
      });
    }

    /* the main thread also runs */
    int last = NperThread * (AmntOfThreads-1);
    double p_up, p_lo, vs = v, ws = w;
    bool REPEAT;
    for (int i = last; i < N; i++) {
      // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
      vs = v; ws = w;
      if (truncated) {
        if(sv_or_sw) {
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
      } else {
        if (sv) vs += sv*onenorm();
        if (sw) ws += sw*(oneuni()-0.5);
      }
      q[i] = rdiffusion_UPbound(bound, a, vs, ws, t0, st0);
      resp[i] = q[i] > 0 ? 2 : 1;
      if (resp[i] == 1) q[i] = fabs(q[i]);
    }

    /* join threads */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } /* END "WITH PARALLELIZATION" */

  /* START "WITH _NO_ PARALLELIZATION" */
  else {

    double p_up, p_lo, vs = v, ws = w;
    bool REPEAT;
    for (int i = 0; i != N; i++) {
			if (i % 1024 == 0) R_CheckUserInterrupt();
      vs = v; ws = w;
      if (truncated) {
        if(sv_or_sw) {
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
      } else {
        if (sv) vs += sv*onenorm();
        if (sw) ws += sw*(oneuni()-0.5);
      }
      q[i] = rdiffusion_UPbound(bound, a, vs, ws, t0, st0);
      resp[i] = q[i] > 0 ? 2 : 1;
      if (resp[i] == 1) q[i] = fabs(q[i]);
    }

  } /* END "WITH _NO_ PARALLELIZATION" */

}


void method3_one(int N, double a, double v, double w, double t0, double sv, double sw, double st0, int R, double bound, double err, int K, int epsFLAG, int NThreads, double *q, int *resp) {

  bool truncated = R_FINITE(bound), sv_or_sw = (sv > 0 || sw > 0);

  /* END "WITH PARALLELIZATION" */
  if (NThreads) {

    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* starting threads while ... */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        double p_up, p_lo, vs = v, ws = w;
        bool REPEAT;
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
          if(sv_or_sw) {

            REPEAT = true;
            while (REPEAT) {
              vs = v; ws = w;
              if (sv) vs += sv * onenorm();
              if (sw) ws += sw * (oneuni()-0.5);
              if (truncated) { // truncated
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
          double tau = 0;
          if (t0) {
            tau = t0;
            if (st0) tau += st0*oneuni();
          }
          q[i] = tau + rwiener_diag2(R-1, bound-t0, a, vs, ws, err, K, epsFLAG);
          resp[i] = R;
        }
      });
    }

    /* the main thread also runs */
    int last = NperThread * (AmntOfThreads-1);
    double p_up, p_lo, vs = v, ws = w;
    bool REPEAT;
    for (int i = last; i < N; i++) {
      // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
      if(sv_or_sw) {
        REPEAT = true;
        while (REPEAT) {
          vs = v; ws = w;
          if (sv) vs += sv * onenorm();
          if (sw) ws += sw * (oneuni()-0.5);
          if (truncated) { // truncated
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
      double tau = 0;
      if (t0) {
        tau = t0;
        if (st0) tau += st0*oneuni();
      }
      q[i] = tau + rwiener_diag2(R-1, bound-t0, a, vs, ws, err, K, epsFLAG);
      resp[i] = R;
    }

    /* join threads */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } /* END "WITH PARALLELIZATION" */

  /* START "WITH _NO_ PARALLELIZATION" */
  else {

    double p_up, p_lo, vs = v, ws = w;
    bool REPEAT;
    for (int i = 0; i != N; i++) {
			if (i % 1024 == 0) R_CheckUserInterrupt();
      if(sv_or_sw) {
        REPEAT = true;
        while (REPEAT) {
          vs = v; ws = w;
          if (sv) vs += sv * onenorm();
          if (sw) ws += sw * (oneuni()-0.5);
          if (truncated) { // truncated
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
      double tau = 0;
      if (t0) {
        tau = t0;
        if (st0) tau += st0*oneuni();
      }
      q[i] = tau + rwiener_diag2(R-1, bound-t0, a, vs, ws, err, K, epsFLAG);
      resp[i] = R;
    }

  } /* END "WITH _NO_ PARALLELIZATION" */

}

void method3_both(int N, double a, double v, double w, double t0, double sv, double sw, double st0, int R, double bound, double err, int K, int epsFLAG, int NThreads, double *q, int *resp) {

  bool truncated = R_FINITE(bound), sv_or_sw = (sv > 0 || sw > 0);

  /* END "WITH PARALLELIZATION" */
  if (NThreads) {

    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* starting threads while ... */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        double p_up, p_lo, vs = v, ws = w, u;
        int up_or_down;
      	bool REPEAT;
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
          vs = v; ws = w;
          if (truncated) { // truncated
            if(sv_or_sw) {
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
          } else { // not truncated
            if (sv) vs += sv * onenorm();
        		if (sw) ws += sw * (oneuni()-0.5);
        		// p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
        		p_lo = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
        		up_or_down = oneuni() < p_lo ? 0 : 1;
          }
          double tau = 0;
          if (t0) {
            tau = t0;
            if (st0) tau += st0*oneuni();
          }
      		q[i] = tau + rwiener_diag2(up_or_down, bound-t0, a, vs, ws, err, K, epsFLAG);
      		resp[i] = up_or_down + 1;
        }
      });
    }

    /* the main thread also runs */
    int last = NperThread * (AmntOfThreads-1);
    double p_up, p_lo, vs = v, ws = w, u;
    int up_or_down;
  	bool REPEAT;
    for (int i = last; i < N; i++) {
      // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
      vs = v; ws = w;
      if (truncated) { // truncated
        if(sv_or_sw) {
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
      } else { // not truncated
        if (sv) {vs += sv * onenorm();}
    		if (sw) {ws += sw * (oneuni()-0.5);}
    		// p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
    		p_lo = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
    		up_or_down = oneuni() < p_lo ? 0 : 1;
      }
      double tau = 0;
      if (t0) {
        tau = t0;
        if (st0) tau += st0*oneuni();
      }
  		q[i] = tau + rwiener_diag2(up_or_down, bound-t0, a, vs, ws, err, K, epsFLAG);
  		resp[i] = up_or_down + 1;
    }

    /* join threads */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } /* END "WITH PARALLELIZATION" */

  /* START "WITH _NO_ PARALLELIZATION" */
  else {

    double p_up, p_lo, vs = v, ws = w, u;
    int up_or_down;
  	bool REPEAT;
  	for (int i = 0; i != N; i++) {
			if (i % 1024 == 0) R_CheckUserInterrupt();
  		vs = v; ws = w;
      if (truncated) { // truncated
        if(sv_or_sw) {
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
      } else { // not truncated
        if (sv) {vs += sv * onenorm();}
    		if (sw) {ws += sw * (oneuni()-0.5);}
    		// p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
    		p_lo = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
    		up_or_down = oneuni() < p_lo ? 0 : 1;
      }
      double tau = 0;
      if (t0) {
        tau = t0;
        if (st0) tau += st0*oneuni();
      }
  		q[i] = tau + rwiener_diag2(up_or_down, bound-t0, a, vs, ws, err, K, epsFLAG);
  		resp[i] = up_or_down + 1;
  	}

  } /* END "WITH _NO_ PARALLELIZATION" */

}


void method4_one(int N, double a, double v, double w, double t0, double sv, double sw, double st0, int R, double bound, double err, int K, int epsFLAG, int NThreads, double *q, int *resp) {

	bool truncated = R_FINITE(bound), sv_or_sw = (sv > 0 || sw > 0);

  /* END "WITH PARALLELIZATION" */
  if (NThreads) {

    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* starting threads while ... */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        double P, vs = v, ws = w;
        bool REPEAT;
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
          vs = v; ws = w;
          if(sv_or_sw) {
            REPEAT = true;
            while (REPEAT) {
              vs = v, ws = w;
              if (sv) vs += sv * onenorm();
              if (sw) ws += sw * (oneuni()-0.5);
              if (R == 2) {vs = -vs; ws = 1 - ws;}
              if (truncated) { // truncated
                P = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
              } else { // not truncated
                P = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
              }
              REPEAT = oneuni() > P;
            }
          } else {
            if (R == 2) {
              vs = -vs;
              ws = 1 - ws;
            }
          }
					ars_archiv ars_store_t;
					initialize_ars(a, vs, ws, 0.0, 0.0, bound-t0, ars_store_t);
					double tau = 0;
					if (t0) {
					  tau = t0;
					  if (st0) tau += st0*oneuni();
					}
					q[i] = tau + make_rwiener2(ars_store_t, nullptr, bound-t0, a, vs, ws, 0.0, 0.0, err, K, epsFLAG, 0);
          resp[i] = R;
        }
      });
    }

    /* the main thread also runs */
    int last = NperThread * (AmntOfThreads-1);
    double P, vs = v, ws = w;
    bool REPEAT;
		for (int i = last; i < N; i++) {
		  // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
      vs = v; ws = w;
      if(sv_or_sw) {
        REPEAT = true;
        while (REPEAT) {
          vs = v, ws = w;
          if (sv) vs += sv * onenorm();
          if (sw) ws += sw * (oneuni()-0.5);
          if (R == 2) {vs = -vs; ws = 1 - ws;}
          if (truncated) { // truncated
            P = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
          } else { // not truncated
            P = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
          }
          REPEAT = oneuni() > P;
        }
      } else {
        if (R == 2) {
          vs = -vs;
          ws = 1 - ws;
        }
      }
			ars_archiv ars_store1;
			initialize_ars(a, vs, ws, 0.0, 0.0, bound-t0, ars_store1);
			double tau = 0;
			if (t0) {
			  tau = t0;
			  if (st0) tau += st0*oneuni();
			}
			q[i] = tau + make_rwiener2(ars_store1, nullptr, bound-t0, a, vs, ws, 0.0, 0.0, err, K, epsFLAG, 0);
      resp[i] = R;
    }

    /* join threads */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } /* END "WITH PARALLELIZATION" */

  /* START "WITH _NO_ PARALLELIZATION" */
  else {

		double P, vs = v, ws = w;
    bool REPEAT;
    for (int i = 0; i != N; i++) {
			if (i % 1024 == 0) R_CheckUserInterrupt();
      vs = v; ws = w;
      if(sv_or_sw) {
        REPEAT = true;
        while (REPEAT) {
          vs = v, ws = w;
          if (sv) vs += sv * onenorm();
          if (sw) ws += sw * (oneuni()-0.5);
          if (R == 2) {vs = -vs; ws = 1 - ws;}
          if (truncated) { // truncated
            P = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
          } else { // not truncated
            P = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
          }
          REPEAT = oneuni() > P;
        }
      } else {
        if (R == 2) {
          vs = -vs;
          ws = 1 - ws;
        }
      }
			ars_archiv ars_store1;
			initialize_ars(a, vs, ws, 0.0, 0.0, bound-t0, ars_store1);
			double tau = 0;
			if (t0) {
			  tau = t0;
			  if (st0) tau += st0*oneuni();
			}
			q[i] = tau + make_rwiener2(ars_store1, nullptr, bound-t0, a, vs, ws, 0.0, 0.0, err, K, epsFLAG, 0);
      resp[i] = R;
    }

  } /* END "WITH _NO_ PARALLELIZATION" */

}

void method4_both(int N, double a, double v, double w, double t0, double sv, double sw, double st0, int R, double bound, double err, int K, int epsFLAG, int NThreads, double *q, int *resp) {

	bool truncated = R_FINITE(bound), sv_or_sw = (sv > 0 || sw > 0);

  /* END "WITH PARALLELIZATION" */
  if (NThreads) {

    /* prepare threads */
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(AmntOfThreads-1);

    /* starting threads while ... */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j] = std::thread([=]() {
        double p_up, p_lo, vs = v, ws = w, u;
        int up_or_down;
      	bool REPEAT;
        for (int i = j*NperThread; i < (j+1)*NperThread; i++) {
          // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
          vs = v; ws = w;
          if (truncated) { // truncated
            if(sv_or_sw) {
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
          } else { // not truncated
            if (sv) vs += sv * onenorm();
        		if (sw) ws += sw * (oneuni()-0.5);
        		// p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
        		p_lo = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
        		up_or_down = oneuni() < p_lo ? 0 : 1;
          }
					ars_archiv ars_store_t;
					if (up_or_down == 1) {
						initialize_ars(a, -vs, 1-ws, 0.0, 0.0, bound-t0, ars_store_t);
					  double tau = 0;
					  if (t0) {
					    tau = t0;
					    if (st0) tau += st0*oneuni();
					  }
						q[i] = tau + make_rwiener2(ars_store_t, nullptr, bound-t0, a, -vs, 1-ws, 0.0, 0.0, err, K, epsFLAG, 0);
					}
					else {
						initialize_ars(a, vs, ws, 0.0, 0.0, bound-t0, ars_store_t);
					  double tau = 0;
					  if (t0) {
					    tau = t0;
					    if (st0) tau += st0*oneuni();
					  }
						q[i] = tau + make_rwiener2(ars_store_t, nullptr, bound-t0, a, vs, ws, 0.0, 0.0, err, K, epsFLAG, 0);
					}
      		resp[i] = up_or_down + 1;
        }
      });
    }

    /* the main thread also runs */
    int last = NperThread * (AmntOfThreads-1);
    double p_up, p_lo, vs = v, ws = w, u;
    int up_or_down;
  	bool REPEAT;
    for (int i = last; i < N; i++) {
      // if (i % 1024 == 0) R_CheckUserInterruptGuarded();
      vs = v; ws = w;
      if (truncated) { // truncated
        if(sv_or_sw) {
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
      } else { // not truncated
        if (sv) {vs += sv * onenorm();}
    		if (sw) {ws += sw * (oneuni()-0.5);}
    		// p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
    		p_lo = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
    		up_or_down = oneuni() < p_lo ? 0 : 1;
      }
			ars_archiv ars_store_t;
			if (up_or_down == 1) {
				initialize_ars(a, -vs, 1-ws, 0.0, 0.0, bound-t0, ars_store_t);
			  double tau = 0;
			  if (t0) {
			    tau = t0;
			    if (st0) tau += st0*oneuni();
			  }
				q[i] = tau + make_rwiener2(ars_store_t, nullptr, bound-t0, a, -vs, 1-ws, 0.0, 0.0, err, K, epsFLAG, 0);
			}
			else {
				initialize_ars(a, vs, ws, 0.0, 0.0, bound-t0, ars_store_t);
			  double tau = 0;
			  if (t0) {
			    tau = t0;
			    if (st0) tau += st0*oneuni();
			  }
				q[i] = tau + make_rwiener2(ars_store_t, nullptr, bound-t0, a, vs, ws, 0.0, 0.0, err, K, epsFLAG, 0);
			}
			resp[i] = up_or_down + 1;
    }

    /* join threads */
    for (int j = 0; j < AmntOfThreads-1; j++) {
      threads[j].join();
    }

  } /* END "WITH PARALLELIZATION" */

  /* START "WITH _NO_ PARALLELIZATION" */
  else {

    double p_up, p_lo, vs = v, ws = w, u;
    int up_or_down;
  	bool REPEAT;
  	for (int i = 0; i != N; i++) {
			if (i % 1024 == 0) R_CheckUserInterrupt();
  		vs = v; ws = w;
      if (truncated) { // truncated
        if(sv_or_sw) {
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
      } else { // not truncated
        if (sv) {vs += sv * onenorm();}
    		if (sw) {ws += sw * (oneuni()-0.5);}
    		// p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
    		p_lo = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
    		up_or_down = oneuni() < p_lo ? 0 : 1;
      }
			ars_archiv ars_store_t;
			if (up_or_down == 1) {
				initialize_ars(a, -vs, 1-ws, 0.0, 0.0, bound-t0, ars_store_t);
			  double tau = 0;
			  if (t0) {
			    tau = t0;
			    if (st0) tau += st0*oneuni();
			  }
				q[i] = tau + make_rwiener2(ars_store_t, nullptr, bound-t0, a, -vs, 1-ws, 0.0, 0.0, err, K, epsFLAG, 0);
			}
			else {
				initialize_ars(a, vs, ws, 0.0, 0.0, bound-t0, ars_store_t);
			  double tau = 0;
			  if (t0) {
			    tau = t0;
			    if (st0) tau += st0*oneuni();
			  }
				q[i] = tau + make_rwiener2(ars_store_t, nullptr, bound-t0, a, vs, ws, 0.0, 0.0, err, K, epsFLAG, 0);
			}
			resp[i] = up_or_down + 1;
  	}

  } /* END "WITH _NO_ PARALLELIZATION" */

}


// R=0 is lower bound and R=1 is upper bound
void run_make_rwiener(int choice, int N, double a, double v, double w, double t0, double sv, double sw, double st0, int R, double bound, double err, int K, int epsFLAG, int NThreads, double *q, int *resp, ars_archiv *ars_store1, ars_archiv *ars_store2, int use_store) {
	// double vs, ws;
	// printf("h1 size = %d\n", static_cast<int>(ars_store1->hstore.size()));
	// printf("h2 size = %d\n", static_cast<int>(ars_store2->hstore.size()));
	switch (choice) {
		case 1:
			if (R != 0) { // one-sided

				// // ars_archiv ars_store;
				// if(R == 2) {
				// 	v = -v;
				// 	w = 1-w;
				// }
				// // initialize_ars(a, v, w, sw, sv, bound, ars_store);
				// if (!use_store) initialize_ars(a, v, w, sw, sv, bound, *ars_store1);
				// for (int i = 0; i != N; i++) {
				// 	// q[i] = make_rwiener2(ars_store, bound, a, v, w, sw, sv, err, K, epsFLAG);
				// 	q[i] = make_rwiener2(*ars_store1, bound, a, v, w, sw, sv, err, K, epsFLAG);
				// 	resp[i] = R;
				// }
				method1_one(N, a, v, w, t0, sv, sw, st0, R, bound, err, K, epsFLAG, NThreads, q, resp, ars_store1, ars_store2, use_store);

			} else { // R = 0 -- both sides

				// // ars_archiv ars_store1;
				// // ars_archiv ars_store2;
				// double p_up; //, p_lo;
				// if (std::isfinite(bound)) { // truncated
				// 	if (sv || sw) {
				// 		double temp_u, temp_l;
	      //     pdiff(0, bound, 1.0, a, v, 0, w, sw, sv, 0, err, K, epsFLAG, &temp_u);
				// 		pdiff(0, bound, -1.0, a, v, 0, w, sw, sv, 0, err, K, epsFLAG, &temp_l);
				// 		p_up = temp_u / (temp_u + temp_l);
				// 		// p_lo = 1-p_up;
				// 	} else {
	      //     double temp_u = exp(pwiener(bound, a, -v, 1-w, err, K, epsFLAG));
				// 		double temp_l = exp(pwiener(bound, a, v, w, err, K, epsFLAG));
				// 		p_up = temp_u / (temp_u + temp_l);
				// 		// p_lo = 1-p_up;
				// 	}
				// } else { // not truncated
				// 	if (sv || sw) {
	      //     pdiff(0, bound, 1.0, a, v, 0, w, sw, sv, 0, err, K, epsFLAG, &p_up);
				// 		// p_lo = 1-p_up;
				// 	} else {
	      //     // p_up = exp(pwiener(bound, a, -v, 1-w, err, K, epsFLAG));
				// 		p_up = (1-exp(2*v*a*w))/(exp(-2*v*a*(1-w))-exp(2*v*a*w));
				// 		// p_lo = 1-p_up;
				// 	}
				// }
				// int cnt_up = 0;
				// for (int i = 0; i !=N; i++) {
				// 	cnt_up += oneuni()<=p_up ? 1 : 0;
				// }
				//
				// // printf("ars h1 size = %d\n", static_cast<int>(ars_store1->hstore.size()));
				// if (!use_store) initialize_ars(a, -v, 1-w, sw, sv, bound, *ars_store1);
				// for (int i = 0; i != cnt_up; i++) {
				// 	q[i] = make_rwiener2(*ars_store1, ars_store1, bound, a, -v, 1-w, sw, sv, err, K, epsFLAG, 0);
				// 	resp[i] = 2;
				// }
				//
				// // printf("ars h2 size = %d\n", static_cast<int>(ars_store2->hstore.size()));
				// if (!use_store) initialize_ars(a, v, w, sw, sv, bound, *ars_store2);
				// for (int i = cnt_up; i != N; i++) {
				// 	q[i] = make_rwiener2(*ars_store2, ars_store2, bound, a, v, w, sw, sv, err, K, epsFLAG, 0);
				// 	resp[i] = 1;
				// }
				method1_both(N, a, v, w, t0, sv, sw, st0, R, bound, err, K, epsFLAG, NThreads, q, resp, ars_store1, ars_store2, use_store);

			}
			break;

		case 2:
			if(R != 0) { // one-sided
				// double P;
				// bool REPEAT;
				// for (int i = 0; i != N; i++) {
				// 	vs = v; ws = w;
				// 	if(sv || sw) {
				// 		REPEAT = true;
				// 		while (REPEAT) {
				// 			vs = v, ws = w;
				// 			if (sv) vs += sv * onenorm();
				// 			if (sw) ws += sw * (oneuni()-0.5);
				// 			if (R == 2) {vs = -vs; ws = 1 - ws;}
				// 			if (std::isfinite(bound)) P = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
				// 			if (std::isinf(bound)) P = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
				// 			REPEAT = oneuni() > P;
				// 		}
				// 	} else {
				// 		if (R == 2) {
				// 			vs = -vs;
				// 			ws = 1 - ws;
				// 		}
				// 	}
				// 	q[i] = -rdiffusion_lower_trunc(bound, a, vs, ws);
				// 	resp[i] = R;
				// }
        method2_one(N, a, v, w, t0, sv, sw, st0, R, bound, err, K, epsFLAG, NThreads, q, resp);

			} else { // R = 0 -- both sides
			// 	if (std::isfinite(bound)) { // truncated
			// 		double p_up, p_lo;
			// 		bool REPEAT;
			// 		for (int i = 0; i != N; i++) {
			// 			vs = v; ws = w;
			// 			if(sv || sw) {
			// 				REPEAT = true;
			// 				while (REPEAT) {
			// 					vs = v; ws = w;
			// 					if (sv) vs += sv*onenorm();
			// 					if (sw) ws += sw*(oneuni()-0.5);
			// 					p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
			// 					p_up = exp(pwiener(bound, a, -vs, 1-ws, err, K, epsFLAG));
			// 					REPEAT = oneuni() > (p_up + p_lo);
			// 				}
			// 			}
			// 			q[i] = rdiffusion_UPbound(bound, a, vs, ws);
			// 			resp[i] = q[i] > 0 ? 2 : 1;
			// 			if (resp[i] == 1) q[i] = fabs(q[i]);
			// 		}
			// 	} else { // not truncated
			// 		for (int i = 0; i != N; i++) {
			// 			vs = v; ws = w;
			// 			if (sv) vs += sv*onenorm();
			// 			if (sw) ws += sw*(oneuni()-0.5);
			// 			q[i] = rdiffusion_UPbound(bound, a, vs, ws);
			// 			resp[i] = q[i] > 0 ? 2 : 1;
			// 			if (resp[i] == 1) q[i] = fabs(q[i]);
			// 		}
			//  }
      method2_both(N, a, v, w, t0, sv, sw, st0, R, bound, err, K, epsFLAG, NThreads, q, resp);
			}
			break;

		case 3:
			// double p_lo, p_up;
			if (R != 0) { // one-sided
				// vs = v; ws = w;
				// bool REPEAT;
				// for (int i = 0; i != N; i++) {
				// 	if(sv || sw) {
				// 		REPEAT = true;
				// 		while (REPEAT) {
				// 			vs = v; ws = w;
				// 			if (sv) vs += sv * onenorm();
				// 			if (sw) ws += sw * (oneuni()-0.5);
				// 			if (std::isfinite(bound)) { // truncated
				// 				if (R == 1) {
				// 					p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
				// 					REPEAT = oneuni() > p_lo;
				// 				}
				// 				if (R == 2) {
				// 					p_up = exp(pwiener(bound, a, -vs, 1-ws, err, K, epsFLAG));
				// 					REPEAT = oneuni() > p_up;
				// 				}
				// 			} else { // not truncated
				// 				// p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
				// 				p_lo = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
				// 				if (R == 1) REPEAT = oneuni() > p_lo;
				// 				if (R == 2) REPEAT = oneuni() > 1-p_lo;
				// 			}
				// 		}
				// 	}
				// 	q[i] = rwiener_diag2(R-1, bound, a, vs, ws, err, K, epsFLAG);
				// 	resp[i] = R;
				// }
        method3_one(N, a, v, w, t0, sv, sw, st0, R, bound, err, K, epsFLAG, NThreads, q, resp);

			} else { // R = 0 -- both sides

				// int up_or_down;
				// if (std::isfinite(bound)) { // truncated
				// 	bool REPEAT;
				// 	double u;
				// 	for (int i = 0; i != N; i++) {
				// 		vs = v; ws = w;
				// 		if(sv || sw) {
				// 			REPEAT = true;
				// 			while (REPEAT) {
				// 				vs = v; ws = w;
				// 				if (sv) vs += sv * onenorm();
				// 				if (sw) ws += sw * (oneuni()-0.5);
				// 				p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
				// 				p_up = exp(pwiener(bound, a, -vs, 1-ws, err, K, epsFLAG));
				// 				u = oneuni();
				// 				if (u <= p_lo) {
				// 					REPEAT = false;
				// 					up_or_down = 0;
				// 				} else if (u >= 1-p_up) {
				// 					REPEAT = false;
				// 					up_or_down = 1;
				// 				} else {
				// 					REPEAT = true;
				// 				}
				// 			}
				// 		} else {
				// 			p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
				// 			p_up = exp(pwiener(bound, a, -vs, 1-ws, err, K, epsFLAG));
				// 			up_or_down = oneuni() <= p_up / (p_up + p_lo) ? 1 : 0;
				// 		}
				// 		q[i] = rwiener_diag2(up_or_down, bound, a, vs, ws, err, K, epsFLAG);
				// 		resp[i] = up_or_down + 1;
				// 	}
				// } else { // not truncated
				// 	for (int i = 0; i != N; i++) {
				// 		vs = v; ws = w;
				// 		if (sv) vs += sv * onenorm();
				// 		if (sw) ws += sw * (oneuni()-0.5);
				// 		// p_lo = exp(pwiener(bound, a, vs, ws, err, K, epsFLAG));
				// 		p_lo = (1-exp(-2*vs*a*(1-ws)))/(exp(2*vs*a*ws)-exp(-2*vs*a*(1-ws)));
				// 		up_or_down = oneuni() < p_lo ? 0 : 1;
				// 		q[i] = rwiener_diag2(up_or_down, bound, a, vs, ws, err, K, epsFLAG);
				// 		resp[i] = up_or_down + 1;
				// 	}
				// }
        method3_both(N, a, v, w, t0, sv, sw, st0, R, bound, err, K, epsFLAG, NThreads, q, resp);

			}
			break;

		case 4:

			if (R != 0) { // one-sided

				method4_one(N, a, v, w, t0, sv, sw, st0, R, bound, err, K, epsFLAG, NThreads, q, resp);

			} else { // R = 0 -- both sides

				method4_both(N, a, v, w, t0, sv, sw, st0, R, bound, err, K, epsFLAG, NThreads, q, resp);

			}
			break;

	}
}
