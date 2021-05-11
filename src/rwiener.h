


#ifndef RWIENER_H
#define RWIENER_H

#include "tools.h"
#include "cdf_fncs.h"
#include "pdf_fncs.h"
#include "adaptive_rejection_samp.h"
#include "fncs_seven.h"
#include "cubature.h"
#include <vector>
#include <Rinternals.h>
#include <algorithm>    // std::sort


void run_make_rwiener(int, int, double, double, double, double, double, int, double, double, int, int, double*, int*);


#endif
