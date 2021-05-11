
#ifndef ADAPTIVE_REJECTION_SAMP_H
#define ADAPTIVE_REJECTION_SAMP_H

#include "tools.h"
#include <Rinternals.h>

double fun_upper(int, double, std::vector<piece>);
void generate_intervals(int&, double, std::vector<point>, std::vector<piece>&, std::vector<piece>&, std::vector<double>&);
bool update_intervals(int, double, point, std::vector<point>&, std::vector<piece>&, std::vector<piece>&, std::vector<double>&);
double fun_lower(int, double, std::vector<point>, std::vector<piece>);
double inverse_distribution(int, double, std::vector<piece>, std::vector<double>, double, bool&);



#endif
