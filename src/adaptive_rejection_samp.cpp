//implements adaptive rejection sampling

#include "adaptive_rejection_samp.h"


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
		t -= log(fabs(upper[i].slope));

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
			t -= log(fabs(upper[j].slope));
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
		t -= log(fabs(sl));
		s[k - 1] = logsum(t, s[k - 2]);
		sum = s[k - 1];
	}
	int j = 0;
	double temp = log(xstar) + sum;
	while (temp > s[j]) j++;
	if (j > k - 1) { Rprintf("Wie das?\n"); }


	double sl = upper[j].slope;
	double help = log(fabs(sl)); int sign = sl > 0 ? 1 : -1;
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
		Rprintf("\nnanu j=%d; k-1=%d; t=%g; upper[j]=%g; upper[j+1]=%g; s[j-1]=%g; upper slope=%g; upper absc=%g; temp=%g; fun_upper[j]=%g; fun_upper[j+1]%g\n",
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
