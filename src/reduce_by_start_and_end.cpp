#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
long reduce_by_start_and_end(IntegerVector s, IntegerVector e) {
	int n = s.size();

	if(n == 1) {
		return e[0] - s[0] + 1;
	}

	long w = 0;
	int prev_s = s[0];
	int prev_e = e[0];
	int updated = 0;

	for(int i = 1; i < n; i ++ ){
		if(s[i] - prev_e > 1) {
			w += prev_e - prev_s + 1;
			prev_s = s[i];
			prev_e = e[i];
			updated = 1;
		} else if(e[i] > prev_e) {
			prev_e = e[i];
			updated = 1;
		} else {
			updated = 0;
		}

		if(i == n - 1 & updated) {
			w += prev_e - prev_s + 1;
		}
	}

	return w;
}