#include "Rcpp.h"

template<class Functor>
Rcpp::NumericMatrix sparse_aggregate(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, Rcpp::IntegerVector groupings, int ngroups, int nrows, Functor fun) {
    Rcpp::NumericMatrix output(nrows, ngroups);

    for (size_t c = 0, cend = groupings.size(); c < cend; ++c) {
        auto g = groupings[c];
        if (g == NA_INTEGER) {
            continue;
        }

        int start = p[c], end = p[c + 1];
        auto col = output.column(g);
        for (int idx = start; idx < end; ++idx) {
            fun(col[i[idx]],  x[idx]);
        }
    }

    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix sparse_aggregate_sum(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, Rcpp::IntegerVector groupings, int ngroups, int nrows) {
    return sparse_aggregate(x, i, p, groupings, ngroups, nrows, [](double& y, double incoming) -> void { y += incoming; });
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix sparse_aggregate_detected(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, Rcpp::IntegerVector groupings, int ngroups, int nrows) {
    return sparse_aggregate(x, i, p, groupings, ngroups, nrows, [](double& y, double incoming) -> void { y += (incoming > 0); });
}
