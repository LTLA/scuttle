#include "Rtatami.h"
#include "scran_aggregate/scran_aggregate.hpp"

//[[Rcpp::export(rng=false)]]
SEXP aggregate_across_cells(SEXP x, Rcpp::IntegerVector groups, int num_groups, bool do_sum, bool do_detected, int num_threads) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto NC = mat->ncol();
    const auto NR = mat->nrow();

    if (!sanisizer::is_equal(groups.size(), NC)) {
        throw std::runtime_error("length of 'groups' should be equal to the number of columns in 'x'");
    }

    const int* gptr = groups.begin();

    Rcpp::List output(2);
    scran_aggregate::AggregateAcrossCellsBuffers<double, int> buffers;

    if (do_sum) {
        Rcpp::NumericMatrix sums(NR, num_groups);
        buffers.sums.reserve(num_groups);
        double* osum = sums.begin();
        for (int i = 0; i < num_groups; ++i) {
            buffers.sums.push_back(osum + sanisizer::product_unsafe<std::size_t>(NR, i));
        }
        output[0] = sums;
    } else {
        output[0] = R_NilValue;
    }

    if (do_detected) {
        Rcpp::IntegerMatrix detected(NR, num_groups);
        buffers.detected.reserve(num_groups);
        int* odetected = detected.begin();
        for (int i = 0; i < num_groups; ++i) {
            buffers.detected.push_back(odetected + sanisizer::product_unsafe<std::size_t>(NR, i));
        }
        output[1] = detected;
    } else {
        output[1] = R_NilValue;
    }

    scran_aggregate::AggregateAcrossCellsOptions opt;
    opt.num_threads = num_threads;
    scran_aggregate::aggregate_across_cells(*mat, gptr, buffers, opt);

    return output;
}
