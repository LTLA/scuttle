#include "Rtatami.h"

#include "tatami_stats/tatami_stats.hpp"
#include "scuttle/downsample.h"

#include <vector>
#include <cmath>
#include <iostream>

//[[Rcpp::export]]
Rcpp::RObject downsample(Rcpp::RObject input, double prop_global, Rcpp::Nullable<Rcpp::NumericVector> prop_column, int num_threads) {
    Rtatami::BoundNumericPointer ptr(input);

    // Truncating and replacing negative values with zero.
    tatami::DelayedUnaryIsometricOperation<double, double, int> mat(
        std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
            ptr->ptr,
            std::make_shared<tatami::DelayedUnaryIsometricRoundHelper<double, double, int> >()
        ),
        std::make_shared<tatami::DelayedUnaryIsometricSubstituteLessThanScalarHelper<double, double, int, double> >(0., 0.)
    );

    tatami_stats::sums::Options opt;
    opt.num_threads = num_threads;
    auto colsums = tatami_stats::sums::by_column(mat, opt);

    const int ngenes = mat.nrow();
    const int ncells = mat.ncol();

    std::vector<double> partials;
    double required = 0;
    Rcpp::NumericVector prop_col;
    if (!prop_column.isNull()) {
        prop_col = Rcpp::NumericVector(prop_column);
        if (prop_col.size() != ncells) {
            throw std::runtime_error("length of 'prop' should equal number of columns");
        }
        for (auto x : colsums) {
            if (scuttle::too_large_for_integer_precision(x)) {
                throw std::runtime_error("column sums are too large to maintain integer precision");
            }
        }

    } else {
        partials.resize(colsums.size());
        double last = 0;
        for (int c = 0; c < ncells; ++c) {
            auto idx = ncells - c - 1;
            partials[idx] = last + colsums[idx];
            last = partials[idx];
        }

        if (scuttle::too_large_for_integer_precision(last)) {
            throw std::runtime_error("matrix sum is too large to maintain integer precision");
        }
        required = std::round(last * prop_global);
    }

    Rcpp::List output(ncells);
    std::vector<double> down_x;
    down_x.reserve(ngenes);
    std::vector<int> down_i;
    down_i.reserve(ngenes);

    // Can't easily parallelize as we're using a global RNG... oh well.
    // I suppose we could do so by switching to PCG32; for global sampling,
    // we can hierarchically sample the number in each column before sampling within each column.
    if (mat.is_sparse()) {
        std::vector<int> work_i(ngenes);
        std::vector<double> work_x(ngenes);
        auto ext = tatami::consecutive_extractor<true>(mat, false, 0, ncells);

        for (int c = 0; c < ncells; ++c) {
            double current_total = (prop_column.isNull() ? partials[c] : colsums[c]);
            double current_required = (prop_column.isNull() ? required : std::round(current_total * prop_col[c]));

            auto range = ext->fetch(work_x.data(), work_i.data());
            scuttle::downsample(range.value, range.value + range.number, work_x.data(), current_total, current_required);

            down_i.clear();
            down_x.clear();
            for (int i = 0; i < range.number; ++i) {
                if (work_x[i] > 0) {
                    down_i.push_back(range.index[i]);
                    down_x.push_back(work_x[i]);
                }
            }

            output[c] = Rcpp::List::create(
                Rcpp::NumericVector(down_x.begin(), down_x.end()),
                Rcpp::IntegerVector(down_i.begin(), down_i.end())
            );

            if (prop_column.isNull()) {
                required = current_required;
            }
        }

    } else {
        std::vector<double> work_x(ngenes);
        auto ext = tatami::consecutive_extractor<false>(mat, false, 0, ncells);

        for (int c = 0; c < ncells; ++c) {
            double current_total = (prop_column.isNull() ? partials[c] : colsums[c]);
            double current_required = (prop_column.isNull() ? required : std::round(current_total * prop_col[c]));

            const auto ptr = ext->fetch(work_x.data());
            scuttle::downsample(ptr, ptr + ngenes, work_x.data(), current_total, current_required);

            down_i.clear();
            down_x.clear();
            for (int i = 0; i < ngenes; ++i) {
                if (work_x[i] > 0) {
                    down_i.push_back(i);
                    down_x.push_back(work_x[i]);
                }
            }

            output[c] = Rcpp::List::create(
                Rcpp::NumericVector(down_x.begin(), down_x.end()),
                Rcpp::IntegerVector(down_i.begin(), down_i.end())
            );

            if (prop_column.isNull()) {
                required = current_required;
            }
        }
    }

    return output;
}
