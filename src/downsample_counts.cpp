#include "Rtatami.h"

#include "tatami_stats/tatami_stats.hpp"

#include <vector>
#include <cmath>

//[[Rcpp::export]]
Rcpp::RObject downsample(Rcpp::RObject input, double prop_global, Rcpp::Nullable<Rcpp::NumericVector> prop_column, int num_threads) {
    Rtatami::BoundNumericPointer ptr(input);
    const auto& mat = *(ptr->ptr);
    const int ngenes = mat.nrow();
    const int ncells = mat.ncol();

    tatami_stats::sums::Options opt;
    opt.num_threads = num_threads;
    auto colsums = tatami_stats::sums::by_column(mat, opt);

    // Computing cumulative sums to mitigate problems from loss of precision at very large counts.
    std::vector<double> partials;
    double required = 0;
    Rcpp::NumericVector prop_col;
    if (!prop_column.isNull()) {
        prop_col = Rcpp::NumericVector(prop_column);
        if (prop_col.size() != ncells) {
            throw std::runtime_error("length of 'prop' should equal number of columns");
        }
    } else {
        partials.resize(colsums.size());
        double last = 0;
        for (int c = 0; c < ncells; ++c) {
            auto idx = ncells - c - 1;
            partials[idx] = last + colsums[idx];
            last = partials[idx];
        }
        required = last * prop_global;
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
            auto range = ext->fetch(work_x.data(), work_i.data());
            down_i.clear();
            down_x.clear();

            const double current_total = (prop_column.isNull() ? partials[c] : colsums[c]);
            const double current_required = (prop_column.isNull() ? required : current_total * prop_col[c]);
            double accumulated_total = 0, accumulated_used = 0;

            for (int i = 0; i < range.number; ++i) {
                accumulated_total += range.value[i];
                if (accumulated_total >= current_total || accumulated_used >= current_required) {
                    break;
                }

                const auto out = Rf_rhyper(range.value[i], current_total - accumulated_total, current_required - accumulated_used);
                accumulated_used += out;
                if (out > 0) {
                    down_i.push_back(range.index[i]);
                    down_x.push_back(out);
                }
            }

            output[c] = Rcpp::List::create(
                Rcpp::NumericVector(down_x.begin(), down_x.end()),
                Rcpp::IntegerVector(down_i.begin(), down_i.end())
            );

            if (prop_column.isNull()) {
                required -= accumulated_used;
            }
        }

    } else {
        std::vector<double> work_x(ngenes);
        auto ext = tatami::consecutive_extractor<false>(mat, false, 0, ncells);

        for (int c = 0; c < ncells; ++c) {
            const auto ptr = ext->fetch(work_x.data());
            down_i.clear();
            down_x.clear();

            const double current_total = (prop_column.isNull() ? partials[c] : colsums[c]);
            const double current_required = (prop_column.isNull() ? required : current_total * prop_col[c]);
            double accumulated_total = 0, accumulated_used = 0;

            for (int i = 0; i < ngenes; ++i) {
                accumulated_total += ptr[i];
                if (accumulated_total >= current_total || accumulated_used >= current_required) {
                    break;
                }

                const auto out = Rf_rhyper(ptr[i], current_total - accumulated_total, current_required - accumulated_used);
                accumulated_used += out;
                if (out > 0) {
                    down_i.push_back(i);
                    down_x.push_back(out);
                }
            }

            output[c] = Rcpp::List::create(
                Rcpp::NumericVector(down_x.begin(), down_x.end()),
                Rcpp::IntegerVector(down_i.begin(), down_i.end())
            );

            if (prop_column.isNull()) {
                required -= accumulated_used;
            }
        }
    }

    return output;
}
