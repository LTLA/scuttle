#include "Rtatami.h"

#include "tatami_stats/tatami_stats.hpp"
#include "scuttle/downsample.h"

#include <vector>
#include <cmath>
#include <stdexcept>

static std::shared_ptr<tatami::Matrix<double, int> > sanitize_matrix(std::shared_ptr<tatami::Matrix<double, int> > ptr, bool already_integer) {
    std::shared_ptr<tatami::Matrix<double, int> > output;

    // Rounding values if we're not sure that they're already integer. 
    if (!already_integer) {
        output = std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
            std::move(ptr),
            std::make_shared<tatami::DelayedUnaryIsometricRoundHelper<double, double, int> >()
        );
    } else {
        output = std::move(ptr);
    }

    // Replacing negative values with zero.
    return std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
        std::move(output),
        std::make_shared<tatami::DelayedUnaryIsometricSubstituteLessThanScalarHelper<double, double, int, double> >(0., 0.)
    );
}

/**********************************/

template<typename OutputThing_, bool sparse_output_>
Rcpp::RObject downsample_global_internal(const tatami::Matrix<double, int>& mat, double total, double required) {
    const auto ngenes = mat.nrow();
    const auto ncells = mat.ncol();

    auto output = [&](){
        if constexpr(sparse_output_) {
            return Rcpp::List(ncells);
        } else {
            return OutputThing_(ngenes, ncells);
        }
    }();

    typename std::conditional<sparse_output_, std::vector<double>, bool>::type down_x;
    typename std::conditional<sparse_output_, std::vector<int>, bool>::type down_i;
    if constexpr(sparse_output_) {
        down_x.reserve(ngenes);
        down_i.reserve(ngenes);
    }

    // Can't easily parallelize as we're using a global RNG... oh well.
    // I suppose we could do so by switching to PCG32; for global sampling,
    // we can hierarchically sample the number in each column before sampling within each column.
    if (mat.is_sparse()) {
        std::vector<int> work_i(ngenes);
        std::vector<double> work_x(ngenes);
        auto ext = tatami::consecutive_extractor<true>(mat, false, 0, ncells);

        for (int c = 0; c < ncells; ++c) {
            auto range = ext->fetch(work_x.data(), work_i.data());
            scuttle::downsample(range.value, range.value + range.number, work_x.data(), total, required);

            if constexpr(sparse_output_) {
                down_i.clear();
                down_x.clear();
                for (int i = 0; i < range.number; ++i) {
                    if (work_x[i] > 0) {
                        down_i.push_back(range.index[i]);
                        down_x.push_back(work_x[i]);
                    }
                }

                output[c] = Rcpp::List::create(
                    OutputThing_(down_x.begin(), down_x.end()),
                    Rcpp::IntegerVector(down_i.begin(), down_i.end())
                );
            } else {
                for (int i = 0; i < range.number; ++i) {
                    output(range.index[i], c) = work_x[i];
                }
            }
        }

    } else {
        std::vector<double> work_x(ngenes);
        auto ext = tatami::consecutive_extractor<false>(mat, false, 0, ncells);

        for (int c = 0; c < ncells; ++c) {
            const auto ptr = ext->fetch(work_x.data());
            scuttle::downsample(ptr, ptr + ngenes, work_x.data(), total, required);

            if constexpr(sparse_output_) {
                down_i.clear();
                down_x.clear();
                for (int i = 0; i < ngenes; ++i) {
                    if (work_x[i] > 0) {
                        down_i.push_back(i);
                        down_x.push_back(work_x[i]);
                    }
                }

                output[c] = Rcpp::List::create(
                    OutputThing_(down_x.begin(), down_x.end()),
                    Rcpp::IntegerVector(down_i.begin(), down_i.end())
                );
            } else {
                auto col = output.column(c);
                std::copy(work_x.begin(), work_x.end(), col.begin());
            }
        }
    }

    return output;
}

//[[Rcpp::export]]
Rcpp::RObject downsample_global(Rcpp::RObject input, double prop_global, int num_threads, bool already_integer, bool output_sparse, bool output_integer) {
    Rtatami::BoundNumericPointer ptr(input);
    auto sanitized = sanitize_matrix(ptr->ptr, already_integer);
    const auto& mat = *sanitized;

    // Just computing the total sum in a decently fast manner.
    // Probably could replace with a sums::total() function to switch between row and column sums, but that would be pretty niche.
    tatami_stats::sums::Options opt;
    opt.num_threads = num_threads;
    auto colsums = tatami_stats::sums::by_column(mat, opt);

    const double total = std::accumulate(colsums.begin(), colsums.end(), 0.0);
    if (scuttle::too_large_for_integer_precision(total)) {
        throw std::runtime_error("matrix sum is too large to maintain integer precision");
    }
    if (prop_global < 0 || prop_global > 1) {
        throw std::runtime_error("'prop' should be in [0, 1]");
    }
    const double required = std::round(total * prop_global);

    if (output_sparse) {
        if (output_integer) {
            return downsample_global_internal<Rcpp::IntegerVector, true>(mat, total, required);
        } else {
            return downsample_global_internal<Rcpp::NumericVector, true>(mat, total, required);
        }
    } else {
        if (output_integer) {
            return downsample_global_internal<Rcpp::IntegerMatrix, false>(mat, total, required);
        } else {
            return downsample_global_internal<Rcpp::NumericMatrix, false>(mat, total, required);
        }
    }
}

/**********************************/

template<typename OutputThing_, bool sparse_output_>
Rcpp::RObject downsample_column_internal(const tatami::Matrix<double, int>& mat, Rcpp::NumericVector prop_column) {
    const auto ngenes = mat.nrow();
    const auto ncells = mat.ncol();

    auto output = [&](){
        if constexpr(sparse_output_) {
            return Rcpp::List(ncells);
        } else {
            return OutputThing_(ngenes, ncells);
        }
    }();

    typename std::conditional<sparse_output_, std::vector<double>, bool>::type down_x;
    typename std::conditional<sparse_output_, std::vector<int>, bool>::type down_i;
    if constexpr(sparse_output_) {
        down_x.reserve(ngenes);
        down_i.reserve(ngenes);
    }

    // Can't easily parallelize as we're using a global RNG... oh well.
    // I suppose we could do so by switching to PCG32; for global sampling,
    // we can hierarchically sample the number in each column before sampling within each column.
    if (mat.is_sparse()) {
        std::vector<int> work_i(ngenes);
        std::vector<double> work_x(ngenes);
        auto ext = tatami::consecutive_extractor<true>(mat, false, 0, ncells);

        for (int c = 0; c < ncells; ++c) {
            const auto range = ext->fetch(work_x.data(), work_i.data());
            double current_total = std::accumulate(range.value, range.value + range.number, 0.);
            if (scuttle::too_large_for_integer_precision(current_total)) {
                throw std::runtime_error("column sum is too large to maintain integer precision");
            }

            double current_required = std::round(current_total * prop_column[c]);
            scuttle::downsample(range.value, range.value + range.number, work_x.data(), current_total, current_required);

            if constexpr(sparse_output_) {
                down_i.clear();
                down_x.clear();
                for (int i = 0; i < range.number; ++i) {
                    if (work_x[i] > 0) {
                        down_i.push_back(range.index[i]);
                        down_x.push_back(work_x[i]);
                    }
                }

                output[c] = Rcpp::List::create(
                    OutputThing_(down_x.begin(), down_x.end()),
                    Rcpp::IntegerVector(down_i.begin(), down_i.end())
                );
            } else {
                for (int i = 0; i < range.number; ++i) {
                    output(range.index[i], c) = work_x[i];
                }
            }
        }

    } else {
        std::vector<double> work_x(ngenes);
        auto ext = tatami::consecutive_extractor<false>(mat, false, 0, ncells);

        for (int c = 0; c < ncells; ++c) {
            const auto ptr = ext->fetch(work_x.data());
            auto current_total = std::accumulate(ptr, ptr + ngenes, 0.);
            if (scuttle::too_large_for_integer_precision(current_total)) {
                throw std::runtime_error("column sum is too large to maintain integer precision");
            }

            double current_required = std::round(current_total * prop_column[c]);
            scuttle::downsample(ptr, ptr + ngenes, work_x.data(), current_total, current_required);

            if constexpr(sparse_output_) {
                down_i.clear();
                down_x.clear();
                for (int i = 0; i < ngenes; ++i) {
                    if (work_x[i] > 0) {
                        down_i.push_back(i);
                        down_x.push_back(work_x[i]);
                    }
                }

                output[c] = Rcpp::List::create(
                    OutputThing_(down_x.begin(), down_x.end()),
                    Rcpp::IntegerVector(down_i.begin(), down_i.end())
                );
            } else {
                auto col = output.column(c);
                std::copy(work_x.begin(), work_x.end(), col.begin());
            }
        }
    }

    return output;
}

//[[Rcpp::export]]
Rcpp::RObject downsample_column(Rcpp::RObject input, Rcpp::NumericVector prop_col, bool already_integer, bool output_sparse, bool output_integer) {
    Rtatami::BoundNumericPointer ptr(input);
    auto sanitized = sanitize_matrix(ptr->ptr, already_integer);
    const auto& mat = *sanitized;

    if (prop_col.size() != mat.ncol()) {
        throw std::runtime_error("length of 'prop' should be equal to the number of columns");
    }
    for (auto p : prop_col) {
        if (p < 0 || p > 1) {
            throw std::runtime_error("'prop' should be in [0, 1]");
        }
    }

    if (output_sparse) {
        if (output_integer) {
            return downsample_column_internal<Rcpp::IntegerVector, true>(mat, prop_col);
        } else {
            return downsample_column_internal<Rcpp::NumericVector, true>(mat, prop_col);
        }
    } else {
        if (output_integer) {
            return downsample_column_internal<Rcpp::IntegerMatrix, false>(mat, prop_col);
        } else {
            return downsample_column_internal<Rcpp::NumericMatrix, false>(mat, prop_col);
        }
    }
}
