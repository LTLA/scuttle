#include "Rcpp.h"
#include "beachmat3/beachmat.h"

#include <vector>
#include <algorithm>

template<typename T, class IT>
void compute_cumsum (T* it, size_t n, const Rcpp::IntegerVector& top, IT out) {
    const size_t ntop=top.size();
    if (ntop == 0) {
        return;
    }

    std::partial_sort(it, it + std::min(n, size_t(top[ntop-1])), it + n, std::greater<T>());
    size_t x = 0;
    T accumulated = 0;

    for (size_t target_index : top) {
        while (x < target_index && x < n) { // '<' as top contains 1-based indices.
            accumulated += *(it+x);
            ++x;
        }
        (*out)=accumulated;
        ++out;
    }

    return;
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix cumulative_prop(Rcpp::RObject input, Rcpp::IntegerVector top) {
    auto mat = beachmat::read_lin_block(input);
    const size_t ngenes = mat->get_nrow();
    const size_t ncells = mat->get_ncol();
    const size_t ntop = top.size();    
    Rcpp::NumericMatrix output(ntop, ncells);

    if (mat->is_sparse()) {
        auto smat = beachmat::promote_to_sparse(mat);
        std::vector<int> work_i(ngenes);
        std::vector<double> work_x(ngenes);

        for (size_t c = 0; c < ncells; ++c) {
            auto indices = smat->get_col(c, work_x.data(), work_i.data());
            if (indices.x != work_x.data()) {
                std::copy(indices.x, indices.x + indices.n, work_x.begin());
            }
            auto outcol = output.column(c);
            compute_cumsum(work_x.data(), indices.n, top, outcol.begin());
        }

    } else {
        std::vector<int> work(ngenes);

        for (size_t c = 0; c < ncells; ++c) {
            auto ptr = mat->get_col(c, work.data());
            if (ptr != work.data()) {
                std::copy(ptr, ptr + ngenes, work.begin());
            }
            auto outcol = output.column(c);
            compute_cumsum(work.data(), ngenes, top, outcol.begin());
        }
    }

    return output;
}
