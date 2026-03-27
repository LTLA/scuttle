#include "Rcpp.h"
#include "Rtatami.h"

#include <vector>
#include <algorithm>
#include <cstddef>
#include <cstdint>

template<typename Type_, class Iterator_>
void compute_cumsum(Type_* it, std::size_t n, const Rcpp::IntegerVector& top, Iterator_ out) {
    const auto ntop = top.size();
    if (ntop == 0) {
        return;
    }

    std::partial_sort(it, it + std::min(n, static_cast<std::size_t>(top[ntop-1])), it + n, std::greater<Type_>());
    std::size_t x = 0;
    Type_ accumulated = 0;

    for (std::size_t target_index : top) {
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
    Rtatami::BoundNumericPointer ptr(input);
    const auto& mat = *(ptr->ptr);

    const int ngenes = mat.nrow();
    const int ncells = mat.ncol();
    const int ntop = top.size();    
    Rcpp::NumericMatrix output(ntop, ncells);

    // Not really worth parallelizing it or seting up sparse extractors... we're going to throw it out anyway.
    auto ext = tatami::consecutive_extractor<false>(mat, false, 0, ncells);
    std::vector<double> work(ngenes);
    for (int c = 0; c < ncells; ++c) {
        auto ptr = ext->fetch(work.data());
        tatami::copy_n(ptr, ngenes, work.data());
        auto outcol = output.column(c);
        compute_cumsum(work.data(), ngenes, top, outcol.begin());
    }

    return output;
}
