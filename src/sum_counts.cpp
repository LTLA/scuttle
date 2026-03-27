#include "Rtatami.h"

#include <vector>

// [[Rcpp::export(rng=false)]]
Rcpp::RObject sum_row_counts (Rcpp::RObject counts, Rcpp::IntegerVector genes, Rcpp::IntegerVector runs) {
    Rtatami::BoundNumericPointer ptr(counts);
    const auto& mat = *(ptr->ptr);
    const auto ncells = mat.ncol();
    auto ext = tatami::consecutive_extractor<false>(mat, false, 0, mat.ncol());

    std::vector<double> holding(mat.nrow());
    const auto nsummations = runs.size();
    Rcpp::NumericMatrix output(nsummations, ncells);

    // Not really worth parallelizing it or seting up sparse extractors...
    // we're going to throw it out anyway, in favor of scrapper::aggregateAcrossGenes.
    for (int c = 0; c < ncells; ++c) {
        auto it = ext->fetch(holding.data());
        auto outcol = output.column(c);
        auto oIt = outcol.begin();
        auto gIt = genes.begin();

        for (auto s : runs) {
            while (s > 0) {
                (*oIt) += *(it + *gIt - 1); // get back to zero indexing.
                --s;
                ++gIt;
            }
            ++oIt;
        }
    }

    return output;
} 
