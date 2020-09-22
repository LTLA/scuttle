#include "Rcpp.h"
#include "beachmat3/beachmat.h"
#include <vector>

// [[Rcpp::export(rng=false)]]
Rcpp::RObject sum_row_counts (Rcpp::RObject counts, Rcpp::IntegerVector genes, Rcpp::IntegerVector runs) {
    auto mat = beachmat::read_lin_block(counts);

    const size_t ncells=mat->get_ncol();
    std::vector<double> holding(mat->get_nrow());
    const size_t nsummations=runs.size();
    Rcpp::NumericMatrix output(nsummations, ncells);

    for (size_t c=0; c<ncells; ++c) {
        auto it = mat->get_col(c, holding.data());
        auto outcol=output.column(c);
        auto oIt=outcol.begin();
        auto gIt=genes.begin();

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
