#include "Rcpp.h"
#include "beachmat3/beachmat.h"
#include <vector>

// 'counts' is an abundance matrix
// 'genes' is an integer vector that includes indices for each gene (feature), ordered by group (features from the first group come first, followed by the next group, etc)
// 'runs' specifies the number of genes in each group, allowing us to split the features in 'genes' into distinct groups for summing
// [[Rcpp::export(rng=false)]]
Rcpp::RObject sum_row_counts (Rcpp::RObject counts, Rcpp::IntegerVector genes, Rcpp::IntegerVector runs, bool average = false, bool na_rm = false) {
    // Read the matrix from the 'counts' R object into beachmat format for efficient processing
    auto mat = beachmat::read_lin_block(counts);

    const size_t ncells=mat->get_ncol(); // Get the number of columns (cells) in the matrix
    std::vector<double> holding(mat->get_nrow()); // Vector to hold values for each column of the abundance matrix
    const size_t nsummations=runs.size(); // Get the number of gene groups (runs) to process
    Rcpp::NumericMatrix output(nsummations, ncells); // Initialize an output matrix for storing results.

    // Loop over each column of the abundance matrix
    for (size_t c=0; c<ncells; ++c) {
        auto it = mat->get_col(c, holding.data()); // Get the values of the current column and store them in holding
        auto outcol=output.column(c); // Access the current output column to store results for this cell
        auto oIt=outcol.begin(); // Initialize iterator for the elements/cells of single output column.
        auto gIt=genes.begin(); // Initialize iterator for the gene indices.

        // Loop over each grouping of genes
        for (auto s : runs) {
            int run_count = 0; // Track number of genes in the current run
            double run_sum = 0.0; // Track the sum of gene abundance in the current run
            
            // This loop iterates over the specified genes for the current run in the current cell (column), summing their abundance values
            while (s > 0) {
                // Get the abundance of single gene in this column.
                // '*gIt' is the 1-based index for the gene; subtract 1 for 0-based indexing used in C++
                double value = *(it + *gIt - 1);
                
                // For possible average calculation, get the number of genes that have non-NA value.
                // If the value is non-NA, we increase the number of detected genes by one.
                if ( !(std::isnan(value) || value == NA_INTEGER) ) {
                    ++run_count;
                }
                // Replace NA value with 0. This can be done as we are summing values; it does not affect the result.
                if ( na_rm && (std::isnan(value) || value == NA_INTEGER) ) {
                    value = 0;
                }
                
                run_sum += value; // Add the gene abundance to the current group sum
                --s; // Decrement 's' to track how many genes are left to process in this group
                ++gIt; // Move to the next gene index in the 'genes' vector.
            }
            // Calculate average abundance of this group
            if (average && run_count > 0) {
                run_sum /= run_count;
            }
            
            (*oIt) = run_sum; // Add the abundance to the result matrix in specified cell
            ++oIt; // Move to next run/group, i.e., next cell of result matrix in specified column
        }
    }

    return output;
}
