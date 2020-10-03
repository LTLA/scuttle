#include "Rcpp.h"
#include "beachmat3/beachmat.h"
#include "scuttle/downsample_vector.h"
#include <deque>

// [[Rcpp::export]]
Rcpp::RObject downsample_column(Rcpp::RObject input, Rcpp::NumericVector prop) {
    auto mat = beachmat::read_lin_block(input);
    const size_t ngenes = mat->get_nrow();
    const size_t ncells = mat->get_ncol();
    Rcpp::RNGScope _rng; 

    // Always return a sparse matrix for memory efficiency: we assume that
    // downsampling will increase sparsity sufficiently that, e.g., file-backed
    // matrices can be represented fully in memory.
    std::deque<std::pair<std::pair<int, int>, double> > store;

    if (mat->is_sparse()) {
        auto smat = beachmat::promote_to_sparse(mat);
        std::vector<int> work_i(ngenes);
        std::vector<double> work_x(ngenes);

        auto pIt=prop.begin();
        for (size_t i=0; i<ncells; ++i, ++pIt) {
            // note that the sparse_index contract for 'indices' guarantees
            // that elements are added to 'store' in an already-sorted manner.
            auto indices = smat->get_col(i, work_x.data(), work_i.data());

            scuttle::downsample_vector(indices.x, indices.x + indices.n, work_x.begin(), *pIt);
            for (size_t j = 0; j < indices.n; ++j) {
                if (work_x[j]) {
                    store.push_back(std::make_pair(std::make_pair(i, indices.i[j]), work_x[j]));
                }
            }
        }

    } else {
        std::vector<double> workspace(ngenes);
        auto pIt=prop.begin();
        for (size_t i=0; i<ncells; ++i, ++pIt) {
            auto ptr = mat->get_col(i, workspace.data());
            scuttle::downsample_vector(ptr, ptr + ngenes, workspace.data(), *pIt);

            for (size_t j = 0; j < ngenes; ++j) {
                if (workspace[j]) {
                    store.push_back(std::make_pair(std::make_pair(i, j), workspace[j]));
                }
            }
        }
    }

    return beachmat::as_gCMatrix<Rcpp::NumericVector>(ngenes, ncells, store);
}

// [[Rcpp::export]]
Rcpp::RObject downsample_matrix(Rcpp::RObject rmat, double total, double required) {
    auto mat = beachmat::read_lin_block(rmat);
    const size_t ngenes=mat->get_nrow();
    const size_t ncells=mat->get_ncol();

    Rcpp::RNGScope _rng; 
    scuttle::downsample_vector_part downsampler(total, required, false);
    double subtotal = 0;

    // Always return a sparse matrix for memory efficiency: we assume that
    // downsampling will increase sparsity sufficiently that, e.g., file-backed
    // matrices can be represented fully in memory.
    std::deque<std::pair<std::pair<int, int>, double> > store;

    if (mat->is_sparse()) {
        auto smat = beachmat::promote_to_sparse(mat);
        std::vector<int> work_i(ngenes);
        std::vector<double> work_x(ngenes);

        for (size_t i=0; i<ncells; ++i) {
            // note that the sparse_index contract for 'indices' guarantees
            // that elements are added to 'store' in an already-sorted manner.
            auto indices = smat->get_col(i, work_x.data(), work_i.data());

            // do this before downsampling; contents in ptr may be overwritten!
            subtotal += std::accumulate(indices.x, indices.x + indices.n, 0.0); 

            downsampler(indices.x, indices.x + indices.n, work_x.begin());
            for (size_t j = 0; j < indices.n; ++j) {
                if (work_x[j]) {
                    store.push_back(std::make_pair(std::make_pair(i, indices.i[j]), work_x[j]));
                }
            }
        }
    } else {
        std::vector<double> workspace(ngenes);
        for (size_t i=0; i<ncells; ++i) {
            auto ptr = mat->get_col(i, workspace.data());

            // do this before downsampling; contents in ptr may be overwritten!
            subtotal += std::accumulate(ptr, ptr + ngenes, 0.0); 

            downsampler(ptr, ptr + ngenes, workspace.data());

            for (size_t j = 0; j < ngenes; ++j) {
                if (workspace[j]) {
                    store.push_back(std::make_pair(std::make_pair(i, j), workspace[j]));
                }
            }
        }
    }

    return Rcpp::List::create(
        beachmat::as_gCMatrix<Rcpp::NumericVector>(ngenes, ncells, store),
        Rcpp::NumericVector::create(subtotal)
    );
}
