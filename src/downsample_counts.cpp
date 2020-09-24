#include "Rcpp.h"
#include "beachmat3/beachmat.h"
#include "scuttle/downsample_vector.h"

// [[Rcpp::export]]
Rcpp::RObject downsample_column(Rcpp::RObject input, Rcpp::NumericVector prop) {
    auto mat = beachmat::read_lin_block(input);
    const size_t ngenes = mat->get_nrow();
    const size_t ncells = mat->get_ncol();
    Rcpp::RNGScope _rng; 

    if (mat->is_sparse()) {
        auto smat = beachmat::promote_to_sparse(mat);
        std::vector<int> work_i(ngenes);
        std::vector<double> work_x(ngenes);
        std::map<std::pair<int, int>, double> store;

        auto pIt=prop.begin();
        for (size_t i=0; i<ncells; ++i, ++pIt) {
            auto indices = smat->get_col(i, work_x.data(), work_i.data());
            scuttle::downsample_vector(indices.x, indices.x + indices.n, work_x.begin(), *pIt);
            for (size_t j = 0; j < indices.n; ++j) {
                if (work_x[j]) {
                    store[std::make_pair(i, indices.i[j])] = work_x[j];
                }
            }
        }

        return beachmat::as_gCMatrix<Rcpp::NumericVector>(ngenes, ncells, store);

    } else {
        Rcpp::NumericMatrix output(ngenes, ncells);

        auto pIt=prop.begin();
        for (size_t i=0; i<ncells; ++i, ++pIt) {
            auto curcol = output.column(i);
            auto ptr = mat->get_col(i, curcol.begin());
            scuttle::downsample_vector(ptr, ptr + ngenes, curcol.begin(), *pIt);
        }

        return output;
    }
}

// [[Rcpp::export]]
Rcpp::RObject downsample_matrix(Rcpp::RObject rmat, double total, double required) {
    auto mat = beachmat::read_lin_block(rmat);
    const size_t ngenes=mat->get_nrow();
    const size_t ncells=mat->get_ncol();

    Rcpp::RNGScope _rng; 
    scuttle::downsample_vector_part downsampler(total, required, false);
    Rcpp::RObject output;
    double subtotal = 0;

    if (mat->is_sparse()) {
        auto smat = beachmat::promote_to_sparse(mat);
        std::vector<int> work_i(ngenes);
        std::vector<double> work_x(ngenes);
        std::map<std::pair<int, int>, double> store;

        for (size_t i=0; i<ncells; ++i) {
            auto indices = smat->get_col(i, work_x.data(), work_i.data());

            // do this before downsampling; contents in ptr may be overwritten!
            subtotal += std::accumulate(indices.x, indices.x + indices.n, 0.0); 

            downsampler(indices.x, indices.x + indices.n, work_x.begin());
            for (size_t j = 0; j < indices.n; ++j) {
                if (work_x[j]) {
                    store[std::make_pair(i, indices.i[j])] = work_x[j];
                }
            }
        }

        output = beachmat::as_gCMatrix<Rcpp::NumericVector>(ngenes, ncells, store);

    } else {
        Rcpp::NumericMatrix outmat(ngenes, ncells);
        for (size_t i=0; i<ncells; ++i) {
            auto curcol = outmat.column(i);
            auto ptr = mat->get_col(i, curcol.begin());

            // do this before downsampling; contents in ptr may be overwritten!
            subtotal += std::accumulate(ptr, ptr + ngenes, 0.0); 

            downsampler(ptr, ptr + ngenes, curcol.begin());
        }

        output = outmat;
    }

    return Rcpp::List::create(output, Rcpp::NumericVector::create(subtotal));
}
