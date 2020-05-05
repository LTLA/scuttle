#include "Rcpp.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "scuttle/downsample_vector.h"

template <class M, class O>
Rcpp::RObject downsample_column_internal(Rcpp::RObject input, Rcpp::NumericVector prop) {
    auto mat=beachmat::create_matrix<M>(input);
    const size_t ngenes=mat->get_nrow();
    const size_t ncells=mat->get_ncol();
    typename M::vector tmp(ngenes);

    auto otype=beachmat::output_param(mat.get());
    auto output=beachmat::create_output<O>(ngenes, ncells, otype);

    auto pIt=prop.begin();
    Rcpp::RNGScope _rng; 

    for (size_t i=0; i<ncells; ++i, ++pIt) {
        mat->get_col(i, tmp.begin());
        scuttle::downsample_vector(tmp.begin(), tmp.end(), tmp.begin(), *pIt);
        output->set_col(i, tmp.begin());
    }

    return output->yield();
}

// [[Rcpp::export]]
Rcpp::RObject downsample_column(Rcpp::RObject input, Rcpp::NumericVector prop) {
    int rtype=beachmat::find_sexp_type(input);
    if (rtype==INTSXP) {
        return downsample_column_internal<beachmat::integer_matrix, beachmat::integer_output>(input, prop);
    } else {
        return downsample_column_internal<beachmat::numeric_matrix, beachmat::numeric_output>(input, prop);
    }
}


template <class M, class O>
Rcpp::RObject downsample_matrix_internal(Rcpp::RObject input, double total, double prop) {
    auto mat=beachmat::create_matrix<M>(input);
    const size_t ngenes=mat->get_nrow();
    const size_t ncells=mat->get_ncol();
    typename M::vector tmp(ngenes);

    auto otype=beachmat::output_param(mat.get());
    auto output=beachmat::create_output<O>(ngenes, ncells, otype);

    scuttle::downsample_vector_part downsampler(total, prop);
    Rcpp::RNGScope _rng; 

    for (size_t i=0; i<ncells; ++i) {
        mat->get_col(i, tmp.begin());
        downsampler(tmp.begin(), tmp.end(), tmp.begin());
        output->set_col(i, tmp.begin());
    }

    return output->yield();
}

//[[Rcpp::export]]
Rcpp::RObject downsample_matrix(Rcpp::RObject rmat, double total, double prop) {
    int rtype=beachmat::find_sexp_type(rmat);
    if (rtype==INTSXP) {
        return downsample_matrix_internal<beachmat::integer_matrix, beachmat::integer_output>(rmat, total, prop);
    } else {
        return downsample_matrix_internal<beachmat::numeric_matrix, beachmat::numeric_output>(rmat, total, prop);
    }
}
