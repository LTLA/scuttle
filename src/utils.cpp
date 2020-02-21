#include "utils.h"

#include <sstream>
#include <stdexcept>

// Checking subsetting vector.

Rcpp::IntegerVector process_subset_vector(Rcpp::RObject subset, int upper, bool zero_indexed) {
    if (subset.sexp_type()!=INTSXP) { 
        throw std::runtime_error("subset vector must be an integer vector");
    }  
    Rcpp::IntegerVector sout(subset);

    if (!zero_indexed) {
        sout=Rcpp::clone(sout);
        for (auto& s : sout) { --s; }
    }

    for (auto sIt=sout.begin(); sIt!=sout.end(); ++sIt) {
        if (*sIt < 0 || *sIt >= upper) { 
            throw std::runtime_error("subset indices out of range");
        }
    }
    return sout;
}
