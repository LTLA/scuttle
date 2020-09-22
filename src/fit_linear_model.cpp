#include "Rcpp.h"
#include "beachmat3/beachmat.h"
#include "scuttle/linear_model_fit.h"

#include <stdexcept>
#include <algorithm>
#include <vector>

// [[Rcpp::export(rng=false)]]
Rcpp::RObject fit_linear_model (Rcpp::NumericMatrix qr, Rcpp::NumericVector qraux, Rcpp::RObject exprs, bool get_coefs) {
    // Setting up for QR-based multiplication.
    scuttle::linear_model_fit fitter(qr, qraux);
    const int ncoefs = fitter.get_ncoefs();
    const int ncells = fitter.get_nobs();

    auto emat = beachmat::read_lin_block(exprs);
    if (ncells != static_cast<int>(emat->get_ncol())) {
        throw std::runtime_error("number of rows of QR matrix not equal to number of cells");
    } else if (ncells == 0) {
        throw std::runtime_error("cannot compute variance for zero cells");
    }
    const size_t ngenes = emat->get_nrow();

    // Setting up output objects.
    Rcpp::NumericVector means(ngenes), vars(ngenes);
    auto mIt = means.begin(), vIt = vars.begin();
    std::vector<double> tmp(ncells);
    Rcpp::NumericMatrix coefs((get_coefs ? ncoefs : 0), (get_coefs ? ngenes : 0));
    auto cIt = coefs.begin();

    // Running through each gene and reporting its variance and mean.
    for (size_t s = 0; s < ngenes; ++s) {
        auto ptr = emat->get_row(s, tmp.data());
        if (ptr != tmp.data()) { // forcing a copy, just in case.
            std::copy(ptr, ptr + ncells, tmp.data());
        }
        (*mIt) = std::accumulate(tmp.begin(), tmp.end(), 0.0)/ncells;
        ++mIt;

        fitter.multiply(tmp.data());
        double& curvar = (*vIt);
        for (auto tIt = tmp.begin() + ncoefs; tIt != tmp.end(); ++tIt) { // only using the residual effects.
            curvar += (*tIt) * (*tIt);
        }
        curvar /= ncells - ncoefs;
        ++vIt;

        if (get_coefs) { 
            fitter.solve(tmp.data());
            std::copy(tmp.begin(), tmp.begin() + ncoefs, cIt);
            cIt += ncoefs;
        }
    }
    
    if (get_coefs) {
        return Rcpp::List::create(coefs, means, vars);
    } else {
        return Rcpp::List::create(R_NilValue, means, vars);
    }
}
