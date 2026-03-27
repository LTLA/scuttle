#include "Rcpp.h"
#include "Rtatami.h"
#include "scuttle/linear_model_fit.h"

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <cstddef>

// [[Rcpp::export(rng=false)]]
Rcpp::RObject fit_linear_model (Rcpp::NumericMatrix qr, Rcpp::NumericVector qraux, Rcpp::RObject exprs, bool get_coefs, int nthreads) {
    scuttle::linear_model_fit fitter(qr, qraux);
    const int ncoefs = fitter.get_ncoefs();
    const int ncells = fitter.get_nobs();

    Rtatami::BoundNumericPointer eptr(exprs);
    const auto& emat = *(eptr->ptr);
    if (ncells != emat.ncol()) {
        throw std::runtime_error("number of rows of QR matrix not equal to number of cells");
    } else if (ncells == 0) {
        throw std::runtime_error("cannot compute variance for zero cells");
    }
    const int ngenes = emat.nrow();

    Rcpp::NumericVector means(ngenes), vars(ngenes);
    double* mptr = means.begin();
    double* vptr = vars.begin();
    Rcpp::NumericMatrix coefs((get_coefs ? ncoefs : 0), (get_coefs ? ngenes : 0));
    double* cptr = coefs.begin();

    tatami::parallelize([&](int, int start, int length) -> void {
        auto ext = tatami::consecutive_extractor<false>(emat, true, start, length);
        std::vector<double> buffer(ncells);

        for (int g = start, end = start + length; g < end; ++g) {
            auto ptr = ext->fetch(buffer.data());
            mptr[g] = std::accumulate(ptr, ptr + ncells, 0.0) / ncells; // compute mean BEFORE in-place mutation by multiply().

            tatami::copy_n(ptr, ncells, buffer.data());
            fitter.multiply(buffer.data());

            double curvar = 0;
            for (int c = ncoefs; c < ncells; ++c) { // only using the residual effects.
                curvar += buffer[c] * buffer[c];
            }
            curvar /= ncells - ncoefs;
            vptr[g] = curvar;

            if (get_coefs) { 
                fitter.solve(buffer.data());
                std::copy(buffer.begin(), buffer.begin() + ncoefs, cptr + sanisizer::product_unsafe<std::size_t>(g, ncoefs));
            }
        }
    }, ngenes, nthreads);

    Rcpp::List output(3);
    if (get_coefs) {
        output[0] = coefs;
    } else {
        output[0] = R_NilValue;
    }
    output[1] = means;
    output[2] = vars;

    return output;
}
