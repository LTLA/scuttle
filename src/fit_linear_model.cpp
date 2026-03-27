#include "Rcpp.h"
#include "Rtatami.h"
#include "scuttle/linear_model_fit.h"

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <cstddef>
#include <optional>

// [[Rcpp::export(rng=false)]]
Rcpp::RObject fit_linear_model (Rcpp::NumericMatrix qr, Rcpp::NumericVector qraux, Rcpp::RObject exprs, bool get_coefs, int nthreads) {
    const int nobs = sanisizer::cast<int>(qr.nrow());
    const int ncoefs = sanisizer::cast<int>(qr.ncol());

    const double* qrptr = qr.begin();
    const double* qxptr = qraux.begin();
    if (qraux.size() != ncoefs) {
        throw std::runtime_error("QR auxiliary vector should be of length 'ncol(Q)'");
    }

    constexpr char trans = 'T';
    constexpr int ncol = 1;
    constexpr char side = 'L';
    constexpr char uplo = 'U', xtrans = 'N', diag = 'N';

    int lwork = -1;
    {
        std::vector<double> work(nobs);
        double tmpwork = 0;

        int info = 0;
        F77_CALL(dormqr)(
            &side,
            &trans,
            &nobs,
            &ncol,
            &ncoefs,
            qrptr,
            &nobs,
            qxptr,
            work.data(),
            &nobs,
            &tmpwork,
            &lwork,
            &info
            FCONE
            FCONE
        );
        if (info) {
            throw std::runtime_error("workspace query failed for 'dormqr'");
        }

        lwork = tmpwork + 0.5; // rounding up.
    }

    Rtatami::BoundNumericPointer eptr(exprs);
    const auto& emat = *(eptr->ptr);
    if (nobs != emat.ncol()) {
        throw std::runtime_error("number of rows of QR matrix not equal to number of cells");
    } else if (nobs == 0) {
        throw std::runtime_error("cannot compute variance for zero cells");
    }
    const int ngenes = emat.nrow();

    Rcpp::NumericVector means(ngenes), vars(ngenes);
    double* mptr = means.begin();
    double* vptr = vars.begin();
    Rcpp::NumericMatrix coefs((get_coefs ? ncoefs : 0), (get_coefs ? ngenes : 0));
    double* cptr = coefs.begin();

    tatami::parallelize([&](int t, int start, int length) -> void {
        std::vector<double> work(lwork);
        auto ext = tatami::consecutive_extractor<false>(emat, true, start, length);
        std::vector<double> buffer(nobs);

        for (int g = start, end = start + length; g < end; ++g) {
            auto ptr = ext->fetch(buffer.data());
            mptr[g] = std::accumulate(ptr, ptr + nobs, 0.0) / nobs; // compute mean BEFORE in-place mutation by multiply().
            tatami::copy_n(ptr, nobs, buffer.data());

            int info = 0;
            F77_CALL(dormqr)(
                &side,
                &trans,
                &nobs,
                &ncol,
                &ncoefs,
                qrptr,
                &nobs, 
                qxptr,
                buffer.data(),
                &nobs,
                work.data(),
                &lwork,
                &info
                FCONE
                FCONE
            ); 
            if (info) { 
                throw std::runtime_error("residual calculations failed for 'dormqr'");
            }

            double curvar = 0;
            for (int c = ncoefs; c < nobs; ++c) { // only using the residual effects.
                curvar += buffer[c] * buffer[c];
            }
            curvar /= nobs - ncoefs;
            vptr[g] = curvar;

            if (get_coefs) { 
                F77_CALL(dtrtrs)(
                    &uplo,
                    &xtrans,
                    &diag,
                    &ncoefs,
                    &ncol,
                    qrptr,
                    &nobs, 
                    buffer.data(),
                    &nobs,
                    &info
                    FCONE
                    FCONE
                    FCONE
                );
                if (info) { 
                    throw std::runtime_error("coefficient calculations failed for 'dtrtrs'");
                }

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
