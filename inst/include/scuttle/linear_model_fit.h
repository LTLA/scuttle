#ifndef LINEAR_MODEL_FIT_H 
#define LINEAR_MODEL_FIT_H

#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#ifndef FCONE
#define FCONE
#endif
#include "Rcpp.h"

#include <stdexcept>

namespace scuttle {

class QR_multiplier {
public:
    QR_multiplier(Rcpp::NumericMatrix qr, Rcpp::NumericVector qraux, const char tr='T') :
        QR(qr), AUX(qraux), qrptr(QR.begin()), qxptr(AUX.begin()), nobs(QR.nrow()), ncoef(QR.ncol()), trans(tr) 
    {
        if (AUX.size()!=ncoef) { 
            throw std::runtime_error("QR auxiliary vector should be of length 'ncol(Q)'"); 
        }

        work.resize(nobs);
        double tmpwork=0;

        F77_CALL(dormqr)(&side, &trans, &nobs, &ncol, &ncoef, qrptr, &nobs, 
            qxptr, work.data(), &nobs, &tmpwork, &lwork, &info FCONE FCONE); 
        if (info) { 
            throw std::runtime_error("workspace query failed for 'dormqr'");
        }

        lwork=int(tmpwork+0.5);
        work.resize(lwork);
        return;
    }

    void multiply(double* rhs) {
        F77_CALL(dormqr)(&side, &trans, &nobs, &ncol, &ncoef, qrptr, &nobs, 
            qxptr, rhs, &nobs, work.data(), &lwork, &info FCONE FCONE); 
        if (info) { 
            throw std::runtime_error("residual calculations failed for 'dormqr'");
        }
        return;
    }

    int get_nobs() const {
        return nobs;
    }

    int get_ncoefs() const {
        return ncoef;
    }
protected:
    Rcpp::NumericMatrix QR;
    Rcpp::NumericVector AUX;
    const double* qrptr, * qxptr;
    const int nobs, ncoef;
    const char trans;

    int info=0, lwork=-1;
    std::vector<double> work;

    const int ncol=1;
    const char side='L';
};

class linear_model_fit : public QR_multiplier { 
public:
    linear_model_fit(Rcpp::NumericMatrix qr, Rcpp::NumericVector qraux) : QR_multiplier(qr, qraux, 'T') {}

    void solve(double* rhs) {
        F77_CALL(dtrtrs)(&uplo, &xtrans, &diag, &(this->ncoef), &(this->ncol), this->qrptr, &(this->nobs), 
            rhs, &(this->nobs), &(this->info));

        if (this->info) { 
            throw std::runtime_error("coefficient calculations failed for 'dtrtrs'");
        }
        return;
    }
private:
    const char uplo='U', xtrans='N', diag='N';
};

}

#endif
