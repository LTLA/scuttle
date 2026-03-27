#include "Rcpp.h"
#include "Rtatami.h"
#include "tatami_stats/tatami_stats.hpp"

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <deque>

template <typename T>
void quick_rotate (std::deque<T>& dq) {
    dq.push_back(dq.front());
    dq.pop_front();
    return;
}

/*** A function to estimate the pooled size factors and construct the linear equations. ***/

// [[Rcpp::export(rng=false)]]
Rcpp::List pool_size_factors (Rcpp::RObject exprs, Rcpp::NumericVector pseudo_cell, Rcpp::IntegerVector order, Rcpp::IntegerVector pool_sizes) {
    Rtatami::BoundNumericPointer eptr(exprs);
    const auto& emat = *(eptr->ptr);
    const int ngenes = emat.nrow();
    const int ncells = emat.ncol();

    /******************************
     * Applying checks on inputs. *
     ******************************/

    if (ncells == 0) { 
        throw std::runtime_error("at least one cell required for normalization"); 
    }
    if (ngenes == 0) {
        throw std::runtime_error("insufficient features for median calculations");
    }

    // Checking the input sizes.
    const auto nsizes = pool_sizes.size();
    if (nsizes == 0) {
        return Rcpp::List::create(Rcpp::IntegerVector(0), Rcpp::IntegerVector(0), Rcpp::NumericVector(0));
    }

    int last_size = -1, total_size = 0;
    for (auto s : pool_sizes) { 
        if (s < 1 || s > ncells) { 
            throw std::runtime_error("each element of sizes should be within [1, number of cells]"); 
        }
        if (s < last_size) { 
            throw std::runtime_error("sizes should be sorted"); 
        }
        total_size+=s;
        last_size=s;
    }

    if (!sanisizer::is_equal(ngenes, pseudo_cell.size())) { 
        throw std::runtime_error("length of pseudo-cell vector is not the same as the number of cells"); 
    }

    // THis is equivalent to 'order.size() < ncells * 2 - 1', but without the risk of overflow in the RHS.
    if (
        sanisizer::is_less_than(order.size(), ncells) || 
        sanisizer::is_less_than(order.size() - ncells, ncells - 1)
    ) { 
        throw std::runtime_error("ordering vector is too short for number of cells"); 
    }
    for (auto o : order) { 
        if (o < 0 || o >= ncells) { 
            throw std::runtime_error("elements of ordering vector are out of range");
        }
    }

    /*********************************
     * Setting up the storage space. *
     *********************************/

    std::deque<const double*> xptrs(last_size);
    std::deque<double*> work_xptrs(last_size);
    std::vector<double> work_x(sanisizer::product<typename std::vector<double>::size_type>(last_size, ngenes)); // workspace, should not be referenced except by work_xptrs.
    for (int i = 0; i < last_size; ++i) {
        work_xptrs[i] = work_x.data() + sanisizer::product_unsafe<std::size_t>(i, ngenes);
    }

    std::deque<const int*> iptrs;
    std::deque<int*> work_iptrs;
    std::vector<int> work_i;
    std::deque<int> nnzero;

    const bool is_sparse = emat.is_sparse();
    if (is_sparse) {
        iptrs.resize(last_size);
        work_iptrs.resize(last_size);
        work_i.resize(sanisizer::product<typename std::vector<int>::size_type>(last_size, ngenes)); // workspace, should not be referenced except by work_iptrs.
        nnzero.resize(last_size);
        for (int i = 0; i < last_size; ++i) {
            work_iptrs[i] = work_i.data() + sanisizer::product_unsafe<std::size_t>(i, ngenes);
        }
    }

    std::unique_ptr<tatami::OracularSparseExtractor<double, int> > sparse_ext;
    std::unique_ptr<tatami::OracularDenseExtractor<double, int> > dense_ext;
    auto oracle = std::make_unique<tatami::FixedViewOracle<int> >(static_cast<int*>(order.begin()), order.size());
    if (is_sparse) {
        sparse_ext = emat.sparse_column(std::move(oracle));
    } else {
        dense_ext = emat.dense_column(std::move(oracle));
    }
     
    // The first vector is unfilled as it gets dropped and refilled in the first iteration anyway.
    // We also trigger generation of indices so that each element doesn't generate its own indices.
    for (int s = 1; s < last_size; ++s) {
        if (is_sparse) {
            auto idxs = sparse_ext->fetch(work_xptrs[s], work_iptrs[s]);
            xptrs[s] = idxs.value;
            iptrs[s] = idxs.index;
            nnzero[s] = idxs.number;
        } else {
            xptrs[s] = dense_ext->fetch(work_xptrs[s]);
        }
    }

    // Setting up the output vectors and other bits and pieces.
    auto num_length = sanisizer::product<decltype(std::declval<Rcpp::IntegerVector>().size())>(total_size, ncells);
    Rcpp::IntegerVector row_num(num_length), col_num(num_length);
    auto fac_length = sanisizer::product<decltype(std::declval<Rcpp::NumericVector>().size())>(nsizes, ncells);
    Rcpp::NumericVector pool_factor(fac_length);

    std::vector<double> combined(ngenes), ratios(ngenes);

    /******************************
     * Performing the iterations. *
     ******************************/

    auto rowIt = row_num.begin(), colIt = col_num.begin();
    auto orIt = order.begin();
    const bool is_even = bool(ngenes%2==0);
    const int halfway = int(ngenes/2);

    // Running through the sliding windows.
    for (int win=0; win < ncells; ++win) {
        std::fill(combined.begin(), combined.end(), 0);

        /* Rotating; effectively moves the first element of 'collected' to the end.
         * The is the same as shifting the column that we've moved past to the end,
         * and then overwriting it with the next column.
         */
        quick_rotate(xptrs);
        quick_rotate(work_xptrs);

        if (!is_sparse) {
            xptrs.back() = dense_ext->fetch(work_xptrs.back());
        } else {
            quick_rotate(iptrs);
            quick_rotate(work_iptrs);
            quick_rotate(nnzero);

            auto idxs = sparse_ext->fetch(work_xptrs.back(), work_iptrs.back());
            xptrs.back() = idxs.value;
            iptrs.back() = idxs.index;
            nnzero.back() = idxs.number;
        }

        int index = 0;
        int rownum = win; // Setting the row so that all pools with the same SIZE form consecutive equations.
        for (auto psIt = pool_sizes.begin(); psIt != pool_sizes.end(); ++psIt, rownum += ncells) { 
            const int& SIZE = (*psIt);
            std::fill(rowIt, rowIt + SIZE, rownum);
            rowIt += SIZE;
            std::copy(orIt, orIt + SIZE, colIt);
            colIt += SIZE;

            for (; index<SIZE; ++index) {
                auto val = xptrs[index];

                if (is_sparse) {
                    auto n = nnzero[index];
                    auto rows = iptrs[index];
                    for (int i = 0; i < n; ++i, ++val, ++rows) {
                        combined[*rows] += *val;
                    }
                } else { 
                    for (auto cIt = combined.begin(); cIt != combined.end(); ++cIt, ++val) {
                        (*cIt) += (*val);
                    }
                }
            }
           
            // Computing the median ratio against the reference.
            auto rIt = ratios.begin(), cIt = combined.begin();
            for (auto pcIt = pseudo_cell.begin(); pcIt != pseudo_cell.end(); ++pcIt, ++rIt, ++cIt) {
                (*rIt) = *cIt / *pcIt;
            }

            pool_factor[rownum] = tatami_stats::medians::direct(ratios.data(), ratios.size(), false); // zeros in reference already stripped out.
        }
    }    

    return Rcpp::List::create(row_num, col_num, pool_factor);
}
