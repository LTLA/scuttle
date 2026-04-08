#ifndef SCUTTLE_DOWNSAMPLE_H
#define SCUTTLE_DOWNSAMPLE_H

#include "Rcpp.h"

namespace scuttle {

/* Quick and dirty check for whether a value is too large to have retained integer precision.
 * This is important for downsampling where we need to calculate the remaining total/required numbers exactly.
 * We check this by adding 1, then subtracting the original value, and seeing if we get 1 back;
 * just beyond 2^53, addition of 1 will be rounded to the nearest even number, and then multiples of 4, 8, etc. 
 *
 * Technically, this also fails for y = 2^53, which is slightly over-conservative.
 * However, we don't know whether any rounding was performed to get to 2^53, so it's possible that integer precision was already lost to obtain 'y' in the first place.
 * If y = 2^53 - 1, we can be sure that no rounding was performed (otherwise we would have gotten 2^53) and thus no overflow will have occurred.
 */
inline bool too_large_for_integer_precision(double y) {
    double p1 = y + 1; 
    double diff = p1 - y;
    return !(diff == 1); // use !() so that NaN diffs return 'true'.
}

/* This function considers sampling events without replacement from a vector.
 * Here, though, the vector contains frequencies of events rather than the events themselves.
 * 
 * freqIt: An iterator pointing to the start of the frequency vector, should contain integers.
 * freqEnd: An interator pointing to the end of the frequency vector, should contain integers.
 * freqOut: An iterator pointing to an output vector, indicating how many instances of each event have been sampled.
 * remaining_total: Number of instances yet to be considered for sampling, should be integer.
 * remaining_required: Number of events that still need to be sampled, should be integer and no greater than remaining_total.
 * 
 * Note that remaining_total's initial value may not be simply a sum of all values from [freqIt, freqEnd).
 * This is because we allow multiple applications of this function to sample without replacement from a series of vectors.
 * We keep track of 'remaining_total' and 'remaining_required' to ensure correct sampling when moving from one vector to another.
 *
 * It is assumed that R's RNG state is already initialized prior to calling this function, via `Rcpp::RNGScope()` or otherwise.
 */
template<class Input_, class Output_> 
void downsample(Input_ freqIt, Input_ freqEnd, Output_ freqOut, double& remaining_total, double& remaining_required) {        
    while (freqIt != freqEnd && remaining_required) {
        const double current = *freqIt;
        remaining_total -= current;
        const double chosen = Rf_rhyper(current, remaining_total, remaining_required);
        remaining_required -= chosen;

        // It is allowed for this function to modify in place, so *freqIt must be read before *freqOut is modified. 
        *freqOut = chosen;

        ++freqIt;
        ++freqOut;
    }

    // Mopping up the leftovers.
    while (freqIt != freqEnd) {
        ++freqIt;
        *freqOut = 0;
        ++freqOut;
    }
}  

}

#endif
