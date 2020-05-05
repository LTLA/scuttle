#include "Rcpp.h"
#include <cmath>
#include <cstdint>

namespace scuttle {

namespace {

// Using big integers to avoid overflow when summing across a large matrix.

typedef std::uint64_t bigint_t;

/* This function considers sampling events without replacement from a vector.
 * Here, though, the vector contains frequencies of events rather than the events themselves.
 * The sampling scheme is adapted from John D. Cook, https://stackoverflow.com/a/311716/15485.
 * 
 * freqIt: An iterator pointing to the start of the frequency vector.
 * freqEnd: An interator pointing to the end of the frequency vector.
 * freqOut: An iterator pointing to an output vector, indicating how many instances of each event have been sampled.
 * remaining_total: An integer scalar specifying the remaining number of instances yet to be considered for sampling.
 * remaining_required: An integer scalar specifying the number of events that still need to be sampled.
 * 
 * Note that remaining_total's initial value may not be simply a sum of all values from [freqIt, freqEnd).
 * This is because we allow multiple applications of this function to sample without replacement from a series of vectors.
 * We keep track of 'remaining_total' and 'remaining_required' to ensure correct sampling when moving from one vector to another.
 */
template<class IN, class OUT> 
void downsample(IN freqIt, IN freqEnd, OUT freqOut, bigint_t& remaining_total, bigint_t& remaining_required) 
{        
    while (freqIt!=freqEnd && remaining_required) {
        // It is allowed for this function to modify in place,
        // so *freqIt must be read before *freqOut is modified. 
        const int full_count=*freqIt;
        auto& downsampled=(*freqOut=0);

        for (int i=0; i<full_count && remaining_required; ++i) {
            // Deciding whether or not to keep this instance of this event.
            // This is a safe way of computing NUM_YET_TO_SELECT/NUM_YET_TO_PROCESS > runif(1), 
            // avoiding issues with integer division.
            if ( remaining_total*R::unif_rand() < remaining_required) {
                ++downsampled;
                --remaining_required;
            }
            --remaining_total;
        }
     
        ++freqIt;
        ++freqOut;
    }

    // Mopping up the leftovers.
    while (freqIt!=freqEnd) {
        ++freqIt;
        *freqOut=0;
        ++freqOut;
    }

    return;
}  

}

/* Downsampling a vector to a given proportion. */

template<class IN, class OUT> 
void downsample_vector(IN freqIt, IN freqEnd, OUT freqOut, double prop) {
    double total=0;
    {
        auto freqCopy=freqIt;
        while (freqCopy!=freqEnd) {
            total+=*freqCopy;
            ++freqCopy;
        }
    }

    bigint_t remaining_total=std::round(total);
    bigint_t remaining_required=std::round(std::min(1.0, prop)*total);
    downsample(freqIt, freqEnd, freqOut, remaining_total, remaining_required);
    return;
}

/* Progressively downsampling parts of a vector to a given proportion. */

class downsample_vector_part {
public:
    downsample_vector_part(double total, double prop) : 
        remaining_total(std::round(total)), 
        remaining_required(std::round(std::min(1.0, prop) * total)) {}

    template<class IN, class OUT> 
    void operator()(IN freqIt, IN freqEnd, OUT freqOut) { 
        downsample(freqIt, freqEnd, freqOut, remaining_total, remaining_required);
        return;
    }
private:
    bigint_t remaining_total, remaining_required;
};

}
