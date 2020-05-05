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
 * num_total: An integer scalar specifying the total number of all events.
 * num_sample: An integer scalar specifying the number of events to sample without replacement.
 * num_processed: An integer scalar specifying the number of events that have already been considered for selection.
 * num_selected: An integer scalar specifying the number of events that have already been selected.
 * 
 * Note that num_total may not be simply a sum of all values from [freqIt, freqEnd).
 * This is because we allow multiple applications of this function to sample without replacement from a series of vectors.
 * We keep track of 'num_processed' and 'num_selected' to ensure correct sampling when moving from one vector to another.
 */
template<class IN, class OUT> 
void downsample(IN freqIt, IN freqEnd, OUT freqOut, const bigint_t num_total, 
    const bigint_t num_sample, bigint_t& num_processed, bigint_t& num_selected) 
{        
    while (freqIt!=freqEnd && num_selected < num_sample) {
        // It is allowed for this function to modify in place,
        // so *freqIt must be read before *freqOut is modified. 
        const int full_count=*freqIt;
        auto& downsampled=(*freqOut=0);

        for (int i=0; i<full_count && num_sample > num_selected; ++i) {
            // Deciding whether or not to keep this instance of this event.
            // This is a safe way of computing NUM_YET_TO_SELECT/NUM_YET_TO_PROCESS > runif(1), 
            // avoiding issues with integer division.
            if ( (num_total - num_processed)*R::unif_rand() < num_sample - num_selected) {
                ++downsampled;
                ++num_selected;
            }
            ++num_processed;
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

    bigint_t num_total=std::round(total);
    bigint_t num_sample=std::round(prop*total);
    bigint_t num_processed=0;
    bigint_t num_selected=0;

    downsample(freqIt, freqEnd, freqOut, num_total, num_sample, num_processed, num_selected);
    return;
}

/* Progressively downsampling parts of a vector to a given proportion. */

class downsample_vector_part {
public:
    downsample_vector_part(double total, double prop) : num_total(std::round(total)), num_sample(std::round(prop * total)) {}

    template<class IN, class OUT> 
    void operator()(IN freqIt, IN freqEnd, OUT freqOut) { 
        downsample(freqIt, freqEnd, freqOut, num_total, num_sample, num_processed, num_selected);
        return;
    }
private:
    const bigint_t num_total=0, num_sample=0;
    bigint_t num_processed=0, num_selected=0;
};

}
