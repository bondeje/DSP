#ifndef RTFILTER_H_
#define RTFILTER_H_

/* TODO: need to fix so that I can export these functions to the dll*/

#include <RTFilter.h>

static inline double * FilterBank_get_b(FilterBank * fb) {
    return fb->b;
}
static inline size_t FilterBank_get_b_size(FilterBank * fb) {
    return fb->nb;
}
static inline size_t FIR_order(FilterBank * fb) {
    return fb->nb;
}
static inline double * IIRFilterBank_get_b(IIRFilterBank * ifb) {
    return FilterBank_get_b(&ifb->fb);
}
static inline size_t IIRFilterBank_get_b_size(IIRFilterBank * ifb) {
    return FilterBank_get_b_size(&ifb->fb);
}
static inline double * IIRFilterBank_get_a(IIRFilterBank * ifb) {
    return IIRFilterBank_get_b(ifb) + IIRFilterBank_get_b_size(ifb);
}
static inline size_t IIRFilterBank_get_a_size(IIRFilterBank * ifb) {
    return ifb->na;
}
static inline double * RTFIRFilter_get_b(RTFIRFilter * rtff) {
    return rtff->fb.b;
}
static inline size_t RTFIRFilter_get_b_size(RTFIRFilter * rtff) {
    return rtff->fb.nb;
}
static inline double * RTIIRFilter_get_b(RTIIRFilter * rtif) {
    return rtif->ifb.fb.b;
}
static inline size_t RTIIRFilter_get_b_size(RTIIRFilter * rtif) {
    return IIRFilterBank_get_b_size(&rtif->ifb);
}
static inline double * RTIIRFilter_get_a(RTIIRFilter * rtif) {
    return IIRFilterBank_get_a(&rtif->ifb);
}
static inline size_t RTIIRFilter_get_a_size(RTIIRFilter * rtif) {
    return IIRFilterBank_get_a_size(&rtif->ifb);
}

// This converts the epsilon in the standard definition of the Chebyshev filters to the decibel ripples accepted by scipy and Matlab
// for Chebyshev Type I
static inline double passband_ripple_epsilon_to_db(double epsilon) {
    return 10.0 * log10(1 + pow(epsilon, 2));
}

// This converts the ripple from Matlab and scipy chebyshev filters to the epsilon used in the standard definitions
// for Chebyshev Type I
static inline double passband_ripple_db_to_epsilon(double db) {
    return sqrt(pow(10, db / 10.0) - 1);
}

// This converts the epsilon in the standard definition of the Chebyshev filters to the decibel ripples accepted by scipy and Matlab
// for Chebyshev Type I
static inline double stopband_ripple_epsilon_to_db(double epsilon) {
    return 10.0 * log10(pow(1/epsilon, 2) + 1);
}

// This converts the ripple from Matlab and scipy chebyshev filters to the epsilon used in the standard definitions
// for Chebyshev Type I
static inline double stopband_ripple_db_to_epsilon(double db) {
    return 1.0 / sqrt(pow(10, db / 10.0) - 1.0);
}

#endif // RTFILTER_H_