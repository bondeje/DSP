#ifndef RTFILTER_H_
#define RTFILTER_H_

#include <RTFilter.h>

EXPORT inline double * FilterBank_get_b(FilterBank * fb) {
    return fb->b;
}
EXPORT inline size_t FilterBank_get_b_size(FilterBank * fb) {
    return fb->nb;
}
EXPORT inline size_t FIR_order(FilterBank * fb) {
    return fb->nb;
}
EXPORT inline double * IIRFilterBank_get_b(IIRFilterBank * ifb) {
    return FilterBank_get_b(&ifb->fb);
}
EXPORT inline size_t IIRFilterBank_get_b_size(IIRFilterBank * ifb) {
    return FilterBank_get_b_size(&ifb->fb);
}
EXPORT inline double * IIRFilterBank_get_a(IIRFilterBank * ifb) {
    return IIRFilterBank_get_b(ifb) + IIRFilterBank_get_b_size(ifb);
}
EXPORT inline size_t IIRFilterBank_get_a_size(IIRFilterBank * ifb) {
    return ifb->na;
}
EXPORT inline double * RTFIRFilter_get_b(RTFIRFilter * rtff) {
    return rtff->fb.b;
}
EXPORT inline size_t RTFIRFilter_get_b_size(RTFIRFilter * rtff) {
    return rtff->fb.nb;
}
EXPORT inline double * RTIIRFilter_get_b(RTIIRFilter * rtif) {
    return rtif->ifb.fb.b;
}
EXPORT inline size_t RTIIRFilter_get_b_size(RTIIRFilter * rtif) {
    return IIRFilterBank_get_b_size(&rtif->ifb);
}
EXPORT inline double * RTIIRFilter_get_a(RTIIRFilter * rtif) {
    return IIRFilterBank_get_a(&rtif->ifb);
}
EXPORT inline size_t RTIIRFilter_get_a_size(RTIIRFilter * rtif) {
    return IIRFilterBank_get_a_size(&rtif->ifb);
}


#endif // RTFILTER_H_