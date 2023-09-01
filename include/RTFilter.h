#ifndef RTFILTER_H
#define RTFILTER_H

#include <stddef.h>
#include <dsputils.h>

// configuration options
#ifndef DSP_MALLOC
    #define DSP_MALLOC malloc
#endif
#ifndef DSP_FREE
    #define DSP_FREE free
#endif
#ifndef DSP_CALLOC
    #define DSP_CALLOC calloc
#endif
#ifndef DSP_REALLOC
    #define DSP_REALLOC realloc
#endif

typedef struct RTFilter RTFilter;
typedef struct RTIIRFilter RTIIRFilter;
typedef struct RTFIRFilter RTFIRFilter;
typedef struct IIRFilterState IIRFilterState;
typedef struct IIRFilterBank IIRFilterBank;
typedef struct FilterBank FilterBank;

// this is really an interface. You should not be allocating (or you will leak) or instantiating (or you will screw it up) one directly
struct RTFilter {
    int (*update)(RTFilter * rtlf, double sample);
    int (*initialize)(RTFilter * rtlf, double sample);
    void (*del)(RTFilter * rtlf);
    double filtered_value;
    unsigned int flags;
    int initialized;
};

struct FilterBank {
    double * b;
    size_t nb;
};

struct IIRFilterBank {
    FilterBank fb;
    size_t na;
};

struct RTIIRFilter {
    RTFilter rtf;    // isa
    IIRFilterBank ifb;      // hasa, for second order sections execution, this is not a valid IIRFilter layout as it won't have a "FilterBank"
    double * state;    // state of real-time IIR filter
};

struct RTFIRFilter {
    RTFilter rtf;    // isa
    FilterBank fb;          // hasa
    double * state;    // state of real-time IIR filter
};

// define filter flags
#define FILTER_INFINITE 0x0
#define FILTER_FINITE 0x1
#define FILTER_EXECUTION_MASK 0x2
#define FILTER_TF 0x0
#define FILTER_SOS 0x2

#define FILTER_RESET -1
// DO NOT MAKE THIS > 0. Certain initialization may depend on this being <= 0
#define FILTER_INITIALIZED 0

EXPORT void RTFilter_reset(RTFilter * rtf);

EXPORT void RTFilter_init(RTFilter * rtf, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample), \
    int (*update)(RTFilter * rtf, double sample), \
    void (*del)(RTFilter * rtf));

EXPORT RTFilter * RTFilter_new(unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample), \
    int (*update)(RTFilter * rtf, double sample), \
    void (*del)(RTFilter * rtf));

// Note that this does not free the RTFilter itself as one should never actually be allocated. I might just delete this
EXPORT void RTFilter_del(RTFilter * rtf);
// note that this is very different from RTFilter_init
EXPORT void RTFilter_initialize(RTFilter * rtf, double sample);
EXPORT double RTFilter_update(RTFilter * rtf, double sample);

EXPORT void FilterBank_init(FilterBank * fb, double * b, size_t nb);
EXPORT FilterBank * FilterBank_new(double * b, size_t nb);
EXPORT void FilterBank_del(FilterBank * fb);

EXPORT void FilterBank_print(FilterBank * fb);


EXPORT void IIRFilterBank_init(IIRFilterBank * ifb, double * bank, size_t na, size_t nb);

// creates a filter bank with all values set to 0.0
EXPORT IIRFilterBank * IIRFilterBank_new_empty(size_t na, size_t nb);

// creates a copy of the two arrays
EXPORT IIRFilterBank * IIRFilterBank_new(double * a, size_t na, double * b, size_t nb);

// arguments must be this way so that na can be specified, otherwise passing an integer will simply convert to integer and not size_t
EXPORT IIRFilterBank * IIRFilterBank_new_coefs(size_t na, size_t nb, ...);

EXPORT void IIRFilterBank_del(IIRFilterBank * ifb);

EXPORT void IIRFilterBank_print(IIRFilterBank * ifb);

EXPORT void RTFIRFilter_del(RTFilter * rtf);

EXPORT void RTFIRFilter_init(RTFIRFilter * rtff, double * bank, double * state, size_t nb, \
    unsigned int flags, int (*initialize)(RTFilter * rtf, double sample));

EXPORT RTFIRFilter * RTFIRFilter_new_empty(size_t nb, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample));

EXPORT RTFIRFilter * RTFIRFilter_new(double * b, size_t nb, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample));

EXPORT int RTFIRFilter_update(RTFilter * rtf, double sample);

EXPORT int RTFIRFilter_stable_init(RTFilter * rtf, double sample);

// initialization function which does a partial moving average
EXPORT int RTFIRFilter_partial_init(RTFilter * rtf, double sample);

EXPORT void RTFIRFilter_print(RTFIRFilter * rtff);

EXPORT void RTIIRFilter_del(RTFilter * rtf);

EXPORT void RTIIRFilter_init(RTIIRFilter * rtif, double * bank, double * state, size_t na, size_t nb, \
    unsigned int flags, int (*initialize)(RTFilter * rtf, double sample));

EXPORT RTIIRFilter * RTIIRFilter_new_empty(size_t na, size_t nb, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample));

EXPORT RTIIRFilter * RTIIRFilter_new(double * a, size_t na, double * b, size_t nb, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample));

EXPORT int RTIIRFilter_update(RTFilter * rtf, double sample);

EXPORT int RTIIRFilter_stable_init(RTFilter * rtf, double sample);

EXPORT void RTIIRFilter_print(RTIIRFilter * rtif);

EXPORT int moving_average(RTFIRFilter * rtff, size_t window, \
    int (*initialize)(RTFilter * rtf, double sample));

/* static allocation
wl and wu define the cutoff frequences for 1/2 power attenuation.
wl = 0 --> low-pass filter
wu = 0 --> high-pass filter
wl < wu --> band-pass filter
wl > wu --> band-stop filter
wl == wu
wl and wu are in units of the Nyquist frequency 1/(2*T)
*/
EXPORT int butterworth(RTIIRFilter * rtif, size_t order, double wl, double wu, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample));

// This converts the epsilon in the standard definition of the Chebyshev filters to the decibel ripples accepted by scipy and Matlab
// for Chebyshev Type I
EXPORT inline double passband_ripple_epsilon_to_db(double epsilon) {
    return 10.0 * log10(1 + pow(epsilon, 2));
}

// This converts the ripple from Matlab and scipy chebyshev filters to the epsilon used in the standard definitions
// for Chebyshev Type I
EXPORT inline double passband_ripple_db_to_epsilon(double db) {
    return sqrt(pow(10, db / 10.0) - 1);
}

// This converts the epsilon in the standard definition of the Chebyshev filters to the decibel ripples accepted by scipy and Matlab
// for Chebyshev Type I
EXPORT inline double stopband_ripple_epsilon_to_db(double epsilon) {
    return 10.0 * log10(pow(1/epsilon, 2) + 1);
}

// This converts the ripple from Matlab and scipy chebyshev filters to the epsilon used in the standard definitions
// for Chebyshev Type I
EXPORT inline double stopband_ripple_db_to_epsilon(double db) {
    return 1.0 / sqrt(pow(10, db / 10.0) - 1.0);
}

EXPORT int chebyshev1(RTIIRFilter * rtif, size_t order, double ripple, double wl, double wu, \
    unsigned int flags, int (*initialize)(RTFilter * rtf, double sample));

EXPORT int chebyshev2(RTIIRFilter * rtif, size_t order, double ripple, double wl, double wu, 
    unsigned int flags, int (*initialize)(RTFilter * rtf, double sample));

#ifndef __STDC_NO_COMPLEX__
#include <complex.h>

#ifndef DEFAULT_COMPLEX_TOLERANCE
#define DEFAULT_COMPLEX_TOLERANCE 1e-7
#endif

EXPORT int pzg_to_RTIIRFilter(RTIIRFilter * rtif, double complex * poles, size_t N, \
    double complex * zeros, size_t M, double gain, double wl, double wu);

#endif // __STDC_NO_COMPLEX

#endif // RTFILTER_H