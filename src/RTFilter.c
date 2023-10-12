// always compile with -O3
// gcc -pedantic -Wall -O3 -g -DDEBUG_MALLOC RTFilter.c polynomials.c Lpolys.c legendre.c chebyshev.c ../mallocs/debug_malloc.c -o f.exe

// high-level TODO:
//   reconsider struct layouts. Right now, IIRFilters and FIRFilters each get a pointer to update and del, which should really belong to a class object
//   solving this would require some method to solve the diamond problem and have actual inheritance


#ifdef DEBUG_MALLOC
    #include "../../mallocs/debug_malloc.h"
#endif
#include <stddef.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h> 
#include <time.h> // only for some debugging crap
#include <math.h> // requires NAN, isnan
#include <RTFilter_.h>
// move these defines to main or the dll header
#define POLY_FREE DSP_FREE
#define POLY_MALLOC DSP_MALLOC
#define POLY_CALLOC DSP_CALLOC
#define POLY_REALLOC DSP_REALLOC
#include <polynomials.h>
#include <legendre.h>
#include <chebyshev.h>
#include <specialpolys.h>
#ifndef M_PI
    #define M_PI		3.14159265358979323846
#endif
// Jeff's awesome (/s) debug_malloc. To use, add debug_malloc.h to include path, compile debug_malloc.c, and include compiler flags -DDEBUG_MALLOC and optionally -DDEBUG_MALLOC_VERBOSE. The first will just have a summary, the second will also tell you where and how much memory for each allocation and whether it was freed.



// define filter flags
#define FILTER_INFINITE 0x0
#define FILTER_FINITE 0x1
#define FILTER_EXECUTION_MASK 0x2
#define FILTER_TF 0x0
#define FILTER_SOS 0x2

#define FILTER_RESET -1
// DO NOT MAKE THIS > 0. Certain initialization may depend on this being <= 0
#define FILTER_INITIALIZED 0

// TODO: inline this
static inline bool flag_set(unsigned int flags, unsigned int value, unsigned int mask) {
    if (!mask) {
        mask = value;
    }
    return (flags & mask) == value;
}

void RTFilter_reset(RTFilter * rtf) {
    rtf->initialized = FILTER_RESET;
    rtf->filtered_value = 0.0;
}

void RTFilter_init(RTFilter * rtf, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample), \
    int (*update)(RTFilter * rtf, double sample), \
    void (*del)(RTFilter * rtf)) {
    rtf->update = update;
    rtf->initialize = initialize;
    rtf->del = del;
    rtf->flags = flags;
    
    RTFilter_reset(rtf);
}

RTFilter * RTFilter_new(unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample), \
    int (*update)(RTFilter * rtf, double sample), \
    void (*del)(RTFilter * rtf)) {

    RTFilter * rtf = (RTFilter *) DSP_MALLOC(sizeof(RTFilter));
    if (!rtf) {
        return NULL;
    }
    RTFilter_init(rtf, flags, initialize, update, del);
    return rtf;
}

// Note that this does not free the RTFilter itself as one should never actually be allocated. I might just delete this
void RTFilter_del(RTFilter * rtf) {
    if (!rtf) {
        return;
    }
    rtf->del(rtf);
}

// note that this is very different from RTFilter_init
void RTFilter_initialize(RTFilter * rtf, double sample) {
    rtf->initialize(rtf, sample);
}

double RTFilter_update(RTFilter * rtf, double sample) {
    //rtf->update(rtf, sample);
    //return rtf->filtered_value;
    if (rtf->initialized == FILTER_INITIALIZED) {
        rtf->update(rtf, sample);
    }
    else {
        RTFilter_initialize(rtf, sample);
    }
    return rtf->filtered_value;
}

// convenience function for doing a series of updates. There might be a more efficient way to do this with matrix math
int RTFilter_updaten(double * out, RTFilter * rtf, double * samples, size_t n) {
    if (!out || !rtf || !samples) {
        return 1;
    }
    for (size_t i = 0; i < n; i++) {
        out[i] = RTFilter_update(rtf, samples[i]);
    }
    return 0;
}

void FilterBank_print(FilterBank * fb) {
    for (size_t i = 0; i < fb->nb; i++) {
        printf("%.10f ", fb->b[i]);
    }
    printf("\n");
}

void FilterBank_init(FilterBank * fb, double * b, size_t nb) {
    if (!fb) {
        return;
    }
    fb->b = b;
    fb->nb = nb;
}

FilterBank * FilterBank_new(double * b, size_t nb) {
    FilterBank * fb = (FilterBank *) DSP_MALLOC(sizeof(FilterBank));
    if (!fb) {
        return NULL;
    }
    FilterBank_init(fb, b, nb);
    return fb;
}

void FilterBank_del(FilterBank * fb) {
    if (!fb) {
        return;
    }
    DSP_FREE(fb->b);
    fb->b = NULL;
    DSP_FREE(fb);
}

void IIRFilterBank_print(IIRFilterBank * ifb) {
    double * a = IIRFilterBank_get_a(ifb);
    size_t na = IIRFilterBank_get_a_size(ifb);
    printf("a coefficients (%zu):\n", na);
    for (size_t i = 0; i < na; i++) {
        printf("%.10f ", a[i]);
    }
    printf("\nb coefficients (%zu):\n", FilterBank_get_b_size(&ifb->fb));
    FilterBank_print(&ifb->fb);
}

void IIRFilterBank_init(IIRFilterBank * ifb, double * bank, size_t na, size_t nb) {
    if (!ifb) {
        return;
    }
    FilterBank_init((FilterBank *) ifb, bank, nb);
    ifb->na = na;
}

// creates a filter bank with all values set to 0.0
IIRFilterBank * IIRFilterBank_new_empty(size_t na, size_t nb) {
    IIRFilterBank * ifb = DSP_MALLOC(sizeof(IIRFilterBank));
    if (!ifb) {
        return NULL;
    }
    double * bank = (double *) DSP_MALLOC(sizeof(double) * (na + nb));
    if (!bank) {
        DSP_FREE(ifb);
        return NULL;
    }
    for (size_t i = 0; i < na + nb; i++) {
        bank[i] = 0.0;
    }
    IIRFilterBank_init(ifb, bank, na, nb);
    return ifb;
}

// creates a copy of the two arrays
IIRFilterBank * IIRFilterBank_new(double * a, size_t na, double * b, size_t nb) {
    IIRFilterBank * ifb = IIRFilterBank_new_empty(na, nb);
    if (!ifb) {
        return NULL;
    }
    double * b_ = IIRFilterBank_get_b(ifb);
    double * a_ = IIRFilterBank_get_a(ifb);
    for (size_t i = 0; i < na; i++) {
        a_[i] = a[i];   
    }
    for (size_t i = 0; i < nb; i++) {
        b_[i] = b[i];   
    }
    return ifb;
}

// arguments must be this way so that na can be specified, otherwise passing an integer will simply convert to integer and not size_t
IIRFilterBank * IIRFilterBank_new_coefs(size_t na, size_t nb, ...) {
    if (!na || !nb) {
        return NULL;
    }
    IIRFilterBank * ifb = IIRFilterBank_new_empty(na, nb);
    if (!ifb) {
        return NULL;
    }
    double * b = IIRFilterBank_get_b(ifb);
    double * a = IIRFilterBank_get_a(ifb);
    va_list args;
    va_start(args, nb);
    for (size_t i = 0; i < na; i++) {
        a[i] = va_arg(args, double);
    }
    for (size_t i = 0; i < nb; i++) {
        b[i] = va_arg(args, double);
    }
    va_end(args);
    return ifb;
}

void IIRFilterBank_del(IIRFilterBank * ifb) {
    FilterBank_del(&ifb->fb);
}

int RTFIRFilter_update(RTFilter * rtf, double sample) {
    RTFIRFilter * rtff = (RTFIRFilter *)rtf;
    size_t nb =  rtff->fb.nb;
    double * b = rtff->fb.b;

    rtf->filtered_value = sample*b[0] + rtff->state[0];
    size_t i = 0;
    while (i < nb-1) {
        rtff->state[i] = rtff->state[i+1] + b[i+1]*sample;
        i++;
    }
    return 0;
}

// TODO: this might have to get fixed
int RTFIRFilter_stable_init(RTFilter * rtf, double sample) {
    RTFIRFilter * rtff = (RTFIRFilter *)rtf;
    double * b = rtff->fb.b;

    size_t i = rtff->fb.nb-1; // iir order-1
    rtff->state[i] = 0.0;
    double csv = 0.0;
    while (i) { // decrement after check
        csv += b[i];
        i--;
        rtff->state[i] = csv*sample;
    }
    rtf->filtered_value = sample;
    rtf->initialized = FILTER_INITIALIZED;
    return 0;
}

void RTFIRFilter_print(RTFIRFilter * rtff) {
    printf("FIR order %zu, coefficients:", FilterBank_get_b_size(&rtff->fb));
    FilterBank_print(&rtff->fb);
}

// initialization function which does a partial moving average
int RTFIRFilter_partial_init(RTFilter * rtf, double sample) {
    RTFIRFilter * rtff = (RTFIRFilter *)rtf;
    double * b = rtff->fb.b;

    size_t order = rtff->fb.nb; // iir order

    if (rtf->initialized == FILTER_RESET) {
        rtff->state[order-1] = 0.0;
        rtf->filtered_value = sample;
        rtf->initialized = order - 1;
    } else {
        rtf->initialized--;
        rtf->filtered_value = order*(b[0]*sample + rtff->state[0])/(order-rtf->initialized);   
    }
    size_t i = 0;
    while (i < order - 1) {
        rtff->state[i] = rtff->state[i+1] + b[i+1]*sample;
        i++;
    }
    
    if (rtf->initialized == 0) {
        rtf->initialized = FILTER_INITIALIZED;
    }
    return 0;
}

void RTFIRFilter_del(RTFilter * rtf) {
    if (!rtf) {
        return;
    }
    RTFIRFilter * rtff = (RTFIRFilter *)rtf;
    DSP_FREE(rtff->fb.b);
    rtff->fb.b = NULL;
    DSP_FREE(rtff->state);
    rtff->state = NULL;
    DSP_FREE(rtf);
}

void RTFIRFilter_init(RTFIRFilter * rtff, double * bank, double * state, size_t nb, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample)) {
    RTFilter_init(&rtff->rtf, flags, initialize, RTFIRFilter_update, RTFIRFilter_del);
    FilterBank_init(&rtff->fb, bank, nb);
    rtff->state = state;
    for (size_t i = 0; i < FIR_order(&rtff->fb); i++) {
        rtff->state[i] = 0.0;
    }
}

// TODO: finish this implementation
RTFIRFilter * RTFIRFilter_new_empty(size_t nb, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample)) {
    if (!nb) {
        return NULL;
    }
    if (!initialize) {
        initialize = RTFIRFilter_stable_init;
    }
    RTFIRFilter * rtff = (RTFIRFilter *) DSP_MALLOC(sizeof(RTFIRFilter));
    if (!rtff) {
        return NULL;
    }
    // TODO: this is actually kind of ugly, RTFIRFilter_init iterates by calculating the order, but here, we don't have the order because the FilterBank has not been initialized. Find a way around
    double * state = (double *) DSP_MALLOC(sizeof(double)*nb); 
    if (!state) {
        DSP_FREE(rtff);
        return NULL;
    }
    // initialization of rtff->state done by RTFIRFilter_init

    double * bank = (double *) DSP_MALLOC(sizeof(double) * nb);
    if (!bank) {
        DSP_FREE(state);
        DSP_FREE(rtff);
        return NULL;
    }
    for (size_t i = 0; i < nb; i++) {
        bank[i] = 0.0;
    }

    RTFIRFilter_init(rtff, bank, state, nb, flags, initialize);
    
    return (RTFIRFilter*)rtff;
}

// this should copy the filter bank b and not just use it similar to IIR filter
RTFIRFilter * RTFIRFilter_new(double * b, size_t nb, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample)) {
    if (!nb) {
        return NULL;
    }
    if (!initialize) {
        initialize = RTFIRFilter_stable_init;
    }
    RTFIRFilter * rtff = RTFIRFilter_new_empty(nb, flags, initialize);
    if (!rtff) {
        return NULL;
    }
    // RTFIRFilter_new_empty initializes

    // copy b into filterbank
    double * bank = RTFIRFilter_get_b(rtff);
    for (size_t i = 0; i < nb; i++) {
        bank[i] = b[i];
    }
    
    return rtff;
}

// TODO: implement, but it needs to have matrix inversion to be "correct"
RTFIRFilter * MatchedFilter_new(double * signal, size_t nb, double * covariance, size_t dim) {
    return NULL;
}

void RTIIRFilter_get_all(RTFilter * rtf, double ** a, size_t * na, double ** b, size_t * nb) {
    IIRFilterBank * ifb = &((RTIIRFilter *) rtf)->ifb;
    *na = ifb->na;
    *nb = ifb->fb.nb;
    *a = ifb->fb.b + *nb;
    *b = ifb->fb.b;
}

//int test = 1;

int RTIIRFilter_update(RTFilter * rtf, double sample) {
    RTIIRFilter * rtif = (RTIIRFilter *) rtf;
    size_t na = 0, nb = 0;
    double *a = NULL, *b = NULL;
    RTIIRFilter_get_all(rtf, &a, &na, &b, &nb);
    size_t N = ((nb >= na) ? nb : na) - 1;
    /*
    if (!test) {
        printf("na = %zu, nb = %zu\na = ", na, nb);
        for (size_t i = 0; i < na; i++) {
            printf("%.8f ", a[i]);
        }
        printf("\nb = ");
        for (size_t i = 0; i < nb; i++) {
            printf("%.8f ", b[i]);
        }
        printf("\n");
        test = 1;
    }
    */

    double new_filt = (sample*b[0] + rtif->state[0])/a[0];

    double state_value = 0.0;
    size_t i = 0;
    //printf("state ingesting new sample\n");
    while (i < N) {
        //printf("%.5f ", rtif->state[i]);
        state_value = 0.0;
        if (i < nb - 1) {
            state_value += b[i+1] * sample;
        }
        if (i < na - 1) {
            state_value -= a[i+1] * new_filt;
        }
        rtif->state[i] = rtif->state[i+1] + state_value;
        i++;
    }
    //printf("%.5f\n", rtif->state[i]);

    rtf->filtered_value = new_filt;
    
    return 0;
}

int RTIIRFilter_stable_init(RTFilter * rtf, double sample) {
    RTIIRFilter * rtif = (RTIIRFilter *) rtf;
    size_t na = 0, nb = 0;
    double *a = NULL, *b = NULL;
    RTIIRFilter_get_all(rtf, &a, &na, &b, &nb);

    size_t i = ((nb >= na) ? nb : na) - 1; // iir order
    double num = 0.0;
    double den = 0.0;
    for (size_t j = 0; j <= 1; j++) {
        if (j < nb) {
            num += b[j];
        }
        if (j < na) {
            den += a[j];
        }
    }
    double d = num / den; // setting d = 1.0 goes back to the old method
    
    rtif->state[i] = 0.0;
    double csv = 0.0;
    while (i) { // decrement after check
        if (i < nb) {
            csv += b[i];
        }
        if (i < na) {
            csv -= a[i] * d;
        }
        i--;
        rtif->state[i] = csv*sample;
    }
    rtf->filtered_value = d * sample;
    rtf->initialized = FILTER_INITIALIZED;
    return 0;
}

void RTIIRFilter_print(RTIIRFilter * rtif) {
    size_t na = IIRFilterBank_get_a_size(&rtif->ifb);
    size_t nb = IIRFilterBank_get_b_size(&rtif->ifb);
    na = (na >= nb) ? na : nb;
    printf("IIR order %zu:\n", na);
    IIRFilterBank_print(&(rtif->ifb));
}

void RTIIRFilter_del(RTFilter * rtf) {
    if (!rtf) {
        return;
    }
    RTIIRFilter * rtif = (RTIIRFilter *)rtf;
    DSP_FREE(rtif->ifb.fb.b);
    rtif->ifb.fb.b = NULL;
    DSP_FREE(rtif->state);
    rtif->state = NULL;
    DSP_FREE(rtif);
}

void RTIIRFilter_init(RTIIRFilter * rtif, double * bank, double * state, size_t na, size_t nb, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample)) {
    RTFilter_init((RTFilter *)rtif, flags, initialize, RTIIRFilter_update, RTIIRFilter_del);
    IIRFilterBank_init(&rtif->ifb, bank, na, nb);
    rtif->state = state;
    for (size_t i = 0; i < ((nb >= na) ? nb : na); i++) {
        rtif->state[i] = 0.0;
    }
}

RTIIRFilter * RTIIRFilter_new_empty(size_t na, size_t nb, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample)) {
    //printf("in RTIIRFilter_new\n");
    if (!na) {
        return NULL;
    }
    if (!initialize) {
        initialize = RTIIRFilter_stable_init;
    }
    //printf("starting allocations\n");
    RTIIRFilter * rtif = (RTIIRFilter *) DSP_MALLOC(sizeof(RTIIRFilter));
    if (!rtif) {
        return NULL;
    }
    //printf("allocated RTIIRFilter\n");
    double * bank = (double *) DSP_MALLOC(sizeof(double) * (na + nb));
    if (!bank) {
        DSP_FREE(rtif);
        return NULL;
    }
    //printf("allocated bank\n");
    
    size_t Nstate = (na > nb ) ? na : nb;
    // TODO: this is actually kind of ugly, RTFIRFilter_init iterates by calculating the order, but here, we don't have the order because the FilterBank has not been initialized. Find a way around
    double * state = (double *) DSP_MALLOC(sizeof(double) * Nstate);
    if (!state) {
        DSP_FREE(bank);
        DSP_FREE(rtif);
        return NULL;
    }
    //printf("allocations done\n");
    for (size_t i = 0; i < na + nb; i++) {
        bank[i] = 0.0;
    }

    RTIIRFilter_init(rtif, bank, state, na, nb, flags, initialize);
    
    return rtif;
}

RTIIRFilter * RTIIRFilter_new(double * a, size_t na, double * b, size_t nb, unsigned int flags, \
    int (*initialize)(RTFilter * rtf, double sample)) {
    //printf("in RTIIRFilter_new\n");
    if (!na) {
        return NULL;
    }
    RTIIRFilter * rtif = RTIIRFilter_new_empty(na, nb, flags, initialize);
    if (!rtif) {
        return NULL;
    }
    double *a_ = RTIIRFilter_get_a(rtif);
    double *b_ = RTIIRFilter_get_b(rtif);
    //printf("allocations done\n");
    for (size_t i = 0; i < nb; i++) {
        b_[i] = b[i];
    }
    for (size_t i = 0; i < na; i++) {
        a_[i] = a[i];
    }    
    return rtif;
}

int moving_average(RTFIRFilter * rtff, size_t window, int (*initialize)(RTFilter * rtf, double sample)) {
    if (!rtff || !window) {
        return 1;
    }
    if (rtff->fb.nb < window)  {
        return 1;
    }
    if (!initialize) {
        initialize = RTFIRFilter_partial_init;
    }
    double * weights = RTFIRFilter_get_b(rtff);
    for (size_t i = 0; i < window; i++) {
        weights[i] = 1.0/window;
    }
    RTFIRFilter_init(rtff, weights, rtff->state, window, 0, initialize);
    return 0;
}

int dp_lp2lp(IIRFilterBank * ifb, size_t L, double wp, double wlp) { 
    if (wp == wlp) { // do nothing if the source and destination warping frequencies are the same
        return 0;
    }
    size_t N1 = ifb->na;
    size_t N2 = N1 + (L - N1) / 2;
    size_t M1 = ((FilterBank *) ifb)->nb;
    size_t M2 = M1 + (L - M1) / 2;

    ifb->na = L + 1;
    ((FilterBank *) ifb)->nb = L + 1;
    double * pgain = IIRFilterBank_get_b(ifb);
    double * b = pgain + 1;
    double * a = IIRFilterBank_get_a(ifb) + 1;

    double alpha = sin(M_PI * (wp - wlp) / 2.0) / sin(M_PI * (wp + wlp) / 2.0); // really -cos(M_PI / 2.0 * (wl + w0)), but we will use wl == w0
    //printf("alpha = %.4f\n", alpha);
    double alpha2 = pow(alpha, 2);
    size_t i = 0;
    //printf("im %zu, M1 %zu, M2 %zu, in %zu, N1 %zu, N2 %zu\n", im, M1, M2, in, N1, N2);
    while (i <= ((N2 >= M2) ? N2 : M2)) {
        double num = 1.0, den = 1.0;
        if (i < M1) {
            //printf("simple zero %.4f\n", *b_);
            num = 1 - alpha * (*b);
            *b = (*b - alpha) / num;
            b++;
        } else if (i < M2) {
            double B1 = *b;
            double B2 = *(b + 1);
            //printf("complex zero: %.4f, %.4f\n", B1, B2);
            num = 1 - alpha * B1 + alpha2 * B2;
            *b = -(2 * alpha - B1 * (1 + alpha2) + 2 * alpha * B2) / num;
            b++;
            *b = (alpha2 - alpha * B1 + B2) / num;
            //printf("transformed to: %.4f, %.4f\n", *(b_-1), *b_);
            b++;
        }

        if (i < N1) {
            //printf("simple pole %.4f\n", *a_);
            den = 1 - alpha * (*a);
            *a = (*a - alpha) / den;
            a++;
        } else if (i < N2) {
            double A1 = *a;
            double A2 = *(a + 1);
            //printf("complex zero: %.4f, %.4f\n", A1, A2);
            den = 1 - alpha * A1 + alpha2 * A2;
            *a = -(2 * alpha - A1 * (1 + alpha2) + 2 * alpha * A2) / den;
            a++;
            *a = (alpha2 - alpha * A1 + A2) / den;
            //printf("transformed to: %.4f, %.4f\n", *(a_-1), *a_);
            a++;
        }
        (*pgain) *= num / den;
        i++;
    }
    ifb->na = N1;
    ((FilterBank *) ifb)->nb = M1;
    return 0;
}

int dp_lp2hp(IIRFilterBank * ifb, size_t L, double wp, double whp) { 
    size_t N1 = ifb->na;
    size_t N2 = N1 + (L - N1) / 2;
    size_t M1 = ((FilterBank *) ifb)->nb;
    size_t M2 = M1 + (L - M1) / 2;

    ifb->na = L + 1;
    ((FilterBank *) ifb)->nb = L + 1;
    double * pgain = IIRFilterBank_get_b(ifb);
    double * b = pgain + 1;
    double * a = IIRFilterBank_get_a(ifb) + 1;

    double alpha = - cos(M_PI * (wp + whp) / 2.0) / cos(M_PI * (wp - whp) / 2.0); // really -cos(M_PI / 2.0 * (wl + w0)), but we will use wl == w0
    //printf("alpha = %.4f\n", alpha);
    double alpha2 = pow(alpha, 2);
    //double * b_ = b+1;
    //double * a_ = a+1;
    size_t i = 0;
    //printf("im %zu, M1 %zu, M2 %zu, in %zu, N1 %zu, N2 %zu\n", im, M1, M2, in, N1, N2);
    while (i <= ((N2 >= M2) ? N2 : M2)) {
        double num = 1.0, den = 1.0;
        if (i < M1) {
            //printf("simple zero %.4f\n", *b_);
            num = 1 - alpha * (*b);
            *b = (alpha - *b) / num;
            b++;
        } else if (i < M2) {
            double B1 = *b;
            double B2 = *(b + 1);
            //printf("complex zero: %.4f, %.4f\n", B1, B2);
            num = 1 - alpha * B1 + alpha2 * B2;
            *b = (2 * alpha - B1 * (1 + alpha2) + 2 * alpha * B2) / num;
            b++;
            *b = (alpha2 - alpha * B1 + B2) / num;
            //printf("transformed to: %.4f, %.4f\n", *(b_-1), *b_);
            b++;
        }

        if (i < N1) {
            //printf("simple pole %.4f\n", *a_);
            den = 1 - alpha * (*a);
            *a = (alpha - *a) / den;
            a++;
        } else if (i < N2) {
            double A1 = *a;
            double A2 = *(a + 1);
            //printf("complex zero: %.4f, %.4f\n", A1, A2);
            den = 1 - alpha * A1 + alpha2 * A2;
            *a = (2 * alpha - A1 * (1 + alpha2) + 2 * alpha * A2) / den;
            a++;
            *a = (alpha2 - alpha * A1 + A2) / den;
            //printf("transformed to: %.4f, %.4f\n", *(a_-1), *a_);
            a++;
        }
        (*pgain) *= num / den;
        i++;
    }
    ifb->na = N1;
    ((FilterBank *) ifb)->nb = M1;
    return 0;
}

int dp_lp2bp(IIRFilterBank * ifb, size_t L, double wp, double wl, double wu) { 
    size_t N1 = ifb->na;
    size_t N2 = N1 + (L - N1) / 2;
    size_t M1 = ((FilterBank *) ifb)->nb;
    size_t M2 = M1 + (L - M1) / 2;

    // currently 0, 1, 2, 3, ..., L, L+1, ..., 2 * L +1 are occupied by b [0, L + 1) and a [L + 1, 2 * L + 2)
    // need to map this to 0, 1, 3, 5, 7, ..., 2 * L - 1, 2 * L + 1, 2 * L + 2, 2 * L + 4, ..., 4 * L so that it spans the new B [0, 2 * L + 1), [2 * L + 1, 4 * L + 2)
    // set to 0.0 the values between those listed above
    //printf("L = %zu, N1 = %zu, N2 = %zu, M1 = %zu, M2 = %zu\n", L, N1, N2, M1, M2);
    size_t i = 4 * L + 2;
    size_t j = 2 * L + 2;
    double * bank = ifb->fb.b;
    ifb->na = L + 1;
    ((FilterBank *) ifb)->nb = L + 1;
    //printf("before\n");
    //IIRFilterBank_print((FilterBank *) ifb);
    while (j > L + 2) {
        i -= 2;
        j--;
        //printf("j(%zu) --> i(%zu)\n", j, i);
        *(bank + i) = *(bank + j);
        *(bank + i - 1) = 0.0;
    }
    i--;
    j--;
    //printf("j(%zu) --> i(%zu)\n", j, i);
    *(bank + i) = *(bank + j);
    *(bank + i - 1) = 0.0;
    while (j > 2) {
        i -= 2;
        j--;
        //printf("j(%zu) --> i(%zu)\n", j, i);
        *(bank + i) = *(bank + j);
        *(bank + i - 1) = 0.0;
    }
    //printf("i = %zu, and j = %zu\n", i, j);

    ifb->na = 2 * L + 1;
    ((FilterBank *) ifb)->nb = 2 * L + 1;
    //printf("after\n");
    //IIRFilterBank_print((FilterBank *) ifb);
    double * pgain = IIRFilterBank_get_b(ifb);
    double * b = pgain + 1;
    double * a = IIRFilterBank_get_a(ifb) + 1;

    // here...need to make sure it grabs the correct parts of the filter bank.
    double chi = tan(M_PI * wp / 2.0) / tan(M_PI * (wu - wl) / 2.0);
    double alpha = - 2 * cos(M_PI * (wu + wl) / 2.0) / cos(M_PI * (wu - wl) / 2.0) * chi / (chi + 1.0);
    //printf("alpha = %.4f\n", alpha);
    double alpha2 = pow(alpha, 2);
    double beta = (chi - 1.0) / (chi + 1.0);
    double beta2 = pow(beta, 2);
    //double * b_ = b+1;
    //double * a_ = a+1;

    double coefs[] = {  2 * alpha,          \
                        alpha * (1 + beta), \
                        2 * alpha * beta,   \
                        (alpha2 + 2 * beta),\
                        (alpha2 + 1 + beta2)\
                     };

    i = 0;
    //printf("im %zu, M1 %zu, M2 %zu, in %zu, N1 %zu, N2 %zu\n", im, M1, M2, in, N1, N2);
    while (i <= ((N2 >= M2) ? N2 : M2)) {
        double num = 1.0, den = 1.0;
        if (i < M1) {
            //printf("simple zero %.4f\n", *b);
            double B0 = *b;
            num = 1 - beta * B0;
            *b = alpha * (1 - B0) / num;
            b++;
            *b = (beta - B0) / num;
            b++;
        } else if (i < M2) {
            double B1 = *b;
            double B2 = *(b + 2);
            //printf("complex zero: %.4f, %.4f\n", B1, B2);
            num = 1 - beta * B1 + beta2 * B2;
            //*b = (2 * alpha - B1 * alpha * (1 + beta) + 2 * alpha * beta * B2) / num;
            *b = (coefs[0] - B1 * coefs[1] + coefs[2] * B2) / num;
            b++;
            //*b = ((alpha2 + 2 * beta) - B1 * (alpha2 + 1 + beta2) + B2 * (alpha2 + 2 * beta)) / num;
            *b = (coefs[3] - B1 * coefs[4] + B2 * coefs[3]) / num;
            b++;
            //*b = (2 * alpha * beta - B1 * alpha * (1 + beta) + 2 * alpha * B2) / num;
            *b = (coefs[2] - B1 * coefs[1] + coefs[0] * B2) / num;
            b++;
            *b = (beta2 - B1 * beta + B2) / num;
            b++;
        }

        if (i < N1) {
            //printf("simple zero %.4f\n", *b_);
            double A0 = *a;
            den = 1 - beta * A0;
            *a = alpha * (1 - A0) / den;
            a++;
            *a = (beta - A0) / den;
            a++;
        } else if (i < N2) {
            double A1 = *a;
            double A2 = *(a + 2);
            //printf("complex zero: %.4f, %.4f\n", B1, B2);
            den = 1 - beta * A1 + beta2 * A2;
            //*a = (2 * alpha - A1 * alpha * (1 + beta) + 2 * alpha * beta * A2) / den;
            *a = (coefs[0] - A1 * coefs[1] + coefs[2] * A2) / den;
            a++;
            //*a = ((alpha2 + 2 * beta) - A1 * (alpha2 + 1 + beta2) + A2 * (alpha2 + 2 * beta)) / den;
            *a = (coefs[3] - A1 * coefs[4] + A2 * coefs[3]) / den;
            a++;
            //*a = (2 * alpha * beta - A1 * alpha * (1 + beta) + 2 * alpha * A2) / den;
            *a = (coefs[2] - A1 * coefs[1] + coefs[0] * A2) / den;
            a++;
            *a = (beta2 - A1 * beta + A2) / den;
            a++;
        }
        (*pgain) *= num / den;
        i++;
    }
    

    // TODO: should this double?
    ifb->na = 2 * N1;
    ((FilterBank *) ifb)->nb = 2 * M1;
    return 0;
}

int dp_lp2bs(IIRFilterBank * ifb, size_t L, double wp, double wl, double wu) { 
    size_t N1 = ifb->na;
    size_t N2 = N1 + (L - N1) / 2;
    size_t M1 = ((FilterBank *) ifb)->nb;
    size_t M2 = M1 + (L - M1) / 2;

    // currently 0, 1, 2, 3, ..., L, L+1, ..., 2 * L +1 are occupied by b [0, L + 1) and a [L + 1, 2 * L + 2)
    // need to map this to 0, 1, 3, 5, 7, ..., 2 * L - 1, 2 * L + 1, 2 * L + 2, 2 * L + 4, ..., 4 * L so that it spans the new B [0, 2 * L + 1), [2 * L + 1, 4 * L + 2)
    // set to 0.0 the values between those listed above
    printf("L = %zu, N1 = %zu, N2 = %zu, M1 = %zu, M2 = %zu\n", L, N1, N2, M1, M2);
    size_t i = 4 * L + 2;
    size_t j = 2 * L + 2;
    double * bank = ifb->fb.b;
    ifb->na = L + 1;
    ((FilterBank *) ifb)->nb = L + 1;
    //printf("before\n");
    //IIRFilterBank_print((FilterBank *) ifb);
    while (j > L + 2) {
        i -= 2;
        j--;
        //printf("j(%zu) --> i(%zu)\n", j, i);
        *(bank + i) = *(bank + j);
        *(bank + i - 1) = 0.0;
    }
    i--;
    j--;
    //printf("j(%zu) --> i(%zu)\n", j, i);
    *(bank + i) = *(bank + j);
    *(bank + i - 1) = 0.0;
    while (j > 2) {
        i -= 2;
        j--;
        //printf("j(%zu) --> i(%zu)\n", j, i);
        *(bank + i) = *(bank + j);
        *(bank + i - 1) = 0.0;
    }
    //printf("i = %zu, and j = %zu\n", i, j);

    ifb->na = 2 * L + 1;
    ((FilterBank *) ifb)->nb = 2 * L + 1;
    //printf("after\n");
    //IIRFilterBank_print((FilterBank *) ifb);
    double * pgain = IIRFilterBank_get_b(ifb);
    double * b = pgain + 1;
    double * a = IIRFilterBank_get_a(ifb) + 1;

    // here...need to make sure it grabs the correct parts of the filter bank.
    double chi = tan(M_PI * wp / 2.0) * tan(M_PI * (wu - wl) / 2.0);
    double alpha = - 2 * cos(M_PI * (wu + wl) / 2.0) / cos(M_PI * (wu - wl) / 2.0) / (chi + 1.0);
    //printf("alpha = %.4f\n", alpha);
    double alpha2 = pow(alpha, 2);
    double beta = (1.0 - chi) / (chi + 1.0);
    double beta2 = pow(beta, 2);

    double coefs[] = {  2 * alpha,          \
                        alpha * (1 + beta), \
                        2 * alpha * beta,   \
                        (alpha2 + 2 * beta),\
                        (alpha2 + 1 + beta2)\
                     };

    //double * b_ = b+1;
    //double * a_ = a+1;
    i = 0;
    //printf("im %zu, M1 %zu, M2 %zu, in %zu, N1 %zu, N2 %zu\n", im, M1, M2, in, N1, N2);
    while (i <= ((N2 >= M2) ? N2 : M2)) {
        double num = 1.0, den = 1.0;
        if (i < M1) {
            //printf("simple zero %.4f\n", *b);
            double B0 = -(*b);
            num = 1 - beta * B0;
            *b = alpha * (1 - B0) / num;
            b++;
            *b = (beta - B0) / num;
            b++;
        } else if (i < M2) {
            double B1 = -(*b);
            double B2 = *(b + 2);
            //printf("complex zero: %.4f, %.4f\n", B1, B2);
            num = 1 - beta * B1 + beta2 * B2;
            //*b = (2 * alpha - B1 * alpha * (1 + beta) + 2 * alpha * beta * B2) / num;
            *b = (coefs[0] - B1 * coefs[1] + coefs[2] * B2) / num;
            b++;
            //*b = ((alpha2 + 2 * beta) - B1 * (alpha2 + 1 + beta2) + B2 * (alpha2 + 2 * beta)) / num;
            *b = (coefs[3] - B1 * coefs[4] + B2 * coefs[3]) / num;
            b++;
            //*b = (2 * alpha * beta - B1 * alpha * (1 + beta) + 2 * alpha * B2) / num;
            *b = (coefs[2] - B1 * coefs[1] + coefs[0] * B2) / num;
            b++;
            *b = (beta2 - B1 * beta + B2) / num;
            b++;
        }

        if (i < N1) {
            //printf("simple zero %.4f\n", *b_);
            double A0 = -(*a);
            den = 1 - beta * A0;
            *a = alpha * (1 - A0) / den;
            a++;
            *a = (beta - A0) / den;
            a++;
        } else if (i < N2) {
            double A1 = -(*a);
            double A2 = *(a + 2);
            //printf("complex zero: %.4f, %.4f\n", B1, B2);
            den = 1 - beta * A1 + beta2 * A2;
            //*a = (2 * alpha - A1 * alpha * (1 + beta) + 2 * alpha * beta * A2) / den;
            *a = (coefs[0] - A1 * coefs[1] + coefs[2] * A2) / den;
            a++;
            //*a = ((alpha2 + 2 * beta) - A1 * (alpha2 + 1 + beta2) + A2 * (alpha2 + 2 * beta)) / den;
            *a = (coefs[3] - A1 * coefs[4] + A2 * coefs[3]) / den;
            a++;
            //*a = (2 * alpha * beta - A1 * alpha * (1 + beta) + 2 * alpha * A2) / den;
            *a = (coefs[2] - A1 * coefs[1] + coefs[0] * A2) / den;
            a++;
            *a = (beta2 - A1 * beta + A2) / den;
            a++;
        }
        (*pgain) *= num / den;
        i++;
    }
    

    // TODO: should this double?
    ifb->na = 2 * N1;
    ((FilterBank *) ifb)->nb = 2 * M1;
    return 0;
}

// L+1 = number of coefficients in both the numerator and denominator IIRFilterBank components
// wl == wu, fails as invalid
// wl <= 0, performs low-pass transformation, cutoff at wu, since prototype is already built on low-pass, really just does IIRFilterBank transformation
// wu <= 0, performs high-pass transformation, cutoff at wl
// wu > wl, performs band-pass transformation, cutoffs at wl, wu
// wl < wu, performs band-stop transformation, 
// nans in the digital prototype are interpreted as 1
int digital_prototype_to_IIRFilterBank(IIRFilterBank * ifb, size_t L, double wp, double wl, double wu) {
    //printf("in digital_prototype_to_IIRFilterBank\n");
    if (!ifb) {
        return 1;
    }
    if (wl == wu) {
        return 1;
    }
    
    size_t poly_exp = 1;
    size_t na = L + 1;
    size_t N1 = ifb->na;
    size_t N2 = N1 + (na - N1 - 1) / 2;
    size_t nb = L + 1;
    size_t M1 = ((FilterBank *) ifb)->nb;
    size_t M2 = M1 + (nb - M1 - 1) / 2;

    //printf("before transformation\n");
    //IIRFilterBank_print((FilterBank *) ifb);
    if (wl <= 0.0) {
        dp_lp2lp(ifb, L, wp, wu);
    } else if (wu <= 0.0) { // high-pass
        dp_lp2hp(ifb, L, wp, wl);
    } else if (wu > wl) { // band-pass
        dp_lp2bp(ifb, L, wp, wl, wu);
        poly_exp = 2;
        na = 2 * L + 1;
        nb = 2 * L + 1;
    } else if (wu < wl) { // band-stop
        dp_lp2bs(ifb, L, wp, wu, wl);
        poly_exp = 2;
        na = 2 * L + 1;
        nb = 2 * L + 1;
    }// TODO: other transformations. For band-stop and band-pass, set poly_exp = 2;
    //printf("gain after %.8f\n", *pgain);
    //printf("after transformation\n");
    //IIRFilterBank_print((FilterBank *) ifb);

    ifb->na = na;
    ((FilterBank *) ifb)->nb = nb;
    double * pgain = IIRFilterBank_get_b(ifb);
    double * b = pgain;
    double * a = IIRFilterBank_get_a(ifb);
    //printf("gain factor %.16f\n", *pgain);
    if (poly_exp == 1) {
        Polynomial p2;
        double multipliers[3] = {1.0, 1.0, 0.0};
        Polynomial_init(&p2, 1, multipliers, 3);
        
        Polynomial pb;
        Polynomial_init(&pb, 0, b, nb);
        Polynomial pa;
        Polynomial_init(&pa, 0, a, na);

        b++;
        p2.order = 1;
        //printf("p2 order %zu, pb order %zu, p2 ncoefs %zu, pb ncoefs %zu\n", p2.order, pb. order, p2.ncoefs, pb.ncoefs);
        for (size_t i = 0; i < M1; i++) {
            p2.coefs[1] = *b;
            b++;
            Polynomial_mul(&pb, &p2);
        }
        p2.order = 2;
        for (size_t i = M1; i < M2; i++) {
            p2.coefs[1] = *b;
            b++;
            p2.coefs[2] = *b;
            b++;
            Polynomial_mul(&pb, &p2);
        }
        p2.coefs[2] = 0.0;

        a++;
        p2.order = 1;
        //printf("p2 order %zu, pa order %zu, p2 ncoefs %zu, pa ncoefs %zu\n", p2.order, pa. order, p2.ncoefs, pa.ncoefs);
        for (size_t i = 0; i < N1; i++) {
            p2.coefs[1] = *a;
            a++;
            Polynomial_mul(&pa, &p2);
        }
        p2.order = 2;
        //printf("performing order 2 polynomial multiplications\n");
        for (size_t i = N1; i < N2; i++) {
            p2.coefs[1] = *a;
            a++;
            p2.coefs[2] = *a;
            a++;
            Polynomial_mul(&pa, &p2);
        }
        p2.coefs[2] = 0.0;
    } else if (poly_exp == 2) {
        Polynomial p4;
        double multipliers[5] = {1.0, 1.0, 1.0, 0.0, 0.0};
        Polynomial_init(&p4, 2, multipliers, 5);
        
        Polynomial pb;
        Polynomial_init(&pb, 0, b, nb);
        Polynomial pa;
        Polynomial_init(&pa, 0, a, na);

        b++;
        p4.order = 2;
        //printf("p2 order %zu, pb order %zu, p2 ncoefs %zu, pb ncoefs %zu\n", p2.order, pb. order, p2.ncoefs, pb.ncoefs);
        for (size_t i = 0; i < M1; i++) {
            p4.coefs[1] = *b;
            b++;
            p4.coefs[2] = *b;
            b++;
            Polynomial_mul(&pb, &p4);
        }
        p4.order = 4;
        for (size_t i = M1; i < M2; i++) {
            p4.coefs[1] = *b;
            b++;
            p4.coefs[2] = *b;
            b++;
            p4.coefs[3] = *b;
            b++;
            p4.coefs[4] = *b;
            b++;
            Polynomial_mul(&pb, &p4);
        }
        p4.coefs[3] = 0.0;
        p4.coefs[4] = 0.0;

        a++;
        p4.order = 2;
        //printf("p2 order %zu, pb order %zu, p2 ncoefs %zu, pb ncoefs %zu\n", p2.order, pb. order, p2.ncoefs, pb.ncoefs);
        for (size_t i = 0; i < N1; i++) {
            p4.coefs[1] = *a;
            a++;
            p4.coefs[2] = *a;
            a++;
            Polynomial_mul(&pa, &p4);
        }
        p4.order = 4;
        for (size_t i = N1; i < N2; i++) {
            p4.coefs[1] = *a;
            a++;
            p4.coefs[2] = *a;
            a++;
            p4.coefs[3] = *a;
            a++;
            p4.coefs[4] = *a;
            a++;
            Polynomial_mul(&pa, &p4);
        }
        p4.coefs[3] = 0.0;
        p4.coefs[4] = 0.0;
    } else { // error...what the hell do you expect me to do
        return 1;
    }
    //IIRFilterBank_print((FilterBank *) ifb);
    //b = IIRFilterBank_get_b(ifb);
    //double sumb = 0.0;
    //a = IIRFilterBank_get_a(ifb);
    //double suma = 0.0;
    //for (size_t i = 0; i < IIRFilterBank_get_b_size(ifb); i++) {
    //    sumb += b[i];
    //    suma += a[i];
    //}
    //printf("sum b = %.8f\nsum a = %.8f\n", sumb, suma);
    return 0;
}

void butterworth_digital_prototype(IIRFilterBank * ifb, size_t L, double w0) {
    if (!ifb) {
        return;
    }
    size_t N1 = (L & 1) ? 1 : 0;
    //size_t M1 = 0;
    size_t N2 = N1 + (L - N1) / 2;
    ifb->fb.nb = L + 1;
    ifb->na = L + 1;
    //size_t M2 = 0;
    //printf("N1 %zu, N %zu, M1 %zu, M %zu\n", N1, N, M1, M);
    double * pgain = IIRFilterBank_get_b(ifb); // skip first point, which will be for the z^0 component in the filterbank
    double * zerocoefs = pgain + 1;
    double * polecoefs = IIRFilterBank_get_a(ifb); // skip first point, which will be for the z^0 component in the filterbank
    *polecoefs = 1.0;
    polecoefs++;
    // encode the number of real zero and pole coefficients in the filter bank. This violates the class definition of filter bank...but it's not a class
    //IIRFilterBank_print(fb);
    ifb->fb.nb = L;
    //printf("(polecoefs - zerocoefs) / 8 = %d\n", (polecoefs - zerocoefs) / 8);
    ifb->na = N1;
    size_t i = 0;
    double W = tan(M_PI * w0 / 2);
    double W2 = pow(W, 2.0);
    *pgain = 1.0; // gain is always 0 for Butterworth
    //printf("starting loop to calculate coefficients: L %zu, M1 %zu, M %zu, N1 %zu, N %zu\n", *L, M1, M, N1, N);
    while (i < N2) { // changed from L -->N2 because we do not have to worry so much about after N2 for butterworth
        double num = 1.0, den = 1.0;
        //printf("calculating poles coefficient at i=%zu\n", i);
        if (i < N1) {
            num = W;
            //printf("processing an imaginary pole\n");
            den = 1 + W;
            //if (!polecoefs) {
            //    printf("attempting to dereference null pointer in polecoefs at i=%zu < M1 %zu\n", i, N1);
            //}
            *polecoefs = -(1 - W) / den;
            polecoefs++;
            *zerocoefs = 1.0; //nan("");
            zerocoefs++;
        } else if (i < N2) {
            num = W2;
            //printf("processing complex pole and conjugate\n");
            double zi = 2.0 * W * sin( M_PI * (2 * (i - N1) + 1) / (2 * L));
            double z2 = W2;
            den = 1 + zi + z2;
            //if (!polecoefs || !(polecoefs+1)) {
            //    printf("attempting to dereference null pointer in polecoefs at i=%zu < M1 %zu\n", i, Ncomp);
            //}
            *polecoefs = -2.0 * (1 - z2) / den;
            polecoefs++;
            *polecoefs = (1 - zi + z2) / den;
            polecoefs++;
            *zerocoefs = 1.0;
            zerocoefs++;
            *zerocoefs = 1.0;
            zerocoefs++;
        }
        (*pgain) *= num / den;
        i++;
    }
    //printf("computed zero and pole coefficients successfully...returning\n");
    //return 0;
}

/* static allocation
wl and wu define the cutoff frequences for 1/2 power attenuation.
wl = 0 --> low-pass filter
wu = 0 --> high-pass filter
wl < wu --> band-pass filter
wl > wu --> band-stop filter
wl == wu
wl and wu are in units of the Nyquist frequency 1/(2*T)
*/
int butterworth(RTIIRFilter * rtif, size_t order, double wl, double wu, unsigned int flags, int (*initialize)(RTFilter * rtf, double sample)) {
    if (!rtif || wl == wu || !order) {
        return 1;
    }
    size_t L = ((wu > 0.0 && wl > 0.0) ? 2 : 1) * order + 1;
    if (IIRFilterBank_get_a_size(&rtif->ifb) < L || IIRFilterBank_get_b_size(&rtif->ifb) < L)  {
        return 2;
    }
    if (!initialize) {
        rtif->rtf.initialize = RTIIRFilter_stable_init;
    }

    //printf("order %zu, mult * order + 1, %zu\n", order, mult * order + 1);
    double w0 = wu;
    if (w0 < wl) {
        w0 = wl;
    }
    //printf("w0 %.4f, wl %.4f, wu %.4f\n", w0, wl, wu);
    butterworth_digital_prototype(&rtif->ifb, order, w0);
    digital_prototype_to_IIRFilterBank(&rtif->ifb, order, w0, wl, wu);

    return 0;
}

int chebyshev1_digital_prototype(IIRFilterBank * ifb, size_t L, double ripple, double w0) {
    
    if (!ifb) {
        return 1;
    }
    size_t N1 = (L & 1) ? 1 : 0;
    //size_t M1 = 0;
    size_t N2 = N1 + (L - N1) / 2;
    ifb->fb.nb = L + 1;
    ifb->na = L + 1;
    //size_t M2 = 0;
    //printf("N1 %zu, N %zu, M1 %zu, M %zu\n", N1, N, M1, M);
    double * pgain = IIRFilterBank_get_b(ifb); // skip first point, which will be for the z^0 component in the filterbank
    double * zerocoefs = pgain + 1;
    double * polecoefs = IIRFilterBank_get_a(ifb); // skip first point, which will be for the z^0 component in the filterbank
    *polecoefs = 1.0;
    polecoefs++;
    // encode the number of real zero and pole coefficients in the filter bank. This violates the class definition of filter bank...but it's not a class
    //IIRFilterBank_print(fb);
    ifb->fb.nb = L;
    //printf("(polecoefs - zerocoefs) / 8 = %d\n", (polecoefs - zerocoefs) / 8);
    ifb->na = N1;
    size_t i = 0;
    double coshf = cosh(asinh(1.0 / ripple) / L);
    double sinhf = sinh(asinh(1.0 / ripple) / L);
    double W = tan(M_PI * w0 / 2);
    double wsinhf = W * sinhf;
    double W2 = pow(W, 2);
    //*pgain = 1.0; 
    *pgain = (L & 1) ? 1.0 : 1.0/pow(ripple * ripple + 1, .5); // gain for chebyshev1 must be modified since the DC is not necessarily 0
    //printf("starting loop to calculate coefficients: L %zu, M1 %zu, M %zu, N1 %zu, N %zu\n", *L, M1, M, N1, N);
    while (i < N2) { // changed from L -->N2 because we do not have to worry so much about after N2 for butterworth
        double num = 1.0, den = 1.0;
        //printf("calculating poles coefficient at i=%zu\n", i);
        if (i < N1) {
            num = wsinhf;
            //printf("processing an imaginary pole\n");
            den = 1 + wsinhf;
            //if (!polecoefs) {
            //    printf("attempting to dereference null pointer in polecoefs at i=%zu < M1 %zu\n", i, N1);
            //}
            *polecoefs = -(1 - wsinhf) / den;
            polecoefs++;
            *zerocoefs = 1.0; //nan("");
            zerocoefs++;
        } else if (i < N2) { // something is wrong here. The normalization looks OK, but the poles are wrong
            double phase = M_PI * (2 * (i - N1 * 1.0) + 1) / (2 * L);
            double sinphase = sin(phase);
            double cosphase = cos(phase);
            //printf("processing complex pole and conjugate\n");
            double rhok2 =  pow(cosphase * coshf, 2) + pow(sinphase * sinhf, 2);
            double rhokim = 2 * wsinhf * sinphase;
            num = W2 * rhok2;
            den = 1 + rhokim + num;
            //if (!polecoefs || !(polecoefs+1)) {
            //    printf("attempting to dereference null pointer in polecoefs at i=%zu < M1 %zu\n", i, Ncomp);
            //}
            *polecoefs = -2.0 * (1 - num) / den;
            polecoefs++;
            *polecoefs = (1 - rhokim + num) / den;
            polecoefs++;
            *zerocoefs = 1.0;
            zerocoefs++;
            *zerocoefs = 1.0;
            zerocoefs++;
        }
        (*pgain) *= num / den;
        i++;
    }
    //printf("computed zero and pole coefficients successfully...returning\n");
    return 0;
}

int chebyshev1(RTIIRFilter * rtif, size_t order, double ripple, double wl, double wu, unsigned int flags, int (*initialize)(RTFilter * rtf, double sample)) {
    if (!rtif || wl == wu || !order) {
        return 1;
    }
    size_t L = ((wu > 0.0 && wl > 0.0) ? 2 : 1) * order + 1;
    if (IIRFilterBank_get_a_size(&rtif->ifb) < L || IIRFilterBank_get_b_size(&rtif->ifb) < L)  {
        return 2;
    }
    if (!initialize) {
        initialize = RTIIRFilter_stable_init;
    }

    //printf("order %zu, mult * order + 1, %zu\n", order, mult * order + 1);
    double w0 = wu;
    if (w0 < wl) {
        w0 = wl;
    }
    //printf("w0 %.4f, wl %.4f, wu %.4f\n", w0, wl, wu);
    chebyshev1_digital_prototype(&rtif->ifb, order, ripple, w0);
    digital_prototype_to_IIRFilterBank(&rtif->ifb, order, w0, wl, wu);

    return 0;
}

//TODO: this implements chebyshev1...need to test chebyshev1 and implement chebyshev2 here
int chebyshev2_digital_prototype(IIRFilterBank * ifb, size_t L, double ripple, double w0) {
    if (!ifb) {
        return 1;
    }
    size_t N1 = (L & 1) ? 1 : 0;
    size_t M1 = N1;
    size_t N2 = N1 + (L - N1) / 2;
    ifb->fb.nb = L + 1;
    ifb->na = L + 1;
    size_t M2 = M1 + (L - M1) / 2;
    //printf("N1 %zu, N %zu, M1 %zu, M %zu\n", N1, N, M1, M);
    double * pgain = IIRFilterBank_get_b(ifb); // skip first point, which will be for the z^0 component in the filterbank
    double * zerocoefs = pgain + 1;
    double * polecoefs = IIRFilterBank_get_a(ifb); // skip first point, which will be for the z^0 component in the filterbank
    *polecoefs = 1.0;
    polecoefs++;
    // encode the number of real zero and pole coefficients in the filter bank. This violates the class definition of filter bank...but it's not a class
    //IIRFilterBank_print(fb);
    ifb->fb.nb = M1;
    //printf("(polecoefs - zerocoefs) / 8 = %d\n", (polecoefs - zerocoefs) / 8);
    ifb->na = N1;
    size_t i = 0;
    double coshf2 = pow(cosh(asinh(1.0 / ripple) / L), 2); // only the square appears
    double sinhf = sinh(asinh(1.0 / ripple) / L);
    double W = tan(M_PI * w0 / 2);
    double wsinhf = W * sinhf;
    double W2 = pow(W, 2);
    //printf("W = %.8f, sinhf = %.8f, coshf^2 = %.8f\n", W, sinhf, coshf2);
    *pgain = 1.0; // gain is always 0 for Butterworth
    //printf("starting loop to calculate coefficients: L %zu, M1 %zu, M %zu, N1 %zu, N %zu\n", *L, M1, M, N1, N);
    while (i < ((N2 >= M2) ? N2 : M2)) { // in Chebyshev2 case, all zeros are real so M2 > N2
        double num = 1.0, den = 1.0;
        //printf("calculating poles & zeros coefficient at i=%zu\n", i);
        if (i < M1) {
            num = W;
            *zerocoefs = 1.0;
            zerocoefs++;
        } else if (i < M2) {
            double cos2 = pow(cos(M_PI * (2 * (i - M1 * 1.0) + 1) / (2 * L)), 2);
            num = cos2 + W2;
            *zerocoefs = -2 * (cos2 - W2) / num;
            zerocoefs++;
            *zerocoefs = 1.0;
            zerocoefs++;
        }
        if (i < N1) {
            printf("processing an imaginary pole\n");
            double phase = M_PI / 2;
            double sinphase = sin(phase);
            double cosphase = cos(phase);
            double rhok2 =  pow(cosphase, 2) * coshf2 + pow(sinphase * sinhf, 2);
            printf("processing complex pole and conjugate: %.8f + i %.8f\n", cosphase * pow(coshf2, .5)/rhok2, sinphase * sinhf/rhok2);
            den = W + sinhf;
            *polecoefs = -(sinhf - W) / den;
            polecoefs++;
        } else if (i < N2) {
            double phase = M_PI * (2 * (i - N1 * 1.0) + 1) / (2 * L);
            double sinphase = sin(phase);
            double cosphase = cos(phase);
            
            double rhok2 =  pow(cosphase, 2) * coshf2 + pow(sinphase * sinhf, 2);
            double rhokim = 2 * wsinhf * sinphase;
            printf("processing complex pole and conjugate: %.8f + i %.8f\n", cosphase * pow(coshf2, .5)/rhok2, sinphase * sinhf/rhok2);
            den = rhok2 + rhokim + W2;
            *polecoefs = -2.0 * (rhok2 - W2) / den;
            polecoefs++;
            *polecoefs = (rhok2 - rhokim + W2) / den;
            polecoefs++;
            printf("phase = %.8f, rhokim = %.8f, rhok2 = %.8f, A1 = %.8f, A2 = %.8f\n", phase, rhokim, rhok2, *(polecoefs-1), *polecoefs);
        }
        (*pgain) *= num / den;
        i++;
    }
    //printf("computed zero and pole coefficients successfully...returning\n");
    return 0;
}

int chebyshev2(RTIIRFilter * rtif, size_t order, double ripple, double wl, double wu, unsigned int flags, int (*initialize)(RTFilter * rtf, double sample)) {
    if (!rtif || wl == wu || !order) {
        return 1;
    }
    size_t L = ((wu > 0.0 && wl > 0.0) ? 2 : 1) * order + 1;
    if (IIRFilterBank_get_a_size(&rtif->ifb) < L || IIRFilterBank_get_b_size(&rtif->ifb) < L)  {
        return 2;
    }
    if (!initialize) {
        initialize = RTIIRFilter_stable_init;
    }

    //printf("order %zu, mult * order + 1, %zu\n", order, mult * order + 1);
    double w0 = wu;
    if (w0 < wl) {
        w0 = wl;
    }
    //printf("w0 %.4f, wl %.4f, wu %.4f\n", w0, wl, wu);
    chebyshev2_digital_prototype(&rtif->ifb, order, ripple, w0);
    digital_prototype_to_IIRFilterBank(&rtif->ifb, order, w0, wl, wu);

    return 0;
}

int PID(RTIIRFilter * rtif, double kp, double ki, double kd) {
    if (!rtif || RTIIRFilter_get_a_size(rtif) < 2 || RTIIRFilter_get_b_size(rtif) < 3) {
        return 1;
    }
    double * b = RTIIRFilter_get_b(rtif);
    double * a = RTIIRFilter_get_a(rtif);
    a[0] = 1.0;
    a[1] = -1.0;
    b[0] = kp + ki + kd;
    b[1] = -1.0 * (kp + 2 * kd);
    b[2] = kd;
    return 0;
}

// TODO: 
int thiran_digital_prototype(IIRFilterBank * ifb, size_t order, double tau) {
    return 0;
}

// this only makes a low-pass filter
int thiran(RTIIRFilter * rtif, size_t order, double tau, int (*initialize)(RTFilter * rtf, double sample)) {
    if (!rtif || tau == 0 || !order) {
        return 1;
    }
    if (IIRFilterBank_get_a_size(&rtif->ifb) < 1 || IIRFilterBank_get_b_size(&rtif->ifb) < order + 1)  {
        return 2;
    }
    if (!initialize) {
        rtif->rtf.initialize = RTIIRFilter_stable_init;
    }

    //printf("order %zu, mult * order + 1, %zu\n", order, mult * order + 1);
    // need to figure out how to identify the "cutoff frequency and be able to convert it to high-pass & band-pass"
    //printf("w0 %.4f, wl %.4f, wu %.4f\n", w0, wl, wu);

    rtif->ifb.fb.nb = 1;
    double * b = RTIIRFilter_get_b(rtif);
    b[0] = 1.0;
    for (size_t i = order + 1; i < 2 * order + 1; i++) {
        b[0] *= i / (2 * tau + i); // tau is double so not integer division
    }
    rtif->ifb.na = order + 1;
    double * a = RTIIRFilter_get_a(rtif);
    size_t nck = 1;
    for (size_t k = 0; k < order + 1; k++) {
        a[k] = (k & 1 ? -1.0 : 1.0);
        if (k) {
            nck = nck * (order - k + 1) / (k);
        }
        for (size_t i = 0; i < order + 1; i++) {
            a[k] *= (2 * tau + i) / (2 * tau + k + i);
        }
        a[k] *= nck;
    }

    return 0;
}

#ifndef __STDC_NO_COMPLEX__

#define DEFAULT_COMPLEX_TOLERANCE 1e-7

// TODO: not the appropriate module for this function
// TODO: this might have to not just set to 0 but rotate the complex number to align, treating specially the case of 0.0 + i0.0
void flush_complex_to_zero(double complex * arr, size_t N, double tolerance) {
    if (!arr || !N) {
        return;
    }
    if (tolerance <= 0.0) {
        tolerance = DEFAULT_COMPLEX_TOLERANCE;
    }
    //printf("poles after flush to zero\n");
    for (size_t k = 0; k < N; k++) {
        double r = creal(arr[k]);
        if (fabs(r) < tolerance) {
            r = 0.0;
        }
        double i = cimag(arr[k]);
        if (fabs(i) < tolerance) {
            i = 0.0;
        }
        arr[k] = r + I * i;
        //printf("%.4f + i %.4f\n", creal(arr[k]), cimag(arr[k]));
    }
}

double pzg_gain(double complex * poles, size_t N, double complex * zeros, size_t M) {
    //return 1.0;
    size_t nmmin = 0;
    size_t nmmax = 0;
    if (N >= M) {
        nmmin = M;
        nmmax = N;
    } else {
        nmmin = N;
        nmmax = M;
    }
    double complex gain = 1.0;
    for (size_t i = 0; i < nmmin; i++) {
        gain *= zeros[i] / poles[i];
    }
    for (size_t i = nmmin; i < nmmax; i++) {
        gain *= (i < M ? -zeros[i] : 1.0) / (i < N ? -poles[i] : 1.0);
    }
    return 1.0 / sqrt(fabs(creal(gain)));
}

// private
// returns 0 if the array represented by zeros has a real product meaning all elements are either real or 
int pzg_to_RTIIRFilter_check_complete(double complex * arr, size_t N, double tolerance) {
    if (!arr) {
        return 0;
    }
    if (tolerance <= 0.0) {
        tolerance = DEFAULT_COMPLEX_TOLERANCE;
    }
    double complex prod = 1.0;
    for (size_t i = 0; i < N; i++) {
        prod *= arr[i];
    }
    //printf("product of zeros: %.4f + i%.4f\n", creal(prod), cimag(prod));
    return cimag(prod) > tolerance;
}

// private
// returns 0 if succes, in which case arr is in an unknown arrangement, else returns the number of elements in the first quadrant, above real axis.
// results are also sorted such that arr[:*Nimag] have no real component and arr[*Nimag:*Nimag + *Ncomplex] have a real component
int pzg_to_RTIIRFilter_sort_q1(size_t * Nimag, size_t * Ncomplex, double complex * arr, size_t N, double tolerance) {
    *Nimag = 0;
    *Ncomplex = 0;
    if (!arr || !N) {
        return 1;
    } 
    if (tolerance <= 0.0) {
        tolerance = DEFAULT_COMPLEX_TOLERANCE;
    }
    size_t iq1i = 0;
    size_t iq1c = N;
    size_t i = 0;
    while (i < iq1c) {
        if (cimag(arr[i]) >= 0.0) { // TODO: see chebyshev2 for a reason to include these zeros
            double complex temp = arr[i];
            double r = creal(temp);
            if (fabs(r) <= tolerance) {
                arr[i] = arr[iq1i];
                arr[iq1i] = temp;
                iq1i++;
                (*Nimag)++;
            } else if (r > tolerance) {
                iq1c--;
                //printf("complex found: %.4f + i*%.4f\n", r, cimag(arr[i]));
                arr[i] = arr[iq1c];
                arr[iq1c] = temp;
                (*Ncomplex)++;
            }
        }
        i++;
    }
    //printf("Ncomplex %zu\n", *Ncomplex);
    for (size_t i = 0; i < *Ncomplex; i++) {
        double complex temp = arr[(*Nimag) + i];
        arr[(*Nimag) + i] = arr[N - 1 - i];
        arr[N - 1 - i] = temp;
    }
    return 0;
}

int pz_transform_preprocess(double complex * arr, size_t * N) {
    size_t Nimag = 0, Ncomplex = 0;
    if (*N) {
        flush_complex_to_zero(arr, *N, 0.0);
        if (pzg_to_RTIIRFilter_check_complete(arr, *N, 0.0)) {
            //printf("failed RTIIRFilter_check_complete on the poles");
            return 1;
        }
        if (pzg_to_RTIIRFilter_sort_q1(&Nimag, &Ncomplex, arr, *N, 0.0)) {
            //printf("failed RTIIRFilter_sort_q1 on the poles");
            return 1;
        }
    }
    *N = Nimag + 2 * Ncomplex;
    return 0;
}

size_t preprocessed_pz_nimag(double complex * arr, size_t N) {
    size_t Nreal = 0;
    while (Nreal < N && creal(arr[Nreal]) == 0.0) {
        Nreal++;
    }
    return Nreal;
}

// TODO: need an error code to tell user that the filter bank was not allocated to a large enough size
// pass only the poles and zeros that do not have complex conjugates or just one of the conjugates. Those that do not, must be first. 
// preprocess to achieve this with pz_transform_preprocess on both arrays
// the 1+z^-1 terms have coefficients stored as nans in the digital prototype
int pzg_to_digital_prototype(IIRFilterBank * ifb, double w0, double complex * poles, size_t N, double complex * zeros, size_t M, double gain) {
    //printf("in pzg_to_digital_prototype\n");
    if (!ifb) {
        return 1;
    }
    if (gain <= 0.0) {
        gain = 1.0;
    }
    size_t N1 = preprocessed_pz_nimag(poles, N);
    size_t M1 = preprocessed_pz_nimag(zeros, M);
    //printf("N1 %zu, N %zu, M1 %zu, M %zu\n", N1, N, M1, M);
    double * zerocoefs = NULL, * polecoefs = NULL, * pgain = NULL;
    size_t L = 0;
    if (N >= M) {
        //*L = N;
        L = N;
        // this next line probably should not be done;
        IIRFilterBank_init(ifb, IIRFilterBank_get_b(ifb), N + 1, N + 1);
        // add (N-M) zeros with value -1
        pgain = IIRFilterBank_get_b(ifb); // skip first point, which will be for the z^0 component in the filterbank
        zerocoefs = pgain + 1;
        polecoefs = IIRFilterBank_get_a(ifb); // skip first point, which will be for the z^0 component in the filterbank
        *polecoefs = 1.0;
        polecoefs++;
        for (size_t i = 0; i < N-M; i++) {
            *zerocoefs = 1.0; //nan("");
            zerocoefs++;
        }
        // encode the number of real zero and pole coefficients in the filter bank. This violates the class definition of filter bank...but it's not a class
        //IIRFilterBank_print(fb);
        ifb->fb.nb = M1 + (N - M);
        //printf("(polecoefs - zerocoefs) / 8 = %d\n", (polecoefs - zerocoefs) / 8);
        ifb->na = N1;
    } else {
        //*L = M;
        L = M;
        // this next line probably should not be done;
        IIRFilterBank_init(ifb, IIRFilterBank_get_b(ifb), M + 1, M + 1);
        // add (N-M) zeros with value -1
        pgain = IIRFilterBank_get_b(ifb); // skip first point, which will be for the z^0 component in the filterbank
        zerocoefs = pgain + 1;
        polecoefs = IIRFilterBank_get_a(ifb); // skip first point, which will be for the z^0 component in the filterbank
        *polecoefs = 1.0;
        polecoefs++;
        for (size_t i = 0; i < N-M; i++) {
            *polecoefs = 1.0; //nan("");
            polecoefs++;
        }
        // encode the number of real zero and pole coefficients in the filter bank. This violates the class definition of filter bank...but it's not a class
        ifb->fb.nb = M1;
        ifb->na = N1 + (N - M);
    }
    //printf("initialized zero and pole coefficient arrays\n");
    //N = N1 + Ncomp;
    size_t N2 = N - (N - N1) / 2;
    //M = M1 + Mcomp;
    size_t M2 = M - (M - M1) / 2;
    size_t i = 0;
    double W = tan(M_PI * w0 / 2);
    //if (!pgain) {
    //    printf("pgain is NULL pointer\n");
    //}
    *pgain = gain;
    //printf("starting loop to calculate coefficients: L %zu, M1 %zu, M %zu, N1 %zu, N %zu\n", *L, M1, M, N1, N);
    while (i < L) {
        double num = 1.0, den = 1.0;
        //printf("calculating zeros coefficient at i=%zu\n", i);
        if (i < M1) {
            //printf("assigning something to zerocoefs[%zu]\n", (zerocoefs - pgain) / 8);
            double zi = W * cimag(zeros[i]);
            num = 1 + zi;
            //if (!zerocoefs) {
            //    printf("attempting to dereference null pointer in zerocoefs at i=%zu < M1 %zu\n", i, M1);
            //}
            *zerocoefs = -(1 - zi) / num;
            zerocoefs++;
        } else if (i < M2) {
            double zi = 2.0 * W * cimag(zeros[i]);
            double z2 = pow(W * cabs(zeros[i]), 2);
            num = 1 + zi + z2;
            //if (!zerocoefs || !(zerocoefs+1)) {
            //    printf("attempting to dereference null pointer in zerocoefs at i=%zu < M1 %zu\n", i, Mcomp);
            //}
            *zerocoefs = -2.0 * (1 - z2) / num;
            zerocoefs++;
            *zerocoefs = (1 - zi + z2) / num;
            zerocoefs++;
        } else if (N >= M) {
            //printf("multiply by W, %.4f\n", W);
            num = W;
        }
        //printf("calculating poles coefficient at i=%zu\n", i);
        if (i < N1) {
            //printf("processing an imaginary pole\n");
            double zi = W * cimag(poles[i]);
            den = 1 + zi;
            //if (!polecoefs) {
            //    printf("attempting to dereference null pointer in polecoefs at i=%zu < M1 %zu\n", i, N1);
            //}
            *polecoefs = -(1 - zi) / den;
            polecoefs++;
        } else if (i < N2) {
            //printf("processing complex pole and conjugate\n");
            double zi = 2.0 * W * cimag(poles[i]);
            double z2 = pow(W * cabs(poles[i]), 2);
            den = 1 + zi + z2;
            //if (!polecoefs || !(polecoefs+1)) {
            //    printf("attempting to dereference null pointer in polecoefs at i=%zu < M1 %zu\n", i, Ncomp);
            //}
            *polecoefs = -2.0 * (1 - z2) / den;
            polecoefs++;
            *polecoefs = (1 - zi + z2) / den;
            polecoefs++;
        } else if (M > N) {
            //printf("dividing by W, %.4f\n", W);
            den = W;
        }
        (*pgain) *= num / den;
        i++;
    }
    //printf("computed zero and pole coefficients successfully...returning\n");
    return 0;
}

// internal only
int pzg_to_RTIIRFilter(RTIIRFilter * rtif, double complex * poles, size_t N, double complex * zeros, size_t M, double gain, double wl, double wu) {
    if (!rtif || !(poles || zeros) || !(M || N) || wl == wu) {
        return 1;
    }
    
    gain *= pzg_gain(poles, N, zeros, M);
    printf("gain %.4f\n", gain);
    printf("before: N = %zu, M = %zu\n", N, M);
    
    pz_transform_preprocess(poles, &N);
    printf("selected poles\n");
    for (size_t i = 0; i < N; i++) {
        printf("%.4f +i%.4f, ", creal(poles[i]), cimag(poles[i]));
    }
    printf("\nselected zeros\n");
    pz_transform_preprocess(zeros, &M);
    for (size_t i = 0; i < M; i++) {
        printf("%.4f +i%.4f, ", creal(zeros[i]), cimag(zeros[i]));
    }
    printf("\n");

    printf("after: N = %zu, M = %zu\n", N, M);

    size_t L = (N >= M) ? N : M;
    size_t na = ((wu > 0.0 && wl > 0.0) ? 2 : 1) * L + 1;
    size_t nb = na;
    if (RTIIRFilter_get_a_size(rtif) < na || RTIIRFilter_get_b_size(rtif) < nb) {
        return 2;
    }
    //RTIIRFilter_init(rtif, RTIIRFilter_get_b(rtif), rtif->state, na, nb, 0, NULL);
    if (gain <= 0.0) {
        gain = 1.0;
    }

    
    //printf("L %zu, 2 * L + 1 %zu\n", L, 2 * L + 1);
    /*
    if (wu != 0.0 && wl != 0.0) { // band pass or band-stop
        rtif = RTIIRFilter_new_empty(2 * L + 1, 2 * L + 1, 0, RTIIRFilter_stable_init);
    } else {
        rtif = RTIIRFilter_new_empty(L + 1, L + 1, 0, RTIIRFilter_stable_init);
    }
    if (!rtif) {
        return NULL;
    }
    */
    IIRFilterBank ifb = rtif->ifb;
    double w0 = wu;
    if (w0 < wl) {
        w0 = wl;
    }
    //printf("w0 %.4f, wl %.4f, wu %.4f\n", w0, wl, wu);
    pzg_to_digital_prototype(&ifb, w0, poles, N, zeros, M, gain);
    //pzg_to_digital_prototype(IIRFilterBank * ifb, double w0, double complex * poles, size_t N, double complex * zeros, size_t M, double gain)
    
    
    printf("pzg bank\n");
    double * bank = ifb.fb.b;
    for (size_t i = 0; i < 2 * (L + 1); i++) {
        printf("%.8f ", bank[i]);
    }
    printf("\n");
    
    digital_prototype_to_IIRFilterBank(&ifb, L, w0, wl, wu);
    //int digital_prototype_to_IIRFilterBank(IIRFilterBank * ifb, size_t M1, size_t N1, double wl, double wu) {
    return 1;
}

#endif // __STDC_NO_COMPLEX
