#include <filterutils.h>
#include <polynomials.h>
#include <stdlib.h>

#ifndef __STDC_NO_COMPLEX__

int filter_response_ab(double complex * gain, size_t ng, double * a, size_t na, double * b, size_t nb, double * freq) { 
    if (!gain) {
        return 1;
    }
    if (!ng) {
        return 0;
    }
    double * freq_ = NULL;
    if (!freq) {
        freq_ = (double *) DSP_MALLOC(ng * sizeof(double));
        if (!freq_) {
            return 1;
        }
        double df = 1.0 / (ng-1);
        freq_[0] = 0.0;
        for (size_t i = 1; i < ng; i++) {
            freq_[i] = freq_[i-1] + df;
        }
    }

    /* convert transfer function coefficients to polynomials for evaluation */
    Polynomial pa;
    Polynomial_init(&pa, na-1, a, na);
    Polynomial pb;
    Polynomial_init(&pb, nb-1, b, nb);

    /* evaluate transfer function */
    for (size_t i = 0; i < ng; i++) {
        double complex ztom1 = cexp(-I * freq_[i]);
        gain[i] = Polynomial_ceval(&pb, ztom1) / Polynomial_ceval(&pb, ztom1);
    }
    if (!freq) { // freq_ was malloc'd
        DSP_FREE(freq_);
    }
    return 0;
}

int filter_response_pzg(double complex * gain, size_t ng, double complex * zeros, size_t nz, double complex * poles, size_t np, double kgain, double * freq) {
    if (!gain) {
        return 1;
    }
    if (!ng) {
        return 0;
    }
    double * freq_ = NULL;
    if (!freq) {
        freq_ = (double *) DSP_MALLOC(ng * sizeof(double));
        if (!freq_) {
            return 1;
        }
        double df = 1.0 / (ng-1);
        freq_[0] = 0.0;
        for (size_t i = 1; i < ng; i++) {
            freq_[i] = freq_[i-1] + df;
        }
    } else {
        freq_ = freq;
    }

    size_t n = (nz > np ? nz : np);
    /* evaluate transfer function */
    for (size_t i = 0; i < ng; i++) {
        double complex z = cexp(I * freq_[i]);
        gain[i] = kgain;    
        for (size_t j = 0; j < n; j++) { 
            double complex factor = 1.0;
            if (j < nz) {
                factor *= (z - zeros[j]);
            }
            if (j < np) {
                factor /= (z - poles[j]);
            }
            gain[i] *= factor;
        }
    }
    if (!freq) { // freq_ was malloc'd
        DSP_FREE(freq_);
    }
    return 0;
}

#endif // __STDC_NO_COMPLEX__

/* phase in [-pi, pi] */
/* gain and phase are double arrays with at least as much  memory as sizeof(double) * ng */
/* freq are optional but must be at least as big as ng */
int filter_response_ab_noc(double * gain, double * phase, size_t ng, double * a, size_t na, double * b, size_t nb, double * freq) { 
    if (!gain || !phase) {
        return 1;
    }
    if (!ng) {
        return 0;
    }
    int status = 0;
#ifndef __STDC_NO_COMPLEX__
    /* wrapper for filter_gain */
    double complex * gainc = NULL;
    if (!freq) {
        gainc = (double complex *) DSP_MALLOC(ng * (sizeof(double complex) + sizeof(double)));
        if (!gainc) {
            return 1;
        }
        freq = (double *) (gainc + ng); // freq is at the end of the gainc. this might violate strict aliasing
        double df = 1.0 / (ng-1);
        freq[0] = 0.0;
        for (size_t i = 1; i < ng; i++) {
            freq[i] = freq[i-1] + df;
        }
    } else {
        gainc = (double complex *) DSP_MALLOC(ng * sizeof(double complex));
        if (!gainc) {
            return 1;
        }
    }
    
    if (!(status = filter_response_ab(gainc, ng, a, na, b, nb, freq))) {
        /* unwrap gainc */
        for (size_t i = 0; i < ng; i++) {
            gain[i] = cabs(gainc[i]);
            phase[i] = carg(gainc[i]);
        }
    }
    DSP_FREE(gainc); /* never have to free freq since if it was malloc'd it's on top of gainc pointer*/
#else
/* implement complex-free version */
#error complex.h not available and filter_response_ab_noc not yet implemented
#endif
    return status;
}

int filter_response_pzg_noc(double * gain, double * phase, size_t ng, double complex * zeros, size_t nz, double complex * poles, size_t np, double kgain, double * freq) { 
    if (!gain || !phase) {
        return 1;
    }
    if (!ng) {
        return 0;
    }
    int status = 0;
#ifndef __STDC_NO_COMPLEX__
    /* wrapper for filter_gain */
    double complex * gainc = NULL;
    if (!freq) {
        gainc = (double complex *) DSP_MALLOC(ng * (sizeof(double complex) + sizeof(double)));
        if (!gainc) {
            return 1;
        }
        freq = (double *) (gainc + ng); // freq is at the end of the gainc. this might violate strict aliasing
        double df = 1.0 / (ng-1);
        freq[0] = 0.0;
        for (size_t i = 1; i < ng; i++) {
            freq[i] = freq[i-1] + df;
        }
    } else {
        gainc = (double complex *) DSP_MALLOC(ng * sizeof(double complex));
        if (!gainc) {
            return 1;
        }
    }
    
    if (!(status = filter_response_pzg(gainc, ng, zeros, nz, poles, np, kgain, freq))) {
        /* unwrap gainc */
        for (size_t i = 0; i < ng; i++) {
            gain[i] = cabs(gainc[i]);
            phase[i] = carg(gainc[i]);
        }
    }
    DSP_FREE(gainc); /* never have to free freq since if it was malloc'd it's on top of gainc pointer*/
#else
/* implement complex-free version */
#error complex.h not available and filter_response_pzg_noc not yet implemented
#endif
    return status;
}