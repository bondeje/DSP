#include <RTFilter.h>

#ifndef __STDC_NO_COMPLEX__
EXPORT int filter_response_ab(double complex * gain, size_t ng, double * a, size_t na, double * b, size_t nb, double * freq);
EXPORT int filter_response_pzg(double complex * gain, size_t ng, double complex * zeros, size_t nz, double complex * poles, size_t np, double kgain, double * freq);
#endif
EXPORT int filter_response_ab_noc(double * gain, double * phase, size_t ng, double * a, size_t na, double * b, size_t nb, double * freq);
EXPORT int filter_response_pzg_noc(double * gain, double * phase, size_t ng, double complex * zeros, size_t nz, double complex * poles, size_t np, double kgain, double * freq);