#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include <stddef.h>
#include <stdarg.h>
#include <dsputils.h>

#ifdef DEBUG_MALLOC
    #include "../../mallocs/debug_malloc.h"
#endif

#ifndef POLY_FREE
    #define POLY_FREE free
#endif
#ifndef POLY_MALLOC
    #define POLY_MALLOC malloc
#endif
#ifndef POLY_CALLOC
    #define POLY_CALLOC calloc
#endif
#ifndef POLY_REALLOC
    #define POLY_REALLOC realloc
#endif


// only implements polynomials with real coefficients

typedef struct Polynomial Polynomial;

struct Polynomial {
    double * coefs;
    size_t ncoefs; // size of allocation of coefs
    size_t order; // coefs must be an array of len order + 1
};

EXPORT void Polynomial_print(Polynomial * poly);

EXPORT void Polynomial_collapse_order(Polynomial * poly);

EXPORT int Polynomial_clear(Polynomial * poly);

EXPORT void Polynomial_init(Polynomial * poly, size_t order, double * coefs, size_t ncoefs);

// makes a copy of the coefs
EXPORT Polynomial * Polynomial_new(size_t order, double * coefs);

EXPORT Polynomial * Polynomial_new_copy(Polynomial * poly, size_t max_order);

EXPORT void Polynomial_copy(Polynomial * dest, Polynomial * src);

// creates internal coefficient references. delete with Polynomial_del_full
EXPORT Polynomial * Polynomial_new_empty(size_t order);

// must delete with Polynomial_del_full
EXPORT Polynomial * Polynomial_new_coefs(size_t order, ...);

EXPORT Polynomial * Polynomial_new_roots(size_t nroots, ...);

EXPORT void Polynomial_del(Polynomial * poly);

EXPORT double Polynomial_eval(Polynomial * poly, double arg);

//need fail statuses
// Polynomial derivative in place
EXPORT int Polynomial_deriv(Polynomial * poly);

// need fail statuses
EXPORT int Polynomial_resize(Polynomial * poly, size_t min_order);

//need fail statuses
// integrate polynomial in place
EXPORT int Polynomial_int(Polynomial * poly, double offset);

EXPORT double Polynomial_defint(Polynomial * poly, double low, double high);

//need fail statuses
// poly_out is modified in place. poly_out->ncoefs must be at least as large as poly_in->order + 1 >= poly_in
EXPORT int Polynomial_add(Polynomial * poly_out, Polynomial * poly_in);

EXPORT int Polynomial_sadd(Polynomial * poly, double offset);

// poly_out is modified in place. poly_out->ncoefs must be at least as large as poly_in->order + 1 >= poly_in
EXPORT int Polynomial_sub(Polynomial * poly_out, Polynomial * poly_in);

EXPORT int Polynomial_ssub(Polynomial * poly, double offset);

// TODO: need failure ints
// poly_out is modified in place. poly_out->ncoefs must be at least as large as poly_in->order * poly_out->order + 1
EXPORT int Polynomial_mul(Polynomial * poly_out, Polynomial * poly_in);

EXPORT int Polynomial_smul(Polynomial * poly, double multiplier);

EXPORT int Polynomial_argmul(Polynomial * poly, size_t arg_order);

EXPORT int Polynomial_arginv(Polynomial * poly);

// TODO: Polynomial_div

EXPORT int Polynomial_sdiv(Polynomial * poly, double multiplier);

EXPORT int Polynomial_scale_domain(Polynomial * poly, Polynomial * scale);

#ifndef __STDC_NO_COMPLEX__
#include <complex.h>

//Polynomial * Polynomial_new_croots(size_t nroots, ...);

EXPORT double complex Polynomial_ceval(Polynomial * poly, double complex arg);

EXPORT double complex Polynomial_cdefint(Polynomial * poly, double complex low, double complex high);

// get complex zeros of polynomial
EXPORT double complex * Polynomial_croots(Polynomial * poly, double tolerance);

EXPORT double complex * Polynomial_croots_coefs(size_t order, double tolerance, ...);

#endif // __STDC_NO_COMPLEX__

#endif // POLYNOMIALS_H