#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <stddef.h>
#include <polynomials.h>
#ifdef DEBUG_MALLOC
    #include "../../mallocs/debug_malloc.h"
#endif

EXPORT Polynomial * Legendre_new(size_t order);

EXPORT int Legendre_copy(Polynomial * poly, size_t order);

EXPORT Polynomial * Legendre_new_copy(size_t order, size_t max_order);

EXPORT void del_legendre_cache(void);

#endif // LEGENDRE_H