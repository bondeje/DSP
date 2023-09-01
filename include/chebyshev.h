#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include <stddef.h>
#include <polynomials.h>
#ifdef DEBUG_MALLOC
    #include "../../mallocs/debug_malloc.h"
#endif

EXPORT Polynomial * ChebyshevT_new(size_t order);

EXPORT int ChebyshevT_copy(Polynomial * poly, size_t order);

EXPORT Polynomial * ChebyshevT_new_copy(size_t order, size_t max_order);

EXPORT extern void (*ChebyshevT_del)(Polynomial * poly);

EXPORT void del_chebyshevT_cache(void);

EXPORT Polynomial * ChebyshevU_new(size_t order);

EXPORT int ChebyshevU_copy(Polynomial * poly, size_t order);

EXPORT Polynomial * ChebyshevU_new_copy(size_t order, size_t max_order);

EXPORT extern void (*ChebyshevU_del)(Polynomial * poly);

EXPORT void del_chebyshevU_cache(void);

#endif // CHEBYSHEV_H