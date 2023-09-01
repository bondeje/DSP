#ifndef LAGUERRE_H
#define LAGUERRE_H

#include <stddef.h>
#include <polynomials.h>
#ifdef DEBUG_MALLOC
    #include "../../mallocs/debug_malloc.h"
#endif

EXPORT Polynomial * Laguerre_new(size_t order);

EXPORT int Laguerre_copy(Polynomial * poly, size_t order);

EXPORT Polynomial * Laguerre_new_copy(size_t order, size_t max_order);

EXPORT void del_laguerre_cache(void);

#endif // LAGUERRE_H