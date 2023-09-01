// follows the "Probabilist's" Hermite polynomials, not the stupid "Physicist's"

#ifndef HERMITE_H
#define HERMITE_H

#include <stddef.h>
#include <polynomials.h>
#ifdef DEBUG_MALLOC
    #include "../../mallocs/debug_malloc.h"
#endif

EXPORT Polynomial * Hermite_new(size_t order);

EXPORT int Hermite_copy(Polynomial * poly, size_t order);

EXPORT Polynomial * Hermite_new_copy(size_t order, size_t max_order);

EXPORT void del_hermite_cache(void);

#endif // HERMITE_H