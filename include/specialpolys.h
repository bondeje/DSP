#ifndef SPECIALPOLYS_H
#define SPECIALPOLYS_H

#include <stddef.h>
#include <polynomials.h>
#ifdef DEBUG_MALLOC
    #include "../../mallocs/debug_malloc.h"
#endif

EXPORT Polynomial * Lpoly_new(size_t n); // Optimal "L" polynomial from papoulis for maximum roll-off filtering
extern void (*Lpoly_del)(Polynomial * lp);

#endif