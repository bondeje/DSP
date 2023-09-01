#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <polynomials.h>
#include <legendre.h>

#define MIN_LEGENDRE_CACHE_SIZE 2

Polynomial ** legendre_cache = NULL;
size_t legendre_cache_size = 0;
size_t legendre_cache_n = 0;

// forward declaration
Polynomial * Legendre_new_no_cache(size_t order);

int resize_legendre_cache(size_t min_size) {
    Polynomial ** new_legendre_cache = (Polynomial **) POLY_REALLOC(legendre_cache, sizeof(Polynomial*) * min_size);
    if (!new_legendre_cache) {
        return 1;
    }
    for (size_t i = legendre_cache_size; i < min_size; i++) {
        new_legendre_cache[i] = NULL;
    }
    //printf("reallocated legendre cache from %zu to %zu\n", legendre_cache_size, min_size);
    legendre_cache = new_legendre_cache;
    legendre_cache_size = min_size;
    return 0;
}

int init_legendre_cache(size_t size); // forward declaration

int populate_legendre_cache(size_t size) {
    if (!legendre_cache) {
        return init_legendre_cache(size);
    }
    if (legendre_cache_size < size) {
        resize_legendre_cache(size * 2);
    }
    while (legendre_cache_n < size) {
        Polynomial * poly = Legendre_new_no_cache(legendre_cache_n);
        if (!poly) {
            break;
        }
        legendre_cache[legendre_cache_n] = poly;
        legendre_cache_n++;
    }
    //printf("\n... to %zu", legendre_cache_n);
    return legendre_cache_n < size;
}

// need error codes
int init_legendre_cache(size_t size) {
    if (legendre_cache) {
        return 0;
    }
    if (size < MIN_LEGENDRE_CACHE_SIZE) {
        size = MIN_LEGENDRE_CACHE_SIZE;
    }
    legendre_cache = (Polynomial **) POLY_MALLOC(sizeof(Polynomial*) * size);
    if (!legendre_cache) {
        goto fail1;
    }
    legendre_cache_size = size;

    // order n uses n + 1 coefficients
    double base_coefs[2] = {1.0, 0.0};
    if (!(legendre_cache[0] = Polynomial_new(0, base_coefs))) {
        goto fail2;
    }
    base_coefs[0] = 0.0;
    base_coefs[1] = 1.0;
    if (!(legendre_cache[1] = Polynomial_new(1, base_coefs))) {
        goto fail3;
    }
    legendre_cache_n = 2;
    for (size_t i = legendre_cache_n; i < size; i++) {
        legendre_cache[i] = NULL;
    }
    return 0;

fail3:
    Polynomial_del(legendre_cache[0]);
fail2:
    POLY_FREE(legendre_cache);
fail1:
    legendre_cache = NULL;
    legendre_cache_size = 0;
    legendre_cache_n = 0;
    return 1;
}

void del_legendre_cache(void) {
    if (!legendre_cache) {
        return;
    }
    //printf("deleting legendre cache of size %zu\n", legendre_cache_n);
    for (size_t i = 0; i < legendre_cache_n; i++) {
        //printf("deleting legendre polynomial at %p\n", (void*)legendre_cache[i]);
        Polynomial_del(legendre_cache[i]);
        legendre_cache[i] = NULL; // this gave a valgrind error
    }
    POLY_FREE(legendre_cache);
    legendre_cache = NULL;
    legendre_cache_n = 0;
    legendre_cache_size = 0;
}

Polynomial * Legendre_new_no_cache(size_t order) {
    if (populate_legendre_cache(order)) {
        return NULL;
    }
    if (legendre_cache_n > order) {
        return Polynomial_new_copy(legendre_cache[order], order);
    }
    if ((order + 1 > legendre_cache_size) && resize_legendre_cache(order + 1)) {
        return NULL;
    }
    Polynomial * Pn = Polynomial_new_copy(legendre_cache[order-1], order); // this will eventually be output
    if (!Pn) {
        return NULL;
    }
    // Bonnet's recursion formula
    Polynomial_argmul(Pn, 1); // x * P_n(x)
    Polynomial_smul(Pn, (2.0 * order - 1.0) / (order - 1.0));
    Polynomial_sub(Pn, legendre_cache[order - 2]);
    Polynomial_smul(Pn, (order - 1.0) / order);
    return Pn;
}

// copy coefficients from Legendre polynomial of order order into poly without creating a new Polynomial object
int Legendre_copy(Polynomial * poly, size_t order) {
    if (populate_legendre_cache(order + 1)) {
        return 1;
    }
    Polynomial_copy(poly, legendre_cache[order]);
    return 0;
}

Polynomial * Legendre_new_copy(size_t order, size_t max_order) {
    if (populate_legendre_cache(order+1)) {
        return NULL;
    }
    return Polynomial_new_copy(legendre_cache[order], max_order);
}

Polynomial * Legendre_new(size_t order) {
    return Legendre_new_copy(order, order);
}

// Legendre_del is alias for Polynomial_del
void (*Legendre_del)(Polynomial * leg) = Polynomial_del;