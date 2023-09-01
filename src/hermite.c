#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <polynomials.h>
#include <hermite.h>

// TODO: seriously need to rethink how the cache is maintained in here and other Polynomials

#define MIN_HERMITE_CACHE_SIZE 2

Polynomial ** hermite_cache = NULL;
size_t hermite_cache_size = 0;
size_t hermite_cache_n = 0;

// forward declaration
Polynomial * Hermite_new_no_cache(size_t order);

int resize_hermite_cache(size_t min_size) {
    Polynomial ** new_hermite_cache = (Polynomial **) POLY_REALLOC(hermite_cache, sizeof(Polynomial*) * min_size);
    if (!new_hermite_cache) {
        return 1;
    }
    for (size_t i = hermite_cache_size; i < min_size; i++) {
        new_hermite_cache[i] = NULL;
    }
    hermite_cache = new_hermite_cache;
    hermite_cache_size = min_size;
    return 0;
}

int init_hermite_cache(size_t size); // forward declaration

int populate_hermite_cache(size_t size) {
    if (!hermite_cache) {
        if (init_hermite_cache(size)) {
            return 1;
        }
    }
    if (hermite_cache_size < size) {
        resize_hermite_cache(size * 2);
    }
    while (hermite_cache_n < size) {
        Polynomial * poly = Hermite_new_no_cache(hermite_cache_n);
        if (!poly) {
            break;
        }
        hermite_cache[hermite_cache_n] = poly;
        hermite_cache_n++;
    }
    return hermite_cache_n < size;
}

// need error codes
int init_hermite_cache(size_t size) {
    if (hermite_cache) {
        return 0;
    }
    if (size < MIN_HERMITE_CACHE_SIZE) {
        size = MIN_HERMITE_CACHE_SIZE;
    }
    hermite_cache = (Polynomial **) POLY_MALLOC(sizeof(Polynomial*) * size);
    if (!hermite_cache) {
        goto fail1;
    }
    hermite_cache_size = size;

    // order n uses n + 1 coefficients
    double base_coefs[2] = {1.0, 0.0};
    if (!(hermite_cache[0] = Polynomial_new(0, base_coefs))) {
        goto fail2;
    }
    base_coefs[0] = 0.0;
    base_coefs[1] = 1.0;
    if (!(hermite_cache[1] = Polynomial_new(1, base_coefs))) {
        goto fail3;
    }
    hermite_cache_n = 2;
    for (size_t i = hermite_cache_n; i < size; i++) {
        hermite_cache[i] = NULL;
    }
    return 0;

fail3:
    Polynomial_del(hermite_cache[0]);
fail2:
    POLY_FREE(hermite_cache);
fail1:
    hermite_cache = NULL;
    hermite_cache_size = 0;
    hermite_cache_n = 0;
    return 1;
}

void del_hermite_cache(void) {
    //printf("deleting cache\n");
    if (!hermite_cache) {
        return;
    }
    for (size_t i = 0; i < hermite_cache_n; i++) {
        Polynomial_del(hermite_cache[i]);
        hermite_cache[i] = NULL;
    }
    POLY_FREE(hermite_cache);
    hermite_cache = NULL;
    hermite_cache_n = 0;
    hermite_cache_size = 0;
}

Polynomial * Hermite_new_no_cache(size_t order) {
    if (populate_hermite_cache(order)) {
        return NULL;
    }
    if (hermite_cache_n > order) {
        return Polynomial_new_copy(hermite_cache[order], order);
    }
    if ((order + 1 > hermite_cache_size) && resize_hermite_cache(order * 2)) {
        return NULL;
    }
    Polynomial * Hn = Polynomial_new_copy(hermite_cache[order-1], order); // this will eventually be output
    if (!Hn) {
        return NULL;
    }
    Polynomial_argmul(Hn, 1);
    Polynomial_smul(Hn, 1.0 / (order - 1.0));
    Polynomial_sub(Hn, hermite_cache[order - 2]);
    Polynomial_smul(Hn, order - 1.0);
    return Hn;
}

// copy coefficients from Hermite polynomial of order order into poly without creating a new Polynomial object
int Hermite_copy(Polynomial * poly, size_t order) {
    if (populate_hermite_cache(order + 1)) {
        return 1;
    }
    Polynomial_copy(poly, hermite_cache[order]);
    return 0;
}

Polynomial * Hermite_new_copy(size_t order, size_t max_order) {
    if (populate_hermite_cache(order+1)) {
        return NULL;
    }
    return Polynomial_new_copy(hermite_cache[order], max_order);
}

Polynomial * Hermite_new(size_t order) {
    return Hermite_new_copy(order, order);
}

// Hermite_del is alias for Polynomial_del
void (*Hermite_del)(Polynomial * leg) = Polynomial_del;