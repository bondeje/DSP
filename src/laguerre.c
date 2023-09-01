#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <polynomials.h>
#include <laguerre.h>

#define MIN_LAGUERRE_CACHE_SIZE 2

Polynomial ** laguerre_cache = NULL;
size_t laguerre_cache_size = 0;
size_t laguerre_cache_n = 0;

// forward declaration
Polynomial * Laguerre_new_no_cache(size_t order);

int resize_laguerre_cache(size_t min_size) {
    Polynomial ** new_laguerre_cache = (Polynomial **) POLY_REALLOC(laguerre_cache, sizeof(Polynomial*) * min_size);
    if (!new_laguerre_cache) {
        return 1;
    }
    for (size_t i = laguerre_cache_size; i < min_size; i++) {
        new_laguerre_cache[i] = NULL;
    }
    laguerre_cache = new_laguerre_cache;
    laguerre_cache_size = min_size;
    return 0;
}

int init_laguerre_cache(size_t size); // forward declaration

int populate_laguerre_cache(size_t size) {
    if (!laguerre_cache) {
        if (init_laguerre_cache(size)) {
            return 1;
        }
    }
    if (laguerre_cache_size < size) {
        resize_laguerre_cache(size * 2);
    }
    while (laguerre_cache_n < size) {
        Polynomial * poly = Laguerre_new_no_cache(laguerre_cache_n);
        if (!poly) {
            break;
        }
        laguerre_cache[laguerre_cache_n] = poly;
        laguerre_cache_n++;
    }
    return laguerre_cache_n < size;
}

// need error codes
int init_laguerre_cache(size_t size) {
    if (laguerre_cache) { // already initialized
        return 0;
    }
    if (size < MIN_LAGUERRE_CACHE_SIZE) {
        size = MIN_LAGUERRE_CACHE_SIZE;
    }
    laguerre_cache = (Polynomial **) POLY_MALLOC(sizeof(Polynomial*) * size);
    if (!laguerre_cache) {
        goto fail1;
    }
    laguerre_cache_size = size;

    // order n uses n + 1 coefficients
    double base_coefs[2] = {1.0, 0.0};
    if (!(laguerre_cache[0] = Polynomial_new(0, base_coefs))) {
        goto fail2;
    }
    base_coefs[0] = 1.0;
    base_coefs[1] = -1.0;
    if (!(laguerre_cache[1] = Polynomial_new(1, base_coefs))) {
        goto fail3;
    }
    laguerre_cache_n = 2;
    for (size_t i = laguerre_cache_n; i < size; i++) {
        laguerre_cache[i] = NULL;
    }
    return 0;

fail3:
    Polynomial_del(laguerre_cache[0]);
fail2:
    POLY_FREE(laguerre_cache);
fail1:
    laguerre_cache = NULL;
    laguerre_cache_size = 0;
    laguerre_cache_n = 0;
    return 1;
}

void del_laguerre_cache(void) {
    if (!laguerre_cache) {
        return;
    }
    for (size_t i = 0; i < laguerre_cache_n; i++) {
        Polynomial_del(laguerre_cache[i]);
        laguerre_cache[i] = NULL;
    }
    POLY_FREE(laguerre_cache);
    laguerre_cache = NULL;
    laguerre_cache_n = 0;
    laguerre_cache_size = 0;
}

Polynomial * Laguerre_new_no_cache(size_t order) {
    if (populate_laguerre_cache(order)) {
        return NULL;
    }
    if (laguerre_cache_n > order) {
        return Polynomial_new_copy(laguerre_cache[order], order);
    }
    if ((order + 1 > laguerre_cache_size) && resize_laguerre_cache(order + 1)) {
        return NULL;
    }
    Polynomial * Ln = Polynomial_new_copy(laguerre_cache[order-1], order); // this will eventually be output
    if (!Ln) {
        return NULL;
    }
    double coefs[] = {2 * order - 1.0, -1.0};
    Polynomial alpha;
    Polynomial_init(&alpha, 1, coefs, 2);
    Polynomial_smul(Ln, 1.0/(order - 1.0));
    Polynomial_mul(Ln, &alpha);
    Polynomial_sub(Ln, laguerre_cache[order - 2]);
    Polynomial_smul(Ln, (order - 1.0) / order);
    return Ln;
}

// copy coefficients from Laguerre polynomial of order order into poly without creating a new Polynomial object
int Laguerre_copy(Polynomial * poly, size_t order) {
    if (populate_laguerre_cache(order + 1)) {
        return 1;
    }
    Polynomial_copy(poly, laguerre_cache[order]);
    return 0;
}

Polynomial * Laguerre_new_copy(size_t order, size_t max_order) {
    if (populate_laguerre_cache(order+1)) {
        return NULL;
    }
    return Polynomial_new_copy(laguerre_cache[order], max_order);
}

Polynomial * Laguerre_new(size_t order) {
    return Laguerre_new_copy(order, order);
}

// Laguerre_del is alias for Polynomial_del
void (*Laguerre_del)(Polynomial * leg) = Polynomial_del;