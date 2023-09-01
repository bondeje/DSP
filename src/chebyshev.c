#include <stddef.h>
#include <stdlib.h>
#include <polynomials.h>
#include <chebyshev.h>

#define MIN_CHEBYSHEVT_CACHE_SIZE 2

Polynomial ** chebyshevT_cache = NULL;
size_t chebyshevT_cache_size = 0;
size_t chebyshevT_cache_n = 0;

#define MIN_CHEBYSHEVU_CACHE_SIZE 2

Polynomial ** chebyshevU_cache = NULL;
size_t chebyshevU_cache_size = 0;
size_t chebyshevU_cache_n = 0;

// above are global. For each type of Chebyshev polynomial, see below

/* Chebyshev polynomials of the 1st kind. Labeled ChebyshevT */

// forward declaration
Polynomial * ChebyshevT_new_no_cache(size_t order);

int resize_chebyshevT_cache(size_t min_size) {
    Polynomial ** new_chebyshevT_cache = (Polynomial **) POLY_REALLOC(chebyshevT_cache, sizeof(Polynomial*) * min_size);
    if (!new_chebyshevT_cache) {
        return 1;
    }
    for (size_t i = chebyshevT_cache_size; i < min_size; i++) {
        new_chebyshevT_cache[i] = NULL;
    }
    chebyshevT_cache = new_chebyshevT_cache;
    chebyshevT_cache_size = min_size;
    return 0;
}

int init_chebyshevT_cache(size_t size);

int populate_chebyshevT_cache(size_t size) {
    if (!chebyshevT_cache) {
        if (init_chebyshevT_cache(size)) {
            return 1;
        }
    }
    if (chebyshevT_cache_size < size) {
        resize_chebyshevT_cache(size * 2);
    }
    while (chebyshevT_cache_n < size) {
        Polynomial * poly = ChebyshevT_new_no_cache(chebyshevT_cache_n);
        if (!poly) {
            break;
        }
        chebyshevT_cache[chebyshevT_cache_n] = poly;
        chebyshevT_cache_n++;
    }
    return chebyshevT_cache_n < size;
}

// need error codes
int init_chebyshevT_cache(size_t size) {
    if (chebyshevT_cache) {
        return 0;
    }
    if (size < MIN_CHEBYSHEVT_CACHE_SIZE) {
        size = MIN_CHEBYSHEVT_CACHE_SIZE;
    }

    chebyshevT_cache = (Polynomial **) POLY_MALLOC(sizeof(Polynomial*) * size);
    if (!chebyshevT_cache) {
        goto fail1;
    }
    chebyshevT_cache_size = size;

    // order n uses n + 1 coefficients
    double base_coefs[2] = {1.0, 0.0};
    if (!(chebyshevT_cache[0] = Polynomial_new(0, base_coefs))) {
        goto fail2;
    }
    base_coefs[0] = 0.0;
    base_coefs[1] = 1.0;
    if (!(chebyshevT_cache[1] = Polynomial_new(1, base_coefs))) {
        goto fail3;
    }
    chebyshevT_cache_n = 2;
    for (size_t i = chebyshevT_cache_n; i < size; i++) {
        chebyshevT_cache[i] = NULL;
    }
    return 0;

fail3:

    Polynomial_del(chebyshevT_cache[0]);
fail2:

    POLY_FREE(chebyshevT_cache);
fail1:

    chebyshevT_cache = NULL;
    chebyshevT_cache_size = 0;
    chebyshevT_cache_n = 0;
    return 1;
}

void del_chebyshevT_cache(void) {
    if (!chebyshevT_cache) {
        return;
    }
    for (size_t i = 0; i < chebyshevT_cache_n; i++) {
        Polynomial_del(chebyshevT_cache[i]);
        chebyshevT_cache[i] = NULL;
    }
    POLY_FREE(chebyshevT_cache);
    chebyshevT_cache = NULL;
    chebyshevT_cache_n = 0;
    chebyshevT_cache_size = 0;
}

Polynomial * ChebyshevT_new_no_cache(size_t order) {
    
    if (populate_chebyshevT_cache(order)) {
        return NULL;
    }
    if (chebyshevT_cache_n > order) {
        return Polynomial_new_copy(chebyshevT_cache[order], order);
    }
    if ((order + 1 > chebyshevT_cache_size) && resize_chebyshevT_cache(order * 2)) {
        return NULL;
    }
    Polynomial * Tn = Polynomial_new_copy(chebyshevT_cache[order-1], order); // this will eventually be output
    if (!Tn) {
        return NULL;
    }

    // recursion formula
    Polynomial_argmul(Tn, 1); // x * P_n(x)
    Polynomial_smul(Tn, 2);
    Polynomial_sub(Tn, chebyshevT_cache[order-2]);
    return Tn;
}

// copy coefficients from ChebyshevT polynomial of order order into poly without creating a new Polynomial object
int ChebyshevT_copy(Polynomial * poly, size_t order) {
    if (populate_chebyshevT_cache(order + 1)) {
        return 1;
    }
    Polynomial_copy(poly, chebyshevT_cache[order]);
    return 0;
}

Polynomial * ChebyshevT_new_copy(size_t order, size_t max_order) {
    if (populate_chebyshevT_cache(order+1)) {
        return NULL;
    }
    return Polynomial_new_copy(chebyshevT_cache[order], max_order);
}

Polynomial * ChebyshevT_new(size_t order) {
    return ChebyshevT_new_copy(order, order);
}

// ChebyshevT_del is alias for Polynomial_del
void (*ChebyshevT_del)(Polynomial * leg) = Polynomial_del;

/* Chebyshev of the 2nd kind. Labeled ChebyshevU */

// forward declaration
Polynomial * ChebyshevU_new_no_cache(size_t order);

int resize_chebyshevU_cache(size_t min_size) {

    Polynomial ** new_chebyshevU_cache = (Polynomial **) POLY_REALLOC(chebyshevU_cache, sizeof(Polynomial*) * min_size);
    if (!new_chebyshevU_cache) {
        return 1;
    }
    for (size_t i = chebyshevU_cache_size; i < min_size; i++) {
        new_chebyshevU_cache[i] = NULL;
    }
    chebyshevU_cache = new_chebyshevU_cache;
    chebyshevU_cache_size = min_size;
    return 0;
}

int init_chebyshevU_cache(size_t size);

int populate_chebyshevU_cache(size_t size) {
    if (!chebyshevU_cache) {
        if (init_chebyshevU_cache(size)) {
            return 1;
        }
    }
    if (chebyshevU_cache_size < size) {
        resize_chebyshevU_cache(size * 2);
    }
    while (chebyshevU_cache_n < size) {
        Polynomial * poly = ChebyshevU_new_no_cache(chebyshevU_cache_n);
        if (!poly) {
            break;
        }
        chebyshevU_cache[chebyshevU_cache_n] = poly;
        chebyshevU_cache_n++;
    }
    return chebyshevU_cache_n < size;
}

// need error codes
int init_chebyshevU_cache(size_t size) {
    if (chebyshevU_cache) {
        return 0;
    }
    if (size < MIN_CHEBYSHEVU_CACHE_SIZE) {
        size = MIN_CHEBYSHEVU_CACHE_SIZE;
    }
    chebyshevU_cache = (Polynomial **) POLY_MALLOC(sizeof(Polynomial*) * size);
    if (!chebyshevU_cache) {
        goto fail1;
    }
    chebyshevU_cache_size = size;

    // order n uses n + 1 coefficients
    double base_coefs[2] = {1.0, 0.0};
    if (!(chebyshevU_cache[0] = Polynomial_new(0, base_coefs))) {
        goto fail2;
    }
    base_coefs[0] = 0.0;
    base_coefs[1] = 2.0;
    if (!(chebyshevU_cache[1] = Polynomial_new(1, base_coefs))) {
        goto fail3;
    }
    chebyshevU_cache_n = 2;
    for (size_t i = chebyshevU_cache_n; i < size; i++) {
        chebyshevU_cache[i] = NULL;
    }
    return 0;

fail3:
    Polynomial_del(chebyshevU_cache[0]);
fail2:
    POLY_FREE(chebyshevU_cache);
fail1:
    chebyshevU_cache = NULL;
    chebyshevU_cache_size = 0;
    chebyshevU_cache_n = 0;
    return 1;
}

void del_chebyshevU_cache(void) {
    if (!chebyshevU_cache) {
        return;
    }
    for (size_t i = 0; i < chebyshevU_cache_n; i++) {
        Polynomial_del(chebyshevU_cache[i]);
        chebyshevU_cache[i] = NULL;
    }
    POLY_FREE(chebyshevU_cache);
    chebyshevU_cache = NULL;
    chebyshevU_cache_n = 0;
    chebyshevU_cache_size = 0;
}

Polynomial * ChebyshevU_new_no_cache(size_t order) {
    if (populate_chebyshevU_cache(order)) {
        return NULL;
    }
    if (chebyshevU_cache_n > order) {
        return Polynomial_new_copy(chebyshevU_cache[order], order);
    }
    if ((order + 1 > chebyshevU_cache_size) && resize_chebyshevU_cache(order + 1)) {
        return NULL;
    }
    
    Polynomial * Un = Polynomial_new_copy(chebyshevU_cache[order-1], order); // this will eventually be output
    if (!Un) {
        return NULL;
    }

    // Bonnet's recursion formula
    Polynomial_argmul(Un, 1); // x * P_n(x)
    Polynomial_smul(Un, 2.0);
    Polynomial_sub(Un, chebyshevU_cache[order-2]);
    return Un;
}

// copy coefficients from ChebyshevU polynomial of order order into poly without creating a new Polynomial object
int ChebyshevU_copy(Polynomial * poly, size_t order) {
    if (populate_chebyshevU_cache(order + 1)) {
        return 1;
    }
    Polynomial_copy(poly, chebyshevU_cache[order]);
    return 0;
}

Polynomial * ChebyshevU_new_copy(size_t order, size_t max_order) {
    if (populate_chebyshevU_cache(order+1)) {
        return NULL;
    }
    return Polynomial_new_copy(chebyshevU_cache[order], max_order);
}

Polynomial * ChebyshevU_new(size_t order) {
    return ChebyshevU_new_copy(order, order);
}

// ChebyshevU_del is alias for Polynomial_del
void (*ChebyshevU_del)(Polynomial * leg) = Polynomial_del;