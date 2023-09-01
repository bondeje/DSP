#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <polynomials.h>
#include <legendre.h>
#include <specialpolys.h>

Polynomial * Lpoly_new(size_t n) {
    if (!n) {
        return NULL;
    }
    
    Polynomial * Ln = Polynomial_new_empty(n * 2);
    if (!Ln) {
        return NULL;
    }
    
    Polynomial * builder = Polynomial_new_empty(n * 2);
    if (!builder) {
        Polynomial_del(Ln);
        return NULL;
    }
    
    if (n & 1) { // follow odd 
        size_t k = n / 2;
        for (size_t i = 0; i <= k; i++) {
            printf("building with Legendre %zu\n", i);
            if (Legendre_copy(builder, i)) {
                goto fail_cleanup;
            }
            Polynomial_smul(builder, (2.0 * i + 1.0) / (k + 1.0));
            Polynomial_add(Ln, builder);
        }
        Polynomial_copy(builder, Ln); // make builder == Ln
        Polynomial_mul(Ln, builder); // square the sum
        if (Polynomial_int(Ln, 0.0)) {
            goto fail_cleanup;
        }
        Polynomial_ssub(Ln, Polynomial_eval(Ln, -1));

        // create scaling polynomial
        Polynomial_clear(builder);
        Polynomial_sadd(builder, 2.0);
        Polynomial_argmul(builder, 2);
        Polynomial_sadd(builder, -1.0);
        Polynomial_scale_domain(Ln, builder);
        Polynomial_sdiv(Ln, 2.0); // the stupid square root in the sum of Legendre polynomials that gets squared. Seriously, you have a polynomial with integer coefficients and you stick an unnecessary square root in the multiplier of each one?
    } else { // n even
        size_t k = n / 2 - 1;
        for (size_t i = (k & 1); i <= k; i+=2) {
            if (Legendre_copy(builder, i)) {
                goto fail_cleanup;
            }
            Polynomial_smul(builder, (2.0 * i + 1.0));
            Polynomial_add(Ln, builder);
        }
        Polynomial_copy(builder, Ln); // make builder == Ln
        Polynomial_mul(Ln, builder); // square the sum
        // maybe have a constant for (x+1)
        Polynomial_clear(builder);
        Polynomial_sadd(builder, 1.0);
        Polynomial_argmul(builder, 1);
        Polynomial_sadd(builder, 1.0);
        Polynomial_mul(Ln, builder); // (x+1) * [sum_i_k a_i P_i]^2
        if (Polynomial_int(Ln, 0.0)) {
            goto fail_cleanup;
        }
        Polynomial_ssub(Ln, Polynomial_eval(Ln, -1));

        // create scaling polynomial
        Polynomial_clear(builder);
        Polynomial_sadd(builder, 2.0);
        Polynomial_argmul(builder, 2);
        Polynomial_sadd(builder, -1.0);

        Polynomial_scale_domain(Ln, builder);
        Polynomial_sdiv(Ln, (k + 1)*(k + 2)); // the stupid square root in the sum of Legendre polynomials that gets squared. Seriously, you have a polynomial with integer coefficients and you stick an unnecessary square root in the multiplier of each one?
    }

    Polynomial_del(builder);
    
    return Ln;

fail_cleanup:
    Polynomial_del(builder);
    Polynomial_del(Ln);
    return NULL;
}

void (*Lpoly_del)(Polynomial * lp) = Polynomial_del;