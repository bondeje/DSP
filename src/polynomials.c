// TODO: remove allocations from root finding

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef M_PI
    #define M_PI		3.14159265358979323846
#endif
#include <stdarg.h>
#include <stdbool.h>
#include <polynomials.h>

#define POLYNOMIAL_ROOT_TOLERANCE 1e-10

// only implements polynomials with real coefficients
// only implements polynomials with positive exponents

// TODO: inline this
double clenshaw_horner_poly_alpha(double val, size_t order) {
    return val;
}

// TODO: inline this
double clenshaw_horner_poly_beta(double val, size_t order) {
    return 0.0;
}

// for Horner's rule (and ultimately calculate all other bases), set any of alpha, beta, base0, or base1 to NULL
// n <= 2 should never happen as this is linear and Polynomial_eval should short-circuit before calling Clenshaw
double clenshaw_algorithm(double x, double * coefs, size_t n, double (*alpha)(double, size_t), double (*beta)(double, size_t), Polynomial * base0, Polynomial * base1) {
    if (!coefs) {
        return NAN;
    }
    if (!alpha || !beta || !base0 || !base1) {
        alpha = clenshaw_horner_poly_alpha;
        beta = clenshaw_horner_poly_beta;
    }
    n--;
    double bnp1 = coefs[n], bnp2 = 0.0;
    n--;
    // terminate when n == 0
    while (n) {
        double S = coefs[n] + alpha(x, n) * bnp1 + beta(x, n + 1) * bnp2;
        bnp2 = bnp1;
        bnp1 = S;
        n--;
    }
    return Polynomial_eval(base0, x) * (coefs[0] + beta(x, 1) * bnp2) + Polynomial_eval(base1, x) * bnp1;
}

void Polynomial_print(Polynomial * poly) {
    if (!poly || !poly->ncoefs) {
        return;
    }
    size_t order = 1;
    printf("order: %zu\n%.4f", poly->order, poly->coefs[0]);
    while (order < poly->ncoefs) {
        printf(" + %.4f*x^%zu", poly->coefs[order], order);
        order++;
    }
    printf("\n");
}

int Polynomial_clear(Polynomial * poly) {
    for (size_t i = 0; i < poly->ncoefs; i++) {
        poly->coefs[i] = 0.0;
    }
    poly->order = 0;
    return 0;
}

void Polynomial_collapse_order(Polynomial * poly) {
    while (poly->order && poly->coefs[poly->order] == 0.0) {
        poly->order--;
    }
}

void Polynomial_init(Polynomial * poly, size_t order, double * coefs, size_t ncoefs) {
    poly->order = order;
    poly->ncoefs = ncoefs;
    poly->coefs = coefs;
    Polynomial_collapse_order(poly);
}

// makes a copy of the coefs
Polynomial * Polynomial_new(size_t order, double * coefs) {
    Polynomial * poly = (Polynomial *) malloc(sizeof(Polynomial));
    if (!poly) {
        return NULL;
    }
    double * new_coefs = (double *) malloc(sizeof(double) * (order + 1));
    if (!new_coefs) {
        POLY_FREE(poly);
        return NULL;
    }
    for (size_t i = 0; i < order + 1; i++) {
        new_coefs[i] = coefs[i];
    }
    Polynomial_init(poly, order, new_coefs, order + 1);
    return poly;
}

Polynomial * Polynomial_new_copy(Polynomial * poly, size_t max_order) {
    Polynomial * copy;
    if (max_order > poly->order) {
        copy = Polynomial_new_empty(max_order);
        Polynomial_add(copy, poly);
    } else {
        copy = Polynomial_new(poly->order, poly->coefs);
    }
    return copy;
}

void Polynomial_copy(Polynomial * dest, Polynomial * src) {
    Polynomial_resize(dest, src->order);
    Polynomial_clear(dest);
    Polynomial_add(dest, src);
}

// creates internal coefficient references. delete with Polynomial_del
Polynomial * Polynomial_new_empty(size_t order) {
    Polynomial * poly = (Polynomial *) malloc(sizeof(Polynomial));
    if (!poly) {
        return NULL;
    }
    double * new_coefs = (double *) malloc(sizeof(double) * (order + 1));
    for (size_t i = 0; i < order + 1; i++) {
        new_coefs[i] = 0.0;
    }
    if (!new_coefs) {
        POLY_FREE(poly);
        return NULL;
    }
    Polynomial_init(poly, order, new_coefs, order + 1);
    Polynomial_clear(poly);
    return poly;
}

// must delete with Polynomial_del
Polynomial * Polynomial_new_coefs(size_t order, ...) {
    Polynomial * poly = (Polynomial *) malloc(sizeof(Polynomial));
    if (!poly) {
        return NULL;
    }
    double * coefs = (double *) malloc(sizeof(double) * (order + 1));
    if (!coefs) {
        POLY_FREE(poly);
        return NULL;
    }
    va_list args;
    va_start(args, order);
    for (size_t i = 0; i < order + 1; i++) {
        coefs[i] = va_arg(args, double);
    }
    va_end(args);

    Polynomial_init(poly, order, coefs, order + 1);
    return poly;
}

Polynomial * Polynomial_new_roots(size_t nroots, ...) { // provide list of roots r_i and return a new polynomial PRODUCT_i=0^nroots-1 (x - r_i)
    Polynomial * poly = Polynomial_new_coefs(0, 1.0); // make constant polynomial
    Polynomial * mult = Polynomial_new_coefs(1, 0.0, 1.0); // a placeholder for linear terms (x - r_i)
    va_list args;
    va_start(args, nroots);
    for (size_t i = 0; i < nroots; i++) {
        mult->coefs[0] = -1.0 * va_arg(args, double);
        Polynomial_mul(poly, mult);
    }
    va_end(args);
    Polynomial_del(mult);
    return poly;
}

void Polynomial_del(Polynomial * poly) {
    Polynomial_clear(poly);
    POLY_FREE(poly->coefs);
    poly->coefs = NULL;
    POLY_FREE(poly);
}

//TODO: once I have Clenshaw's algorithm up and running (https://en.wikipedia.org/wiki/Clenshaw_algorithm), apply it as a generalization of Horner's method
// uses Horner's method
double Polynomial_eval(Polynomial * poly, double arg) {
    if (!poly) {
        return NAN;
    }
    switch (poly->order) {
        case 0:
            return poly->coefs[0];
        case 1:
            return poly->coefs[0] + arg * poly->coefs[1];
        default: // TODO: test clenshaw_algorithm
            // return clenshaw_algorithm(arg, poly->coefs, poly->order+1, NULL, NULL, NULL, NULL);
            break;
    }
    
    size_t n = poly->order;
    double val = poly->coefs[n];
    while (n) {
        n--;
        val = poly->coefs[n] + arg * val;
    }
    return val;
}

//need fail statuses
// Polynomial derivative in place
int Polynomial_deriv(Polynomial * poly) {
    if (!poly) {
        return 1; // failed
    }
    for (size_t n = 0; n < poly->order; n++) {
        poly->coefs[n] = (n+1) * poly->coefs[n+1];
    }
    poly->coefs[poly->order] = 0.0;
    poly->order -= 1;
    return 0;
}

// need fail statuses
int Polynomial_resize(Polynomial * poly, size_t min_order) {
    //printf("in polynomial resize. ncoefs %zu, min order %zu\n", poly->ncoefs, min_order);
    if (poly->ncoefs > min_order) {
        //printf("no-op\n");
        return 0; // no-op
    }
    //printf("I should have returned...why didn't I?\n");
    size_t ncoefs = min_order + 1;
    //printf("resizing array from order %zu to %zu\n", poly->order, min_order);
    double * coefs = (double *) realloc(poly->coefs, sizeof(double)*ncoefs);
    if (!coefs) {
        return 1;
    }
    for (size_t i = poly->ncoefs; i < ncoefs; i++) {
        coefs[i] = 0.0;
    }
    poly->coefs = coefs;
    //printf("old size: %zu\nnew size: %zu\n", poly->ncoefs, ncoefs);
    poly->ncoefs = ncoefs;
    return 0;
}

//need fail statuses
// integrate polynomial in place
int Polynomial_int(Polynomial * poly, double offset) {
    if (!poly) {
        return 1;
    }
    //printf("resizing polynomial\n");
    Polynomial_resize(poly, poly->order + 1);
    //printf("starting integration, %zu %zu\n", poly->ncoefs, poly->order);
    for (size_t n = poly->order+1; n > 0; n--) {
        poly->coefs[n] = poly->coefs[n-1] / n; // don't have to worry about integer division while coefs are floating
    }
    //printf("...integrated\n");
    poly->coefs[0] = offset;
    poly->order++;
    return 0;
}

double Polynomial_defint(Polynomial * poly, double low, double high) {
    double result;
    Polynomial * poly_int = Polynomial_new_copy(poly, poly->order + 1);
    if (Polynomial_int(poly_int, 0.0)) {
        Polynomial_del(poly_int);
        return NAN;
    }
    result = Polynomial_eval(poly_int, high) - Polynomial_eval(poly_int, low);
    Polynomial_del(poly_int);
    return result;
}

// x^(arg_order) * P_n(x); arg_order >= 0
int Polynomial_argmul(Polynomial * poly, size_t arg_order) {
    if (!poly) {
        return 1;
    }
    if (!arg_order) { // do nothing
        return 0;
    }
    Polynomial_resize(poly, poly->order + arg_order);
    size_t i = 0;
    size_t N = poly->order;
    while (i <= N) {
        poly->coefs[N - i + arg_order] = poly->coefs[N - i];
        i++;
    }
    for (size_t j = 0; j < arg_order; j++) {
        poly->coefs[j] = 0.0;
    }
    poly->order += arg_order;
    Polynomial_collapse_order(poly);
    return 0;
}

int Polynomial_arginv(Polynomial * poly) {
    if (!poly) {
        return 1;
    }
    size_t i = 0;
    size_t N = poly->order;
    while (i < N) {
        double temp = poly->coefs[i];
        poly->coefs[i] = poly->coefs[N];
        poly->coefs[N] = temp;
        i++;
        N--;
    }

    return 0;
}

//need fail statuses
// poly_out is modified in place. poly_out->ncoefs must be at least as large as poly_in->order + 1 >= poly_in
int Polynomial_add(Polynomial * poly_out, Polynomial * poly_in) {
    if (!poly_in || !poly_out) {
        return 1;
    }
    Polynomial_resize(poly_out, poly_in->order);
    for (size_t i = 0; i <= poly_in->order; i++) {
        if (i <= poly_out->order) {
            poly_out->coefs[i] += poly_in->coefs[i];
        } else {
            poly_out->coefs[i] = poly_in->coefs[i];
        }
        
    }
    // increase polynomial order if new coefficients added
    poly_out->order = poly_out->order > poly_in->order ? poly_out->order : poly_in->order;
    Polynomial_collapse_order(poly_out);
    return 0;
}

int Polynomial_sadd(Polynomial * poly, double offset) {
    if (!poly) {
        return 1;
    }
    poly->coefs[0] += offset;
    return 0;
}

// poly_out is modified in place. poly_out->ncoefs must be at least as large as poly_in->order + 1 >= poly_in
int Polynomial_sub(Polynomial * poly_out, Polynomial * poly_in) {
    if (!poly_in || !poly_out || poly_out->ncoefs <= poly_in->order) {
        return 1;
    }
    Polynomial_resize(poly_out, poly_in->order);
    for (size_t i = 0; i <= poly_in->order; i++) {
        if (i <= poly_out->order) {
            poly_out->coefs[i] -= poly_in->coefs[i];
        } else {
            poly_out->coefs[i] = -poly_in->coefs[i];
        }
    }
    // increase polynomial order if new coefficients added
    poly_out->order = poly_out->order > poly_in->order ? poly_out->order : poly_in->order;
    // decrease order if maximum polynomial order is zeroed out
    Polynomial_collapse_order(poly_out);
    return 0;
}

int Polynomial_ssub(Polynomial * poly, double offset) {
    return Polynomial_sadd(poly, -offset);
}

// TODO: need failure ints
// poly_out is modified in place. poly_out->ncoefs must be at least as large as poly_in->order * poly_out->order + 1
int Polynomial_mul(Polynomial * poly_out, Polynomial * poly_in) {
    if (!poly_out || !poly_in) {
        return 1;
    }
    //printf("poly_in order %zu, poly_out order %zu, poly_in ncoefs %zu, poly_out ncoefs %zu\n", poly_in->order, poly_out->order, poly_in->ncoefs, poly_out->ncoefs);
    //printf("in Polynomial_mul\n");
    size_t no = poly_out->order;
    size_t ni = poly_in->order;
    size_t n = no + ni;
    if (Polynomial_resize(poly_out, n)) {
        return 1;
    }
    //Polynomial_print(poly_out);
    //Polynomial_print(poly_in);
    poly_out->coefs[n] = poly_out->coefs[no] * poly_in->coefs[ni];
    poly_out->order = n;
    //Polynomial_print(poly_out);
    while (n) {
        n--;
        size_t j = (n >= no) ? 0 : no - n;
        double val = 0.0;
        //printf("%zu %zu %zu %zu \n", n, no, ni, j);
        while (j <= no && n - (no-j) <= ni) {
            val += poly_out->coefs[no - j] * poly_in->coefs[n - (no-j)]; // no-(max-n)+j + ni - j = no - (no+ni-n) + j + ni - j = n
            j++;
            //printf("%f\t", val);
        }
        poly_out->coefs[n] = val;
        //Polynomial_print(poly_out);
    }
    Polynomial_collapse_order(poly_out);
    return 0;
}

int Polynomial_smul(Polynomial * poly, double multiplier) {
    if (!poly) {
        return 1;
    }
    for (size_t i = 0; i <= poly->order; i++) {
        poly->coefs[i] *= multiplier;
    }
    return 0;
}

// TODO: Polynomial_div, but this is somewhat ill-defined. Probably just to quotient, remainder

int Polynomial_sdiv(Polynomial * poly, double multiplier) {
    return Polynomial_smul(poly, 1.0/multiplier);
}

// P_n(x) where x = P_m(y)
// TODO: there's probably a way to do this without the Polynomial_new_copy. Reduce allocations by 2 each time this is called
int Polynomial_scale_domain(Polynomial * poly, Polynomial * scale) {
    if (!poly || !scale) {
        return 1;
    }
    
    size_t old_order = poly->order;
    size_t new_order = old_order * scale->order;
    Polynomial_resize(poly, new_order);
    Polynomial * poly_copy = Polynomial_new_copy(poly, new_order); // copy for the coefficients
    if (!poly_copy) {
        return 1;
    }
    Polynomial_clear(poly);
    size_t i = 0;
    Polynomial_sadd(poly, poly_copy->coefs[old_order-i]);
    i++;
    // Horner's method on polynomials
    while (i <= old_order) {
        Polynomial_mul(poly, scale);
        Polynomial_sadd(poly, poly_copy->coefs[old_order-i]);
        i++;
    }
    Polynomial_del(poly_copy); // clean up poly_copy
    return 0;
}

#ifndef __STDC_NO_COMPLEX__
    #include <complex.h>
    //#include <cmath.h>

/*
Polynomial * Polynomial_new_croots(size_t nroots, ...) { // provide list of roots r_i and return a new polynomial PRODUCT_i=0^nroots-1 (x - r_i)
    Polynomial * poly = Polynomial_new_coefs(0, 1.0); // make constant polynomial
    Polynomial * mult = Polynomial_new_coefs(1, 0.0, 1.0); // a placeholder for linear terms (x - r_i)
    va_list args;
    va_start(args, nroots);
    for (size_t i = 0; i < nroots; i++) {
        mult->coefs[0] = -1.0 * va_arg(args, double);
        Polynomial_mul(poly, mult);
    }
    va_end(args);
    Polynomial_del(mult);
    return poly;
}
*/

//TODO: once I have Clenshaw's algorithm up and running (https://en.wikipedia.org/wiki/Clenshaw_algorithm), apply it as a generalization of Horner's method
// uses Horner's method
double complex Polynomial_ceval(Polynomial * poly, double complex arg) {
    if (!poly) {
        return NAN;
    }
    if (poly->order == 0) {
        return poly->coefs[0];
    }
    
    size_t n = poly->order;
    double complex val = poly->coefs[n];
    while (n) {
        n--;
        val = poly->coefs[n] + arg * val;
    }
    return val;
}

double complex Polynomial_cdefint(Polynomial * poly, double complex low, double complex high) {
    double complex result;
    Polynomial * poly_int = Polynomial_new_copy(poly, poly->order + 1);
    if (Polynomial_int(poly_int, 0.0)) {
        Polynomial_del(poly_int);
        return NAN;
    }
    result = Polynomial_ceval(poly_int, high) - Polynomial_ceval(poly_int, low);
    Polynomial_del(poly_int);
    return result;
}

double get_max_csquared_error(double complex * x, size_t nx, double complex (*cfunc)(void * arg, double complex x), void * arg) {
    double result = 0.0;

    for (size_t i = 0; i < nx; i++) {
        result = fmax(result, cabs(cfunc(arg, x[i])));
    }

    return result;
}

struct pczsr {
    Polynomial * poly;
    double complex * results;
    size_t n;
};

double complex pcz_store_return(void * ptr, double complex value) {
    struct pczsr * p = (struct pczsr *) ptr;
    p->results[p->n] = Polynomial_ceval(p->poly, value);
    return p->results[p->n++];
}

// TODO: replace with something more appropriate
double get_cmax_cabs(double complex * arr, size_t n) {
    double result = 0.0;
    for (size_t i = 0; i < n; i++) {
        result = fmax(result, cabs(arr[i]));
    }
    return result;
}

/* refactoring

// Aberth method
// for a polynomial of order N, the arrays zeros, dw, and dcw must be of length N, N + 1, and N, respectively
int Polynomial_crootsw(double complex * zeros, Polynomial * poly, double tolerance, double * dw, double complex * dcw) {
    if (!poly || !zeros) {
        return 1; // invalid address
    }
    size_t N = poly->order;
    if (!N) {
        return 2; // no roots for order 0 polynomial
    }
    // handle N == 1 or 2 cases especially
    if (N == 1) {
        zeros[0] = -poly->coefs[0] / poly->coefs[1];
        return zeros;
    } else if (N == 2) {
        double complex a = poly->coefs[2];
        double complex b = poly->coefs[1];
        double complex c = poly->coefs[0];
        zeros[0] = (-b - csqrt(b * b - 4 * a * c)) / (2 * a);
        zeros[1] = (-b + csqrt(b * b - 4 * a * c)) / (2 * a);
        return zeros;
    }
    // else iterate
    if (tolerance <= 0.0) {
        tolerance = POLYNOMIAL_ROOT_TOLERANCE;
    }

    double complex p_dp = 0.0 + 0.0*I;

    struct pczsr state = {poly, dcw, 0};

    // initial estimates. based on Alberth's paper's suggestion of centering on -poly->coefs[N-1]/poly->coefs[N]/N
    // TODO: need a better way to initialize the zero estimates
    double center = -poly->coefs[N-1] / poly->coefs[N] / N; // might have divide by zero here
    double radius = fabs(center) < 1.0 ? 1.0 : fabs(center);
    for (size_t i = 0; i < N; i++) {
        zeros[i] = center + radius * cexp(I * M_PI / N * (2.0 * i + 0.5));
        //zeros[i] = cexp(I * M_PI / N * (2.0 * i));
        dcw[i] = 1.0;
    }

    Polynomial poly_der;
    Polynomial_init(&poly_der, N, dw, N + 1);
    Polynomial_copy(&poly_der, poly);
    Polynomial_deriv(&poly_der);

    //Polynomial_print(poly_der);
    size_t iter = 0;
    while (get_max_csquared_error(zeros, N, pcz_store_return, (void*)&state) > tolerance) {
        for (size_t i = 0; i < N; i++) {
            p_dp = state.results[i] / Polynomial_ceval(poly_der, zeros[i]);
            double complex push = 0.0;
            for (size_t j = 0; j < N; j++) {
                push += (j == i) ? 0.0 : 1.0 / (zeros[i] - zeros[j]);
            }
            dcw[i] = p_dp / (1 - p_dp * push);
        }
        for (size_t i = 0; i < N; i++) {
            zeros[i] -= dcw[i];
        }
        state.n = 0;
        iter++;
    }


    // TODO: outline of handling multiple roots
    //      requisite: need to implement the root finding of the companion matrix.
    //      requisite: implement polynomial quotient/remainder division
    // once we reach this point, 
    // 1) separate all the roots that create zeros in poly_der to create two sets: multiple root candidates & simple roots
    // 2) create a polynomial that is the division of 'poly' with the polynomial created from the set of simple roots, ignore the remainder as it is likely just numerical inaccuracy
    // 3) run the companion matrix method to get the multiple roots
    return 0;
}

// for order N polynomial, zeros must have size N
int Polynomial_croots(double complex * zeros, Polynomial * poly, double tolerance) {
    if (!zeros || !poly) {
        return 1;
    }
    size_t N = poly->order;
    if (!N) {
        return 2;
    }
    double * dw = (double *) POLY_MALLOC(sizeof(double) * (N + 1));
    if (!dw) {
        return 3;
    }
    double complex * dcw = (double complex *) POLY_MALLOC(sizeof(double complex) * N);
    if (!dcw) {
        POLY_FREE(dw);
        return 3;
    }
    int status = Polynomial_crootsw(zeros, poly, tolerance, dw, dcw);

    POLY_FREE(dw);
    POLY_FREE(dwc);
    return status;
}

int Polynomial_croots_coefs(double complex * zeros, size_t order, double tolerance, ...) {
    if (!zeros) {
        return 1;
    }
    if (!order) {
        return 2;
    }
    Polynomial * poly = Polynomial_new_empty(order);
    if (!poly) {
        return 3;
    }
    va_list args;
    va_start(args, tolerance);
    for (size_t i = 0; i < order + 1; i++) {
        poly->coefs[i] = va_arg(args, double);
    }
    va_end(args);
    poly->order = order;
    int status = Polynomial_croots(zeros, poly, tolerance);
    Polynomial_del(poly);
    return status;
}

*/

// get complex zeros of polynomial
// todo, need to handle non-simple zeros better. They converge very slowly
// handles poly->order == 1 or 2 analytically
// Aberth method
// TODO: handle zeros at 0 separately since they should be obvious
double complex * Polynomial_croots(Polynomial * poly, double tolerance) {
    if (!poly) {
        printf("poly not valid addres\n");
        return NULL;
    }
    size_t N = poly->order;
    if (!N) {
        printf("order 0 polynomial does not have roots\n");
        return NULL;
    }
    if (tolerance <= 0.0) {
        tolerance = POLYNOMIAL_ROOT_TOLERANCE;
    }

    double complex * zeros = (double complex *) POLY_MALLOC(sizeof(double complex) * poly->order);
    if (!zeros) {
        return NULL;
    }
    
    if (N == 1) {
        zeros[0] = -poly->coefs[0] / poly->coefs[1];
        return zeros;
    } else if (N == 2) {
        double complex a = poly->coefs[2];
        double complex b = poly->coefs[1];
        double complex c = poly->coefs[0];
        zeros[0] = (-b - csqrt(b * b - 4 * a * c)) / (2 * a);
        zeros[1] = (-b + csqrt(b * b - 4 * a * c)) / (2 * a);
        return zeros;
    }
    // else iterate

    double complex * w = (double complex *) POLY_MALLOC(sizeof(double complex) * poly->order);
    if (!w) {
        POLY_FREE(zeros);
        return NULL;
    }

    double complex p_dp = 0.0 + 0.0*I;

    struct pczsr state = {poly, w, 0};

    // initial estimates. based on Alberth's paper's suggestion of centering on -poly->coefs[N-1]/poly->coefs[N]/N
    // TODO: need a better way to initialize the zero estimates
    double center = -poly->coefs[N-1]/poly->coefs[N]/N; // might have divide by zero here
    double radius = fabs(center) < 1.0 ? 1.0 : fabs(center);
    for (size_t i = 0; i < N; i++) {
        zeros[i] = center + radius * cexp(I * M_PI / N * (2.0 * i + 0.5));
        //zeros[i] = cexp(I * M_PI / N * (2.0 * i));
        w[i] = 1.0;
    }

    Polynomial * poly_der = Polynomial_new_copy(poly, poly->order);
    Polynomial_deriv(poly_der);

    //Polynomial_print(poly_der);
    size_t iter = 0;
    while (get_max_csquared_error(zeros, N, pcz_store_return, (void*)&state) > tolerance) {
    //while (get_cmax_cabs(w, N) > tolerance) { // this did not work very well at all
        for (size_t i = 0; i < N; i++) {
            //state.results[i] = Polynomial_ceval(poly, zeros[i]); // only keep in for get_cmax_cabs condition
            p_dp = state.results[i] / Polynomial_ceval(poly_der, zeros[i]);
            double complex push = 0.0;
            for (size_t j = 0; j < N; j++) {
                push += (j == i) ? 0.0 : 1.0 / (zeros[i] - zeros[j]);
            }
            w[i] = p_dp / (1 - p_dp * push);
        }
        for (size_t i = 0; i < N; i++) {
            zeros[i] -= w[i];
        }
        state.n = 0;
        iter++;
    }
    POLY_FREE(w);
    /*
    printf("values at zeros:\n");
    for (size_t i = 0; i < N; i++) {
        printf("%.7f ", cabs(Polynomial_ceval(poly, zeros[i])));
    }
    printf("\nderiviatives at zeros:\n");
    for (size_t i = 0; i < N; i++) {
        printf("%.7f ", cabs(Polynomial_ceval(poly_der, zeros[i])));
    }
    printf("\ntolerance %.7f, sqrt(tolerance) %.7f\n", tolerance, sqrt(tolerance));
    */


    // TODO: outline of handling multiple roots
    //      requisite: need to implement the root finding of the companion matrix.
    //      requisite: implement polynomial quotient/remainder division
    // once we reach this point, 
    // 1) separate all the roots that create zeros in poly_der to create two sets: multiple root candidates & simple roots
    // 2) create a polynomial that is the division of 'poly' with the polynomial created from the set of simple roots, ignore the remainder as it is likely just numerical inaccuracy
    // 3) run the companion matrix method to get the multiple roots





    Polynomial_del(poly_der);
    //printf("found zeros in %zu iterations\n", iter);
    return zeros;
}

double complex * Polynomial_croots_coefs(size_t order, double tolerance, ...) {
    if (!order) {
        return NULL;
    }
    if (tolerance <= 0.0) {
        tolerance = POLYNOMIAL_ROOT_TOLERANCE;
    }
    double complex * zeros = NULL;
    Polynomial * poly = Polynomial_new_empty(order);
    va_list args;
    va_start(args, tolerance);
    for (size_t i = 0; i < order + 1; i++) {
        poly->coefs[i] = va_arg(args, double);
    }
    va_end(args);
    poly->order = order;

    zeros = Polynomial_croots(poly, tolerance);

    Polynomial_del(poly);
    return zeros;
}

#endif // __STDC_NO_COMPLEX__

