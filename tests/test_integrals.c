#include <math.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <gsl/gsl_errno.h>

#include "common.h"
#include "parser.h"
#include "background.h"
#include "integrals.h"
#include "tools.h"

/* comparing anything above this value (in Mpc) is not sensible */
#ifndef MAX_SEPARATION
#define MAX_SEPARATION 10000
#endif

#ifndef COFFE_H0
#define COFFE_H0 (1. / 2997.92458)
#endif
#define h (0.67)

static int coffe_test_integrals(
    const coffe_integral_array_t *integrals
)
{
    int error_flag = 0;
    const size_t size_files = 8;
    size_t size;
    double *xvalue, *yvalue[size_files];

    gsl_error_handler_t *default_handler =
        gsl_set_error_handler_off();

    /* load the files */
    for (size_t i = 0; i < size_files; ++i){
        const size_t size_name = 256;
        char name[size_name];
        snprintf(
            name,
            size_name,
            COFFE_TEST_DATADIR "/tests/benchmarks/benchmark_integral%zu.dat",
            i
        );
        coffe_read_ncol(
            name,
            2, &size,
            &xvalue,
            &yvalue[i]
        );
    }

    const struct nl_terms terms[] = {
        {.n = 0, .l = 0},
        {.n = 0, .l = 2},
        {.n = 0, .l = 4},
        {.n = 1, .l = 1},
        {.n = 1, .l = 3},
        {.n = 2, .l = 0},
        {.n = 2, .l = 2},
        {.n = 3, .l = 1}
    };

    assert(size_files == sizeof(terms) / sizeof(*terms));

    /* compare the standard integrals (ones computes with 2FAST) */
    for (size_t integral = 0; integral < size_files; ++integral){
        const int n = terms[integral].n, l = terms[integral].l;
        for (size_t i = 0; i < size - 1; ++i){
            const double x = xvalue[i] / (COFFE_H0 * h);
            const double y_expected = yvalue[integral][i];
            double y_obtained = coffe_interp_spline(
                &coffe_find_integral(
                    integrals,
                    n,
                    l,
                    COFFE_INTEGER,
                    COFFE_INTEGER
                )->result,
                x
            );
            if (n > l){
                y_obtained *= pow(COFFE_H0 * h, n - l);
            }
            fprintf(
                stderr,
                "Integral = %zu, n = %d, l = %d, separation = %e, expected = %e, obtained = %e\n",
                integral, n, l, x, y_expected, y_obtained
            );
            /* only compare until MAX_SEPARATION, as above it doesn't matter */
            if (x > 1.5 && x < MAX_SEPARATION){
                weak_assert(
                    approx_equal(y_expected, y_obtained, 1e-2, 0),
                    &error_flag
                );
            }
        }
    }

    /* compare the standard integrals (in flatsky) */
    {
        const size_t size_name = 256;
        char name[size_name];
        snprintf(
            name,
            size_name,
            COFFE_TEST_DATADIR "/tests/benchmarks/benchmark_integral%zu.dat",
            (size_t)9
        );
        double *x_array, *y_array;
        coffe_read_ncol(
            name,
            2,
            &size,
            &x_array,
            &y_array
        );

        for (size_t i = 0; i < size - 1; ++i){
            const double x = x_array[i] / (COFFE_H0 * h);
            const double y_expected = y_array[i];
            const double y_obtained = coffe_interp_spline(
                &coffe_find_integral(
                    integrals,
                    1,
                    -1,
                    COFFE_HALF_INTEGER,
                    COFFE_HALF_INTEGER
                )->result,
                x
            ) * pow(COFFE_H0 * h, 1);

            fprintf(
                stderr,
                "Integral = %zu, n = %d, l = %d, separation = %e, expected = %e, obtained = %e\n",
                (size_t)9, 1, -1, x, y_expected, y_obtained
            );
            /* only compare until MAX_SEPARATION, as above it doesn't matter */
            if (x > 2 && x < 5000){
                weak_assert(
                    approx_equal(y_expected, y_obtained, 1e-3, 0),
                    &error_flag
                );
            }
        }
        free(x_array);
        free(y_array);
    }

    /* memory cleanup */
    free(xvalue);
    for (size_t i = 0; i < size_files; ++i)
        free(yvalue[i]);

    /* the divergent integral isn't computed with 2FAST */
    double *divergent_x, *divergent_y;
    size_t divergent_size;

    coffe_read_ncol(
        COFFE_TEST_DATADIR "/tests/benchmarks/benchmark_integral8.dat",
        2,
        &divergent_size,
        &divergent_x,
        &divergent_y
    );

    for (size_t i = 0; i < divergent_size - 1; ++i){

        const double y_expected = divergent_y[i];
        const double y_obtained = coffe_interp_spline(
            &coffe_find_integral(
                integrals,
                4,
                0,
                COFFE_INTEGER,
                COFFE_INTEGER
            )->result,
            divergent_x[i] / (COFFE_H0 * h)
        ) * pow(COFFE_H0 * h, 4);

        fprintf(
            stderr,
            "Integral 8, separation = %e, expected = %e, obtained = %e\n",
            divergent_x[i], y_expected, y_obtained
        );

        weak_assert(
            approx_equal(
                y_expected,
                y_obtained,
                1e-3,
                0
            ),
            &error_flag
        );
    }

    /* memory cleanup */
    free(divergent_x);
    free(divergent_y);

    /* test renormalization */
    double *ren_x, *ren_y, *ren_z;
    size_t ren_size;

    coffe_read_ncol(
        COFFE_TEST_DATADIR "/tests/benchmarks/benchmark_integral8_renormalization.dat",
        3,
        &ren_size,
        &ren_x,
        &ren_y,
        &ren_z
    );

    for (size_t i = 0; i < ren_size; ++i){

        const double y_expected = ren_z[i];
        const double y_obtained = coffe_interp_spline2d(
            &coffe_find_integral(
                integrals,
                4,
                0,
                COFFE_INTEGER,
                COFFE_INTEGER
            )->renormalization,
            ren_x[i] / (COFFE_H0 * h),
            ren_y[i] / (COFFE_H0 * h)
        ) * pow(COFFE_H0 * h, 4);

        fprintf(
            stderr,
            "Integral 8 renormalization, x1 = %e, x2 = %e, expected = %e, obtained = %e\n",
            ren_x[i] / (COFFE_H0 * h), ren_y[i] / (COFFE_H0 * h), y_expected, y_obtained
        );

        weak_assert(
            approx_equal(
                y_expected,
                y_obtained,
                1e-3,
                0
            ),
            &error_flag
        );
    }

    /* memory cleanup */
    free(ren_x);
    free(ren_y);
    free(ren_z);

    if (!error_flag)
        COFFE_TESTS_PRINT_SUCCESS;

    gsl_set_error_handler(default_handler);

    return error_flag;
}

int main(void)
{
    coffe_parameters_t par;
    coffe_parse_default_parameters(&par);

    par.divergent = 1;
    par.flatsky_nonlocal = 1;

    coffe_background_t bg = {.flag = 0};
    coffe_background_init(&par, &bg);

    coffe_integral_array_t integrals = {.array = NULL, .size = 0};
    coffe_integrals_init(&par, &bg, &integrals);

    const int error_flag = coffe_test_integrals(&integrals);

    coffe_parameters_free(&par);
    coffe_background_free(&bg);
    coffe_integrals_free(&integrals);

    return error_flag;
}
