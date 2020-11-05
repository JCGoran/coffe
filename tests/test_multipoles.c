#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <gsl/gsl_errno.h>

#include "common.h"
#include "parser.h"
#include "background.h"
#include "integrals.h"
#include "functions.h"
#include "multipoles.h"
#include "signal.h"
#include "tools.h"

#ifndef NAMES_MAXSIZE
#define NAMES_MAXSIZE 10
#endif

static void change_signal(
    struct coffe_correlation_contributions *signal,
    char *contrib,
    ...
)
{
    va_list args;
    va_start(args, contrib);
    char *current = contrib;

    while (contrib != NULL){
        if (strcmp(contrib, "den") == 0){
            signal->den = 1;
        }
        if (strcmp(contrib, "rsd") == 0){
            signal->rsd = 1;
        }
        if (strcmp(contrib, "d1") == 0){
            signal->d1 = 1;
        }
        if (strcmp(contrib, "d2") == 0){
            signal->d2 = 1;
        }
        if (strcmp(contrib, "g1") == 0){
            signal->g1 = 1;
        }
        if (strcmp(contrib, "g2") == 0){
            signal->g2 = 1;
        }
        if (strcmp(contrib, "g3") == 0){
            signal->g3 = 1;
        }
        if (strcmp(contrib, "len") == 0){
            signal->len = 1;
        }
        if (strcmp(contrib, "g4") == 0){
            signal->g4 = 1;
        }
        if (strcmp(contrib, "g5") == 0){
            signal->g5 = 1;
        }
        /* special keyword "all" */
        if (strcmp(contrib, "all") == 0){
            signal->den = 1;
            signal->rsd = 1;
            signal->d1 = 1;
            signal->d2 = 1;
            signal->g1 = 1;
            signal->g2 = 1;
            signal->g3 = 1;
            signal->g4 = 1;
            signal->g5 = 1;
            signal->len = 1;
        }

        contrib = va_arg(args, char *);
    }
    va_end(args);
}

static void reset_signal(
    struct coffe_correlation_contributions *signal
)
{
    signal->den = 0;
    signal->rsd = 0;
    signal->d1 = 0;
    signal->d2 = 0;
    signal->g1 = 0;
    signal->g2 = 0;
    signal->g3 = 0;
    signal->g4 = 0;
    signal->g5 = 0;
    signal->len = 0;
}

static int coffe_test_multipoles(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integral_array_t *integral
)
{
    /* no errors initially */
    int error_flag = 0;
    /* disabling GSL's stupid error handler */
    gsl_error_handler_t *default_handler =
        gsl_set_error_handler_off();

    char names[][NAMES_MAXSIZE] = {
        "den", "rsd", "d1", "d2", "g1", "g2", "g3", "g4", "g5", "len", "all"
    };
    const int multipoles[] = {0, 2, 4};

    /* test of individual contributions */
    for (size_t j = 0; j < COFFE_ARRAY_SIZE(names); ++j){
        char *type = names[j];
        change_signal(&par->correlation_contrib, type, NULL);
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(multipoles); ++i){
            const int l = multipoles[i];
            const size_t size_name = 256;
            char name[size_name];
            snprintf(
                name,
                size_name,
                DATADIR "/tests/benchmarks/benchmark_%s_multipoles%d.dat",
                type,
                l
            );

            double *xvalue, *yvalue;
            size_t size;
            /* reading the benchmark file */
            coffe_read_ncol(
                name,
                2,
                &size,
                &xvalue, &yvalue
            );

            for (size_t k = 0; k < size; ++k){
                const double x = xvalue[k] * COFFE_H0;
                const double y_expected = yvalue[k];
                const double y_obtained = coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            NONINTEGRATED, MULTIPOLES
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            SINGLE_INTEGRATED, MULTIPOLES
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            DOUBLE_INTEGRATED, MULTIPOLES
                        );

                fprintf(
                    stderr,
                    "l = %d, separation = %.3f, type %s, expected = %e, obtained = %e\n",
                    l, xvalue[k], type, y_expected, y_obtained
                );

                /**
                    some tests have numerical issues due to floating point precision:
                    this is mainly in cases when the integral is actually very close to zero,
                    such as den-den with l >= 2 (it is exactly 0 in the flat-sky case),
                    so we do NOT enforce the check
                **/
                if (!(
                        (
                            (l == 2) && (
                                strcmp(type, "den") == 0 ||
                                strcmp(type, "d2") == 0 ||
                                strcmp(type, "g1") == 0 ||
                                strcmp(type, "g2") == 0 ||
                                strcmp(type, "g3") == 0
                            )
                        ) ||
                        (
                            (l == 4) && (
                                strcmp(type, "den") == 0 ||
                                strcmp(type, "d1") == 0 ||
                                strcmp(type, "d2") == 0 ||
                                strcmp(type, "g1") == 0 ||
                                strcmp(type, "g2") == 0 ||
                                strcmp(type, "g3") == 0 ||
                                strcmp(type, "g4") == 0 ||
                                strcmp(type, "g5") == 0
                            )
                        )
                    )
                ){
                    weak_assert(
                        approx_equal_const_epsilon(y_expected, y_obtained),
                        &error_flag
                    );
                }
            }
        }
        reset_signal(&par->correlation_contrib);
    }

    /* test of some notable special cases */
    {
        /* density + rsd (standard terms) */
        par->correlation_contrib.den = 1;
        par->correlation_contrib.rsd = 1;
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(multipoles); ++i){
            const int l = multipoles[i];
            const size_t size_name = 256;
            char name[size_name];
            snprintf(
                name,
                size_name,
                DATADIR "/tests/benchmarks/benchmark_std_multipoles%d.dat",
                l
            );

            double *xvalue, *yvalue;
            size_t size;
            /* reading the benchmark file */
            coffe_read_ncol(
                name,
                2,
                &size,
                &xvalue, &yvalue
            );

            for (size_t k = 0; k < size; ++k){
                const double x = xvalue[k] * COFFE_H0;
                const double y_expected = yvalue[k];
                const double y_obtained = coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            NONINTEGRATED, MULTIPOLES
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            SINGLE_INTEGRATED, MULTIPOLES
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            DOUBLE_INTEGRATED, MULTIPOLES
                        );

                fprintf(
                    stderr,
                    "l = %d, separation = %.3f, type den+rsd, expected = %e, obtained = %e\n",
                    l, xvalue[k], y_expected, y_obtained
                );

                weak_assert(
                    approx_equal_const_epsilon(y_expected, y_obtained),
                    &error_flag
                );
            }
        }
        reset_signal(&par->correlation_contrib);
    }
    {
        /* density + rsd + lensing */
        par->correlation_contrib.den = 1;
        par->correlation_contrib.rsd = 1;
        par->correlation_contrib.len = 1;
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(multipoles); ++i){
            const int l = multipoles[i];
            const size_t size_name = 256;
            char name[size_name];
            snprintf(
                name,
                size_name,
                DATADIR "/tests/benchmarks/benchmark_std_len_multipoles%d.dat",
                l
            );

            double *xvalue, *yvalue;
            size_t size;
            /* reading the benchmark file */
            coffe_read_ncol(
                name,
                2,
                &size,
                &xvalue, &yvalue
            );

            for (size_t k = 0; k < size; ++k){
                const double x = xvalue[k] * COFFE_H0;
                const double y_expected = yvalue[k];
                const double y_obtained = coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            NONINTEGRATED, MULTIPOLES
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            SINGLE_INTEGRATED, MULTIPOLES
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            DOUBLE_INTEGRATED, MULTIPOLES
                        );

                fprintf(
                    stderr,
                    "l = %d, separation = %.3f, type den+rsd+len, expected = %e, obtained = %e\n",
                    l, xvalue[k], y_expected, y_obtained
                );

                weak_assert(
                    approx_equal_const_epsilon(y_expected, y_obtained),
                    &error_flag
                );
            }
        }
        reset_signal(&par->correlation_contrib);
    }
    {
        /* density + rsd + d1 */
        par->correlation_contrib.den = 1;
        par->correlation_contrib.rsd = 1;
        par->correlation_contrib.d1 = 1;
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(multipoles); ++i){
            const int l = multipoles[i];
            const size_t size_name = 256;
            char name[size_name];
            snprintf(
                name,
                size_name,
                DATADIR "/tests/benchmarks/benchmark_std_d1_multipoles%d.dat",
                l
            );

            double *xvalue, *yvalue;
            size_t size;
            /* reading the benchmark file */
            coffe_read_ncol(
                name,
                2,
                &size,
                &xvalue, &yvalue
            );

            for (size_t k = 0; k < size; ++k){
                const double x = xvalue[k] * COFFE_H0;
                const double y_expected = yvalue[k];
                const double y_obtained = coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            NONINTEGRATED, MULTIPOLES
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            SINGLE_INTEGRATED, MULTIPOLES
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            DOUBLE_INTEGRATED, MULTIPOLES
                        );

                fprintf(
                    stderr,
                    "l = %d, separation = %.3f, type den+rsd+d1, expected = %e, obtained = %e\n",
                    l, xvalue[k], y_expected, y_obtained
                );

                weak_assert(
                    approx_equal_const_epsilon(y_expected, y_obtained),
                    &error_flag
                );
            }
        }
        reset_signal(&par->correlation_contrib);
    }
    {
        /* flatsky lensing-lensing */
        par->correlation_contrib.len = 1;
        par->flatsky_nonlocal = 1;
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(multipoles); ++i){
            const int l = multipoles[i];
            const size_t size_name = 256;
            char name[size_name];
            snprintf(
                name,
                size_name,
                DATADIR "/tests/benchmarks/benchmark_flatsky_lensing_lensing_multipoles%d.dat",
                l
            );

            double *xvalue, *yvalue;
            size_t size;
            /* reading the benchmark file */
            coffe_read_ncol(
                name,
                2,
                &size,
                &xvalue, &yvalue
            );

            for (size_t k = 0; k < size; ++k){
                const double x = xvalue[k] * COFFE_H0;
                const double y_expected = yvalue[k];
                const double y_obtained = coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            NONINTEGRATED, MULTIPOLES
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            SINGLE_INTEGRATED, MULTIPOLES
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            DOUBLE_INTEGRATED, MULTIPOLES
                        );

                fprintf(
                    stderr,
                    "l = %d, separation = %.3f, type lensing (flatsky), expected = %e, obtained = %e\n",
                    l, xvalue[k], y_expected, y_obtained
                );

                weak_assert(
                    approx_equal_const_epsilon(y_expected, y_obtained),
                    &error_flag
                );
            }
        }
        reset_signal(&par->correlation_contrib);
        par->flatsky_nonlocal = 0;
    }
    {
        /* flatsky density-lensing */
        par->correlation_contrib.den = 1;
        par->correlation_contrib.len = 1;
        par->flatsky_local_nonlocal = 1;
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(multipoles); ++i){
            const int l = multipoles[i];
            const size_t size_name = 256;
            char name[size_name];
            snprintf(
                name,
                size_name,
                DATADIR "/tests/benchmarks/benchmark_flatsky_density_lensing_multipoles%d.dat",
                l
            );

            double *xvalue, *yvalue;
            size_t size;
            /* reading the benchmark file */
            coffe_read_ncol(
                name,
                2,
                &size,
                &xvalue, &yvalue
            );

            for (size_t k = 0; k < size; ++k){
                const double x = xvalue[k] * COFFE_H0;
                const double y_expected = yvalue[k];
                const double y_obtained = coffe_integrate(
                            par, bg, integral,
                            x, 0, l,
                            SINGLE_INTEGRATED, MULTIPOLES
                        );

                fprintf(
                    stderr,
                    "l = %d, separation = %.3f, type density-lensing (flatsky), expected = %e, obtained = %e\n",
                    l, xvalue[k], y_expected, y_obtained
                );

                weak_assert(
                    approx_equal_const_epsilon(y_expected, y_obtained),
                    &error_flag
                );
            }
        }
        reset_signal(&par->correlation_contrib);
        par->flatsky_nonlocal = 0;
    }
    {
        /* odd multipoles */
        par->correlation_contrib.den = 1;
        par->correlation_contrib.rsd = 1;
        par->correlation_contrib.d1 = 1;
        par->correlation_contrib.d2 = 1;
        par->correlation_contrib.g1 = 1;
        par->correlation_contrib.g2 = 1;
        par->correlation_contrib.g3 = 1;
        par->correlation_contrib.len = 1;
        par->correlation_contrib.g4 = 1;
        par->correlation_contrib.g5 = 1;
        const int odd_multipoles[] = {1, 3};
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(odd_multipoles); ++i){
            const int l = odd_multipoles[i];
            /* number of separations */
            const size_t size = 15;

            for (size_t k = 0; k < size; ++k){
                const double step = 20.;
                /* dimensionless separation */
                const double x =  step * (k + 1) * COFFE_H0;
                const double y_expected = 0;
                const double y_obtained = coffe_integrate(
                    par, bg, integral,
                    x, 0, l,
                    NONINTEGRATED, MULTIPOLES
                )
                +
                coffe_integrate(
                    par, bg, integral,
                    x, 0, l,
                    SINGLE_INTEGRATED, MULTIPOLES
                )
                +
                coffe_integrate(
                    par, bg, integral,
                    x, 0, l,
                    DOUBLE_INTEGRATED, MULTIPOLES
                );

                fprintf(
                    stderr,
                    "l = %d, separation = %.3f, type all, expected = %e, obtained = %e\n",
                    l, step * (k + 1), y_expected, y_obtained
                );

                weak_assert(
                    approx_equal(
                        y_expected,
                        y_obtained,
                        5e-4,
                        1e-14
                    ),
                    &error_flag
                );
            }
        }
        reset_signal(&par->correlation_contrib);
        par->flatsky_nonlocal = 0;
    }

    if (!error_flag)
        COFFE_TESTS_PRINT_SUCCESS;

    return error_flag;
}

int main(void)
{
    struct coffe_parameters_t par;
    coffe_parse_default_parameters(&par);

    par.divergent = 1;

    #ifdef _OPENMP
    par.nthreads = omp_get_num_procs();
    #endif

    struct coffe_background_t bg;
    coffe_background_init(&par, &bg);

    struct coffe_integral_array_t integrals;
    par.flatsky_nonlocal = 1;
    par.flatsky_local_nonlocal = 1;
    coffe_integrals_init(&par, &bg, &integrals);
    par.flatsky_nonlocal = 0;
    par.flatsky_local_nonlocal = 0;

    const int error_flag = coffe_test_multipoles(&par, &bg, &integrals);

    coffe_parameters_free(&par);
    coffe_background_free(&bg);
    coffe_integrals_free(&integrals);

    return error_flag;
}
