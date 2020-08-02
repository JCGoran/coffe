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
#include "corrfunc.h"
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

static int coffe_test_corrfunc(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t *integral
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
    const double mu[] = {0.0, 0.2, 0.5, 0.8, 0.95};

    /* test of individual contributions */
    for (size_t j = 0; j < COFFE_ARRAY_SIZE(names); ++j){
        char *type = names[j];
        change_signal(&par->correlation_contrib, type, NULL);
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(mu); ++i){
            const double m = mu[i];
            const size_t size_name = 256;
            char name[size_name];
            snprintf(
                name,
                size_name,
                "tests/benchmarks/benchmark_%s_corrfunc%zu.dat",
                type,
                i
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
                            x, m, 0,
                            NONINTEGRATED, CORRFUNC
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, m, 0,
                            SINGLE_INTEGRATED, CORRFUNC
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, m, 0,
                            DOUBLE_INTEGRATED, CORRFUNC
                        );

                fprintf(
                    stderr,
                    "mu = %.3f, separation = %.3f, type %s, expected = %e, obtained = %e\n",
                    m, xvalue[k], type, y_expected, y_obtained
                );

                weak_assert(
                    approx_equal(y_expected, y_obtained),
                    &error_flag
                );
            }
        }
        reset_signal(&par->correlation_contrib);
    }
    /* test of some notable special cases */
    {
        /* density + rsd (standard terms) */
        par->correlation_contrib.den = 1;
        par->correlation_contrib.rsd = 1;
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(mu); ++i){
            const double m = mu[i];
            const size_t size_name = 256;
            char name[size_name];
            snprintf(
                name,
                size_name,
                "tests/benchmarks/benchmark_std_corrfunc%zu.dat",
                i
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
                            x, m, 0,
                            NONINTEGRATED, CORRFUNC
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, m, 0,
                            SINGLE_INTEGRATED, CORRFUNC
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, m, 0,
                            DOUBLE_INTEGRATED, CORRFUNC
                        );

                fprintf(
                    stderr,
                    "mu = %.3f, separation = %.3f, type den+rsd, expected = %e, obtained = %e\n",
                    m, xvalue[k], y_expected, y_obtained
                );

                weak_assert(
                    approx_equal(y_expected, y_obtained),
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
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(mu); ++i){
            const double m = mu[i];
            const size_t size_name = 256;
            char name[size_name];
            snprintf(
                name,
                size_name,
                "tests/benchmarks/benchmark_std_len_corrfunc%zu.dat",
                i
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
                            x, m, 0,
                            NONINTEGRATED, CORRFUNC
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, m, 0,
                            SINGLE_INTEGRATED, CORRFUNC
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, m, 0,
                            DOUBLE_INTEGRATED, CORRFUNC
                        );

                fprintf(
                    stderr,
                    "mu = %.3f, separation = %.3f, type den+rsd+len, expected = %e, obtained = %e\n",
                    m, xvalue[k], y_expected, y_obtained
                );

                weak_assert(
                    approx_equal(y_expected, y_obtained),
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
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(mu); ++i){
            const double m = mu[i];
            const size_t size_name = 256;
            char name[size_name];
            snprintf(
                name,
                size_name,
                "tests/benchmarks/benchmark_std_d1_corrfunc%zu.dat",
                i
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
                            x, m, 0,
                            NONINTEGRATED, CORRFUNC
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, m, 0,
                            SINGLE_INTEGRATED, CORRFUNC
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, m, 0,
                            DOUBLE_INTEGRATED, CORRFUNC
                        );

                fprintf(
                    stderr,
                    "mu = %.3f, separation = %.3f, type den+rsd+d1, expected = %e, obtained = %e\n",
                    m, xvalue[k], y_expected, y_obtained
                );

                weak_assert(
                    approx_equal(y_expected, y_obtained),
                    &error_flag
                );
            }
        }
        reset_signal(&par->correlation_contrib);
    }
    {
        /* flatsky lensing-lensing */
        par->correlation_contrib.len = 1;
        par->flatsky_lensing_lensing = 1;
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(mu); ++i){
            const double m = mu[i];
            const size_t size_name = 256;
            char name[size_name];
            snprintf(
                name,
                size_name,
                "tests/benchmarks/benchmark_flatsky_lensing_lensing_corrfunc%zu.dat",
                i
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
                            x, m, 0,
                            NONINTEGRATED, CORRFUNC
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, m, 0,
                            SINGLE_INTEGRATED, CORRFUNC
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, m, 0,
                            DOUBLE_INTEGRATED, CORRFUNC
                        );

                fprintf(
                    stderr,
                    "mu = %.3f, separation = %.3f, type len (flatsky), expected = %e, obtained = %e\n",
                    m, xvalue[k], y_expected, y_obtained
                );

                weak_assert(
                    approx_equal(y_expected, y_obtained),
                    &error_flag
                );
            }
        }
        reset_signal(&par->correlation_contrib);
        par->flatsky_lensing_lensing = 0;
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
    par.nonzero_terms[8].n = 4, par.nonzero_terms[8].l = 0;

    {
        free(par.mu);
        const double mu[] = {0.0, 0.2, 0.5, 0.8, 0.95};
        par.mu = coffe_malloc(sizeof(double) * COFFE_ARRAY_SIZE(mu));
        par.mu_len = COFFE_ARRAY_SIZE(mu);
        for (size_t i = 0; i < COFFE_ARRAY_SIZE(mu); ++i)
            par.mu[i] = mu[i];
    }

    #ifdef _OPENMP
    par.nthreads = omp_get_num_procs();
    #endif

    struct coffe_background_t bg;
    coffe_background_init(&par, &bg);

    struct coffe_integrals_t integrals[10];
    par.flatsky_lensing_lensing = 1;
    coffe_integrals_init(&par, &bg, integrals);
    par.flatsky_lensing_lensing = 0;

    const int error_flag = coffe_test_corrfunc(&par, &bg, integrals);

    coffe_parameters_free(&par);
    coffe_background_free(&bg);
    coffe_integrals_free(integrals);

    return error_flag;
}
