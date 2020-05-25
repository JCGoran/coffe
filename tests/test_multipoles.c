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
#define LEN(arr) ((int) (sizeof (arr) / sizeof (arr)[0]))

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
    const struct coffe_parameters_t *par,
    const struct coffe_background_t *bg,
    const struct coffe_integrals_t *integral
)
{
    /* no errors initially */
    int error_flag = 0;
    /* disabling GSL's stupid error handler */
    gsl_error_handler_t *default_handler =
        gsl_set_error_handler_off();

    const char names[][NAMES_MAXSIZE] = {
        "den", "rsd", "d1", "d2", "g1", "g2", "g3", "g4", "g5", "len", "all"
    };
    const int multipoles[] = {0, 2, 4};

    /* test of individual contributions */
    for (int j = 0; j < LEN(names); ++j){
        const char *type = names[j];
        change_signal(&par->correlation_contrib, type, NULL);
        for (int i = 0; i < LEN(multipoles); ++i){
            const int l = multipoles[i];
            const size_t size_name = 256;
            char name[size_name];
            snprintf(
                name,
                size_name,
                "tests/benchmarks/benchmark_%s_multipoles%d.dat",
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
                    "l = %d, separation = %.3f, type %s\n",
                    l, xvalue[k], type
                );
                weak_assert(
                    approx_equal(y_expected, y_obtained),
                    &error_flag
                );
            }
        }
        reset_signal(&par->correlation_contrib);
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

    #ifdef _OPENMP
    par.nthreads = omp_get_num_procs();
    #endif

    struct coffe_background_t bg;
    coffe_background_init(&par, &bg);

    struct coffe_integrals_t integrals[10];
    coffe_integrals_init(&par, &bg, integrals);

    const int error_flag = coffe_test_multipoles(&par, &bg, integrals);

    coffe_parameters_free(&par);
    coffe_background_free(&bg);
    coffe_integrals_free(integrals);

    return error_flag;
}
