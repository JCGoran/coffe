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

    /* compiling with benchmarked values */
    const double xvalue[] = {
        10, 20, 40, 100, 150
    };
    const char names[][NAMES_MAXSIZE] = {
        "den", "rsd", "d1", "d2", "g1", "g2", "g3", "g4", "g5", "len"
    };
    const int multipoles[] = {0, 2, 4};
    const double yvalue[][sizeof(xvalue) / sizeof(*xvalue)] = {
    {
        #include "benchmark_multipoles0_den.dat"
    },
    {
        #include "benchmark_multipoles2_den.dat"
    },
    {
        #include "benchmark_multipoles4_den.dat"
    },
    {
        #include "benchmark_multipoles0_rsd.dat"
    },
    {
        #include "benchmark_multipoles2_rsd.dat"
    },
    {
        #include "benchmark_multipoles4_rsd.dat"
    },
    {
        #include "benchmark_multipoles0_d1.dat"
    },
    {
        #include "benchmark_multipoles2_d1.dat"
    },
    {
        #include "benchmark_multipoles4_d1.dat"
    },
    {
        #include "benchmark_multipoles0_d2.dat"
    },
    {
        #include "benchmark_multipoles2_d2.dat"
    },
    {
        #include "benchmark_multipoles4_d2.dat"
    },
    {
        #include "benchmark_multipoles0_g1.dat"
    },
    {
        #include "benchmark_multipoles2_g1.dat"
    },
    {
        #include "benchmark_multipoles4_g1.dat"
    },
    {
        #include "benchmark_multipoles0_g2.dat"
    },
    {
        #include "benchmark_multipoles2_g2.dat"
    },
    {
        #include "benchmark_multipoles4_g2.dat"
    },
    {
        #include "benchmark_multipoles0_g3.dat"
    },
    {
        #include "benchmark_multipoles2_g3.dat"
    },
    {
        #include "benchmark_multipoles4_g3.dat"
    },
    {
        #include "benchmark_multipoles0_g4.dat"
    },
    {
        #include "benchmark_multipoles2_g4.dat"
    },
    {
        #include "benchmark_multipoles4_g4.dat"
    },
    {
        #include "benchmark_multipoles0_g5.dat"
    },
    {
        #include "benchmark_multipoles2_g5.dat"
    },
    {
        #include "benchmark_multipoles4_g5.dat"
    },
    {
        #include "benchmark_multipoles0_len.dat"
    },
    {
        #include "benchmark_multipoles2_len.dat"
    },
    {
        #include "benchmark_multipoles4_len.dat"
    }

    /*
    {
        #include "benchmark_multipoles0_den_rsd.dat"
    },
    {
        #include "benchmark_multipoles2_den_rsd.dat"
    },
    {
        #include "benchmark_multipoles4_den_rsd.dat"
    },
    {
        #include "benchmark_multipoles0_den_rsd_len.dat"
    },
    {
        #include "benchmark_multipoles2_den_rsd_len.dat"
    },
    {
        #include "benchmark_multipoles4_den_rsd_len.dat"
    },
    {
        #include "benchmark_multipoles0_den.dat"
    }
    {
        #include "benchmark_multipoles0_den.dat"
    }
    {
        #include "benchmark_multipoles0_den.dat"
    }
    */
    };

    int counter = 0;

    for (int j = 0; j < sizeof(names[NAMES_MAXSIZE]) / sizeof(*(names[NAMES_MAXSIZE])); ++j){
        change_signal(&par->correlation_contrib, names[j], NULL);
        for (int i = 0; i < sizeof(multipoles) / sizeof(*multipoles); ++i){
            for (int k = 0; k < sizeof(xvalue) / sizeof(*xvalue); ++k){
                const double x = xvalue[k] * COFFE_H0;
                const double y_expected = yvalue[counter][k];
                const double y_obtained = coffe_integrate(
                            par, bg, integral,
                            x, 0, multipoles[i],
                            NONINTEGRATED, MULTIPOLES
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, 0, multipoles[i],
                            SINGLE_INTEGRATED, MULTIPOLES
                        )
                        +
                        coffe_integrate(
                            par, bg, integral,
                            x, 0, multipoles[i],
                            DOUBLE_INTEGRATED, MULTIPOLES
                        );

                fprintf(
                    stderr,
                    "l = %d, separation = %.3f, type %s\n",
                    multipoles[i], x, names[j]
                );
                weak_assert(
                    approx_equal(y_expected, y_obtained),
                    &error_flag
                );
            }
            ++counter;
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
