#include <math.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "common.h"
#include "parser.h"
#include "background.h"
#include "integrals.h"
#include "tools.h"

static int coffe_test_integrals(
    const struct coffe_integrals_t *integrals
)
{
    int error_flag = 0;
    const size_t size_files = 8;
    size_t size;
    double *xvalue, *yvalue[size_files];

    /* load the files */
    for (size_t i = 0; i < size_files; ++i){
        const size_t size_name = 256;
        char name[size_name];
        snprintf(
            name,
            size_name,
            "tests/benchmarks/benchmark_integral%zu.dat",
            i
        );
        coffe_read_ncol(
            name,
            2, &size,
            &xvalue,
            &yvalue[i]
        );
    }

    /* compare the standard integrals (ones computes with 2FAST) */
    for (size_t integral = 0; integral < size_files; ++integral){
        for (size_t i = 0; i < size - 1; ++i){
            const double x = xvalue[i];
            const double y_expected = yvalue[integral][i];
            const double y_obtained = coffe_interp_spline(&integrals[integral].result, xvalue[i]);
            fprintf(
                stderr,
                "Integral %zu, separation %e, value %e\n",
                integral, x, y_obtained
            );
            weak_assert(
                approx_equal(y_expected, y_obtained),
                &error_flag
            );
        }
    }

    /* the divergent integral isn't computed with 2FAST */
    const double divergent_x[] = {
        #include "BENCHMARK_INTEGRALS8_X.dat"
    };
    const double divergent_y[] = {
        #include "BENCHMARK_INTEGRALS8_Y.dat"
    };
    assert(sizeof(divergent_x) == sizeof(divergent_y));

    for (size_t i = 0; i < sizeof(divergent_x) / sizeof(*divergent_x) - 1; ++i)
        weak_assert(
            approx_equal(
                divergent_y[i],
                coffe_interp_spline(
                    &integrals[8].result, divergent_x[i]
                )
            ),
            &error_flag
        );

    /* test renormalization */
    const double ren_x[] = {
        #include "BENCHMARK_INTEGRALS_RENORMALIZATION_X.dat"
    };
    const double ren_y[] = {
        #include "BENCHMARK_INTEGRALS_RENORMALIZATION_Y.dat"
    };
    const double ren_z[] = {
        #include "BENCHMARK_INTEGRALS_RENORMALIZATION_Z.dat"
    };

    assert(sizeof(ren_x) == sizeof(ren_y));
    assert(sizeof(ren_x) == sizeof(ren_z));

    for (size_t i = 0; i < sizeof(ren_x) / sizeof(*ren_x); ++i)
        weak_assert(
            approx_equal(
                ren_z[i] * pow(COFFE_H0, 4),
                gsl_spline2d_eval(
                    integrals[8].renormalization.spline,
                    ren_x[i] * COFFE_H0, ren_y[i] * COFFE_H0,
                    integrals[8].renormalization.xaccel,
                    integrals[8].renormalization.yaccel
                )
            ),
            &error_flag
        );

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

    const int error_flag = coffe_test_integrals(integrals);

    coffe_parameters_free(&par);
    coffe_background_free(&bg);
    coffe_integrals_free(integrals);

    return error_flag;
}
