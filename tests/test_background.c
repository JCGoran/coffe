#include <math.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "common.h"
#include "parser.h"
#include "background.h"
#include "tools.h"

static int coffe_test_background(
    const coffe_background_t *bg
)
{
    int error_flag = 0;
    size_t size;
    double *z, *a, *Hz, *conformal_Hz, *conformal_Hz_prime, *D1, *f, *comoving_distance;

    /* load the file */
    coffe_read_ncol(
        DATADIR "/tests/benchmarks/benchmark_background.dat",
        8, &size,
        &z, &a, &Hz, &conformal_Hz,
        &conformal_Hz_prime, &D1, &f, &comoving_distance
    );

    for (size_t i = 0; i < size; ++i){

        /* test H(z) */
        weak_assert(
            approx_equal_const_epsilon(
                Hz[i] / COFFE_H0,
                coffe_interp_spline(&bg->Hz, z[i])
            ),
            &error_flag
        );

        /* test conformal H(z) */
        weak_assert(
            approx_equal_const_epsilon(
                conformal_Hz[i] / COFFE_H0,
                coffe_interp_spline(&bg->conformal_Hz, z[i])
            ),
            &error_flag
        );

        /* test conformal_Hz_prime */
        weak_assert(
            approx_equal_const_epsilon(
                conformal_Hz_prime[i] / COFFE_H0 / COFFE_H0,
                coffe_interp_spline(&bg->conformal_Hz_prime, z[i])
            ),
            &error_flag
        );

        /* test D1 */
        weak_assert(
            approx_equal_const_epsilon(
                D1[i] / D1[0],
                coffe_interp_spline(&bg->D1, z[i])
            ),
            &error_flag
        );

        /* test f */
        weak_assert(
            approx_equal_const_epsilon(
                f[i],
                coffe_interp_spline(&bg->f, z[i])
            ),
            &error_flag
        );

        /* test comoving_distance */
        weak_assert(
            approx_equal_const_epsilon(
                comoving_distance[i] * COFFE_H0,
                coffe_interp_spline(&bg->comoving_distance, z[i])
            ),
            &error_flag
        );
    }

    if (!error_flag)
        COFFE_TESTS_PRINT_SUCCESS;

    return error_flag;
}

int main(void)
{
    coffe_parameters_t par;
    coffe_parse_default_parameters(&par);

    par.divergent = 1;

    coffe_background_t bg = {.flag = 0};
    coffe_background_init(&par, &bg);

    const int error_flag = coffe_test_background(&bg);

    coffe_parameters_free(&par);
    coffe_background_free(&bg);

    return error_flag;
}
