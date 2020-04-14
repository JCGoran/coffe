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
    const struct coffe_background_t *bg
)
{
    const double z[] = {
        #include "BENCHMARK_Z.dat"
    };
    /* test H(z) */
    const double Hz[] = {
        #include "BENCHMARK_HZ.dat"
    };
    assert(sizeof(z) == sizeof(Hz));
    for (int i = 0; i < sizeof(z) / sizeof(*z); ++i)
        assert(
            approx_equal(
                Hz[i] / COFFE_H0,
                coffe_interp_spline(&bg->Hz, z[i])
            )
        );

    /* test conformal H(z) */
    const double conformal_Hz[] = {
        #include "BENCHMARK_CONFORMAL_HZ.dat"
    };
    assert(sizeof(z) == sizeof(conformal_Hz));
    for (int i = 0; i < sizeof(z) / sizeof(*z); ++i)
        assert(
            approx_equal(
                conformal_Hz[i] / COFFE_H0,
                coffe_interp_spline(&bg->conformal_Hz, z[i])
            )
        );

    /* test conformal_Hz_prime */
    const double conformal_Hz_prime[] = {
        #include "BENCHMARK_CONFORMAL_HZ_PRIME.dat"
    };
    assert(sizeof(z) == sizeof(conformal_Hz_prime));
    for (int i = 0; i < sizeof(z) / sizeof(*z); ++i)
        assert(
            approx_equal(
                conformal_Hz_prime[i] / COFFE_H0 / COFFE_H0,
                coffe_interp_spline(&bg->conformal_Hz_prime, z[i])
            )
        );

    /* test D1 */
    const double D1[] = {
        #include "BENCHMARK_D1.dat"
    };
    assert(sizeof(z) == sizeof(D1));
    for (int i = 0; i < sizeof(z) / sizeof(*z); ++i)
        assert(
            approx_equal(
                D1[i],
                coffe_interp_spline(&bg->D1, z[i])
            )
        );

    /* test f */
    const double f[] = {
        #include "BENCHMARK_F.dat"
    };
    assert(sizeof(z) == sizeof(f));
    for (int i = 0; i < sizeof(z) / sizeof(*z); ++i)
        assert(
            approx_equal(
                f[i],
                coffe_interp_spline(&bg->f, z[i])
            )
        );

    /* test comoving_distance */
    const double comoving_distance[] = {
        #include "BENCHMARK_COMOVING_DISTANCE.dat"
    };
    assert(sizeof(z) == sizeof(comoving_distance));
    for (int i = 0; i < sizeof(z) / sizeof(*z); ++i)
        assert(
            approx_equal(
                comoving_distance[i] * COFFE_H0,
                coffe_interp_spline(&bg->comoving_distance, z[i])
            )
        );

    COFFE_TESTS_PRINT_SUCCESS;

    return EXIT_SUCCESS;
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
    coffe_test_background(&bg);

    return 0;
}
