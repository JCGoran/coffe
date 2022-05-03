#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_errno.h>

#include "common.h"
#include "parser.h"
#include "background.h"
#include "integrals.h"
#include "functions.h"
#include "tools.h"
#include "covariance.h"


typedef struct legendre_parameters_t
{
    int n, m, a, b;
} legendre_parameters_t;


static double legendre_integral(
    double mu,
    void *p
)
{
    legendre_parameters_t *temp = (legendre_parameters_t *)p;

    return
        gsl_sf_legendre_Pl(temp->n, mu)
       *gsl_sf_legendre_Pl(temp->m, mu)
       *gsl_sf_legendre_Pl(temp->a, mu)
       *gsl_sf_legendre_Pl(temp->b, mu);
}


/**
    testing the integral of 4 Legendre polynomials
**/

static int coffe_test_legendre_integral(int ell_max)
{
    int error_flag = 0;

    for (int i = 0; i < ell_max; ++i){
    for (int j = 0; j < ell_max; ++j){
    for (int k = 0; k < ell_max; ++k){
    for (int l = 0; l < ell_max; ++l){

        const legendre_parameters_t test = {
            .n = i,
            .m = j,
            .a = k,
            .b = l
        };

        const double expected = coffe_integrate_1d(
            &legendre_integral,
            &test,
            -1, 1
        );

        const double obtained = coffe_legendre_integral(i, j, k, l);
        fprintf(
            stderr,
            "n = %d, m = %d, "
            "a = %d, b = %d, "
            "expected = %.3e, obtained = %.3e\n",
            i, j, k, l,
            expected, obtained
        );

        weak_assert(
            approx_equal_const_epsilon(
                expected, obtained
            ),
            &error_flag
        );
    }}}}

    if (!error_flag)
        COFFE_TESTS_PRINT_SUCCESS;

    return error_flag;
}



/**
    testing the various coefficients
**/
static int coffe_test_covariance_coefficients(
    const coffe_parameters_t *par,
    const coffe_background_t *bg
)
{
    const double z_mean_array[] = {0.5, 1.0, 1.5, 2.0};
    int error_flag = 0;

    for (size_t i = 0; i < COFFE_ARRAY_SIZE(z_mean_array); ++i){
        const double z_mean = z_mean_array[i];
        double galaxy_bias1 = 0;
        double galaxy_bias2 = 0;
        double growth_rate = 0;

        /* set bias to nonzero only if there is a density contribution */
        if (par->correlation_contrib.den){
            galaxy_bias1 = coffe_interp_spline(&par->galaxy_bias1, z_mean);
            galaxy_bias2 = coffe_interp_spline(&par->galaxy_bias2, z_mean);
        }

        /* set growth rate to nonzero only if there's an rsd contribution */
        if (par->correlation_contrib.rsd)
            growth_rate = coffe_interp_spline(&bg->f, z_mean);

        const double alpha0_11 = covariance_coefficient(
            galaxy_bias1,
            galaxy_bias1,
            growth_rate,
            0
        );

        const double alpha0_22 = covariance_coefficient(
            galaxy_bias2,
            galaxy_bias2,
            growth_rate,
            0
        );

        const double alpha0_cross = covariance_coefficient(
            galaxy_bias1,
            galaxy_bias2,
            growth_rate,
            0
        );

        const double alpha2_11 = covariance_coefficient(
            galaxy_bias1,
            galaxy_bias1,
            growth_rate,
            2
        );

        const double alpha2_22 = covariance_coefficient(
            galaxy_bias2,
            galaxy_bias2,
            growth_rate,
            2
        );

        const double alpha2_cross = covariance_coefficient(
            galaxy_bias1,
            galaxy_bias2,
            growth_rate,
            2
        );

        const double alpha4 = covariance_coefficient(
            galaxy_bias1,
            galaxy_bias2,
            growth_rate,
            4
        );

        /* b^2 + 2/3 b growth_rate + growth_rate^2/5 */
        const double c0 = pow(galaxy_bias1, 2) + 2 * galaxy_bias1 * growth_rate / 3. + pow(growth_rate, 2) / 5.;
        /* 4/3 b growth_rate + 4/7 growth_rate^2 */
        const double c2 = 4 * galaxy_bias1 * growth_rate / 3. + 4 * pow(growth_rate, 2) / 7.;
        /* 8/35 growth_rate^2 */
        const double c4 = 8 * pow(growth_rate, 2) / 35.;

        weak_assert(
            approx_equal_const_epsilon(
                c0, alpha0_11
            ),
            &error_flag
        );

        weak_assert(
            approx_equal_const_epsilon(
                c2, alpha2_11
            ),
            &error_flag
        );

        weak_assert(
            approx_equal_const_epsilon(
                c4, alpha4
            ),
            &error_flag
        );

    }

    return error_flag;
}



/**
    testing the covariance of multipoles of 2PCF
**/

static int coffe_test_covariance(
    const coffe_covariance_array_t *cov
)
{
    /* no errors initially */
    int error_flag = 0;
    /* disabling GSL's stupid error handler */
    gsl_error_handler_t *default_handler =
        gsl_set_error_handler_off();

    double *r1, *r2, *result;
    size_t size;

    /* TODO un-hardcode this */
    const int multipoles[] = {0, 2, 4};
    const size_t multipoles_size = COFFE_ARRAY_SIZE(multipoles);

    for (size_t mp1 = 0; mp1 < multipoles_size; ++mp1){
    for (size_t mp2 = 0; mp2 < multipoles_size; ++mp2){
        char name[256];
        snprintf(
            name,
            sizeof(name) / sizeof(*name),
            DATADIR "/tests/benchmarks/benchmark_multipoles_covariance_%d%d.dat",
            multipoles[mp1],
            multipoles[mp2]
        );
        coffe_read_ncol(
            name, 3, &size, &r1, &r2, &result
        );

        const size_t sep_size = (size_t)sqrt(size);

        for (size_t i = 0; i < sep_size; ++i){
        for (size_t j = 0; j < sep_size; ++j){
            fprintf(
                stderr,
                "l1 = %d, l2 = %d, "
                "r1 = %f, r2 = %f, "
                "expected = %.3e, obtained = %.3e\n",
                multipoles[mp1], multipoles[mp2],
                r1[i * sep_size + j], r2[i * sep_size + j],
                result[i * sep_size + j],
                coffe_covariance_find(
                    cov, 1.0, multipoles[mp1], multipoles[mp2],
                    r1[i * sep_size + j], r2[i * sep_size + j]
                ).value
            );
            weak_assert(
                approx_equal_const_epsilon(
                    result[i * sep_size + j],
                    coffe_covariance_find(
                        cov, 1.0, multipoles[mp1], multipoles[mp2],
                        r1[i * sep_size + j], r2[i * sep_size + j]
                    ).value
                ),
                &error_flag
            );
        }
        }

        /* memory cleanup */
        free(r1);
        free(r2);
        free(result);
    }
    }

    if (!error_flag)
        COFFE_TESTS_PRINT_SUCCESS;

    return error_flag;
}

int main(void)
{
    coffe_parameters_t par;
    coffe_parse_default_parameters(&par);

    /* need to set parameters manually */
    par.output_type = 4;
    /* in the new code, we can choose to compute just the covariance of a given tracer */
    par.correlation_contrib.rsd = 1;

    par.density1 = (double *)malloc(sizeof(double));
    par.density1[0] = 1e-3;
    par.density1_len = 1;

    par.density2 = (double *)malloc(sizeof(double));
    par.density2[0] = 1e-3;
    par.density2_len = 1;

    par.z_mean = (double *)malloc(sizeof(double));
    par.z_mean[0] = 1.0;
    par.z_mean_len = 1;

    par.deltaz = (double *)malloc(sizeof(double));
    par.deltaz[0] = 0.1;
    par.deltaz_len = 1;

    par.fsky = (double *)malloc(sizeof(double));
    par.fsky[0] = 0.2;
    par.fsky_len = 1;

    par.pixelsize = (double *)malloc(sizeof(double));
    par.pixelsize[0] = 50.0;
    par.pixelsize_len = 1;

    par.sep_len = 6;
    free(par.sep);
    par.sep = coffe_generate_range(50, 350, par.sep_len);

    coffe_background_t bg = {.flag = 0};
    coffe_background_init(&par, &bg);

    coffe_integral_array_t integral = {.array = NULL, .size = 0};
    coffe_integrals_init(&par, &bg, &integral);

    /* can't really integrate just one...*/
    coffe_covariance_array_t covariance = {.array = NULL, .size = 0};
    coffe_covariance_array_t dummy = {.array = NULL, .size = 0};
    coffe_covariance_init(&par, &bg, &covariance, &dummy);

    const int error_flag = coffe_test_covariance(&covariance);

    const int ell_max = 8;

    const int error_flag_legendre = coffe_test_legendre_integral(ell_max);

    const int error_flag_coefficients = coffe_test_covariance_coefficients(&par, &bg);

    coffe_parameters_free(&par);
    coffe_background_free(&bg);
    coffe_integrals_free(&integral);
    coffe_covariance_free(&covariance);

    return error_flag | error_flag_legendre | error_flag_coefficients;
}
