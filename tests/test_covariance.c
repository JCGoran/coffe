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
#include "tools.h"
#include "covariance.h"


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

    par.covariance_density = (double *)malloc(sizeof(double));
    par.covariance_density[0] = 1e-3;
    par.covariance_density_len = 1;

    par.covariance_z_mean = (double *)malloc(sizeof(double));
    par.covariance_z_mean[0] = 1.0;
    par.covariance_z_mean_len = 1;

    par.covariance_deltaz = (double *)malloc(sizeof(double));
    par.covariance_deltaz[0] = 0.1;
    par.covariance_deltaz_len = 1;

    par.covariance_fsky = (double *)malloc(sizeof(double));
    par.covariance_fsky[0] = 0.2;
    par.covariance_fsky_len = 1;

    par.covariance_pixelsize = (double *)malloc(sizeof(double));
    par.covariance_pixelsize[0] = 50.0;
    par.covariance_pixelsize_len = 1;

    par.covariance_step_size = 50.0;
    par.covariance_minimum_separation = 50.0;

    par.sep_len = 6;
    free(par.sep);
    par.sep = coffe_generate_range(50, 350, par.sep_len);

    par.covariance_coords.size =
          par.covariance_z_mean_len
        * par.sep_len
        * par.sep_len
        * par.multipole_values_len
        * par.multipole_values_len;

    par.covariance_coords.array = (coffe_covariance_coords_t *)coffe_malloc(
        sizeof(coffe_covariance_coords_t) * par.covariance_coords.size
    );

    coffe_parse_covariance_from_array(
        par.covariance_z_mean,
        par.covariance_z_mean_len,
        par.multipole_values,
        par.multipole_values_len,
        par.sep,
        par.sep_len,
        &par.covariance_coords
    );

    #ifdef _OPENMP
    par.nthreads = omp_get_num_procs();
    #endif

    coffe_background_t bg;
    coffe_background_init(&par, &bg);

    coffe_integral_array_t integral = {.array = NULL, .size = 0};
    coffe_integrals_init(&par, &bg, &integral);

    /* can't really integrate just one...*/
    coffe_covariance_array_t covariance = {.array = NULL, .size = 0};
    coffe_covariance_array_t dummy = {.array = NULL, .size = 0};
    coffe_covariance_init(&par, &bg, &covariance, &dummy);

    const int error_flag = coffe_test_covariance(&covariance);

    coffe_parameters_free(&par);
    coffe_background_free(&bg);
    coffe_integrals_free(&integral);
    coffe_covariance_free(&covariance);

    return error_flag;
}
