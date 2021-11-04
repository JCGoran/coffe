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
    const struct coffe_covariance_t *covariance
)
{
    /* no errors initially */
    int error_flag = 0;
    /* disabling GSL's stupid error handler */
    gsl_error_handler_t *default_handler =
        gsl_set_error_handler_off();

    double *r1, *r2, *result;
    size_t size;

    for (size_t mp1 = 0; mp1 < covariance->l_len; ++mp1){
        for (size_t mp2 = 0; mp2 < covariance->l_len; ++mp2){
            char name[256];
            snprintf(
                name,
                sizeof(name) / sizeof(*name),
                DATADIR "/tests/benchmarks/benchmark_multipoles_covariance_%d%d.dat",
                covariance->l[mp1],
                covariance->l[mp2]
            );
            coffe_read_ncol(
                name, 3, &size, &r1, &r2, &result
            );

            /* before we do anything, make sure lengths correspond */
            assert(
                covariance->sep_len[0] * covariance->sep_len[0] == size
            );

            const size_t obtained_size = covariance->sep_len[0];

            for (size_t i = 0; i < obtained_size; ++i){
                for (size_t j = 0; j < obtained_size; ++j){
                    fprintf(
                        stderr,
                        "l1 = %d, l2 = %d\n"
                        "r1 = %.3e, r2 = %.3e, expected = %.3e, obtained = %.3e\n",
                        covariance->l[mp1], covariance->l[mp2],
                        r1[i * obtained_size + j], r2[i * obtained_size + j],
                        result[i * obtained_size + j],
                        covariance->result[0][covariance->l_len * mp1 + mp2][obtained_size * j + i]
                    );
                    weak_assert(
                        approx_equal_const_epsilon(
                            result[i * obtained_size + j],
                            covariance->result[0][covariance->l_len * mp1 + mp2][obtained_size * j + i]
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
    struct coffe_parameters_t par;
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

    #ifdef _OPENMP
    par.nthreads = omp_get_num_procs();
    #endif

    struct coffe_background_t bg;
    coffe_background_init(&par, &bg);

    struct coffe_integral_array_t integrals[10];
    coffe_integrals_init(&par, &bg, integrals);

    /* can't really integrate just one...*/
    struct coffe_covariance_t covariance, dummy;
    coffe_covariance_init(&par, &bg, &covariance, &dummy);

    const int error_flag = coffe_test_covariance(&covariance);

    coffe_parameters_free(&par);
    coffe_background_free(&bg);
    coffe_covariance_free(&covariance);

    return error_flag;
}
