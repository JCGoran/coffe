/*
 * This file is part of COFFE
 * Copyright (C) 2019 Goran Jelic-Cizmek
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "common.h"
#include "background.h"
#include "integrals.h"
#include "functions.h"
#include "average_multipoles.h"

#ifdef HAVE_CUBA
#include "cuba.h"
#else
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#endif

struct average_multipoles_params
{
    struct coffe_background_t *bg;
    struct coffe_parameters_t *par;
    struct coffe_integrals_t *integral;
    double sep;
    int l;
};


static int average_multipoles_check_range(
    double **separations,
    size_t *len,
    double zmin,
    double zmax,
    struct coffe_background_t *bg
)
{
    /* checking for min separation */
    qsort(*separations, *len, sizeof(double), coffe_compare_descending);

    double lower_limit = 0.1; /* arbitrary limit */

    size_t counter_neg = 0;
    for (size_t i = 0; i<*len; ++i){
        if ((*separations)[i] >= lower_limit){
            ++counter_neg;
        }
    }
    *separations = (double *)realloc(*separations, sizeof(double)*counter_neg);
    *len = counter_neg;

    qsort(*separations, *len, sizeof(double), coffe_compare_ascending);

    size_t counter = 0;
    double z1, z2;
    for (size_t i = 0; i<*len; ++i){
        double temp_sep = (*separations)[i]*COFFE_H0;
        z1 = interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, zmin) + temp_sep/2.
        );
        z2 = interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, zmax) - temp_sep/2.
        );
        if (z1 < z2 && z1 > zmin && z2 < zmax && gsl_finite(z1) && gsl_finite(z2)){
            ++counter;
        }
    }
    *separations = (double *)realloc(*separations, sizeof(double)*counter);
    if (counter < *len){
        fprintf(
            stderr,
            "WARNING: separations too large for redshift averaged multipoles; "
            "cutting of list at the value %.2f "
            "Mpc/h\n", (*separations)[counter - 1]
        );
    }
    *len = counter;
    return EXIT_SUCCESS;
}



/* integrand of nonintegrated terms for redshift averaged multipoles */

#ifdef HAVE_CUBA
static int average_multipoles_nonintegrated_integrand(
    const int *ndim, const cubareal var[],
    const int *ncomp, cubareal value[],
    void *p
)
#else
static double average_multipoles_nonintegrated_integrand(
    double *var, size_t dim, void *p
)
#endif
{
    struct average_multipoles_params *params = (struct average_multipoles_params *) p;
    struct coffe_parameters_t *par = params->par;
    struct coffe_background_t *bg = params->bg;
    struct coffe_integrals_t *integral = params->integral;
    double sep = params->sep;

    double z1 =
        interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, par->z_min) + sep/2.
        );
    double z2 =
        interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, par->z_max) - sep/2.
        );

    double z = (z2 - z1)*var[0] + z1;
    double mu = 2*var[1] - 1;

#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_nonintegrated(
            par, bg, integral, z, mu, sep
        )
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_nonintegrated(
            par, bg, integral, z, mu, sep
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_nonintegrated(
            par, bg, integral, z, mu, sep
        )
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
    else{
        return functions_nonintegrated(
            par, bg, integral, z, mu, sep
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
#endif
}




/* integrand of single integrated terms for redshift averaged multipoles */

#ifdef HAVE_CUBA
static int average_multipoles_single_integrated_integrand(
    const int *ndim, const cubareal var[],
    const int *ncomp, cubareal value[],
    void *p
)
#else
static double average_multipoles_single_integrated_integrand(
    double *var, size_t dim, void *p
)
#endif
{
    struct average_multipoles_params *params = (struct average_multipoles_params *) p;
    struct coffe_parameters_t *par = params->par;
    struct coffe_background_t *bg = params->bg;
    struct coffe_integrals_t *integral = params->integral;
    double sep = params->sep;

    double z1 =
        interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, par->z_min) + sep/2.
        );
    double z2 =
        interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, par->z_max) - sep/2.
        );

    double z = (z2 - z1)*var[0] + z1;
    double mu = 2*var[1] - 1;
    double x = var[2];

#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_single_integrated(
            par, bg, integral, z, mu, sep, x
        )
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_single_integrated(
            par, bg, integral, z, mu, sep, x
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_single_integrated(
            par, bg, integral, z, mu, sep, x
        )
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
    else{
        return functions_single_integrated(
            par, bg, integral, z, mu, sep, x
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
#endif
}


/* integrand of double integrated terms for redshift averaged multipoles */

#ifdef HAVE_CUBA
static int average_multipoles_double_integrated_integrand(
    const int *ndim, const cubareal var[],
    const int *ncomp, cubareal value[],
    void *p
)
#else
static double average_multipoles_double_integrated_integrand(
    double *var, size_t dim, void *p
)
#endif
{
    struct average_multipoles_params *params = (struct average_multipoles_params *) p;
    struct coffe_parameters_t *par = params->par;
    struct coffe_background_t *bg = params->bg;
    struct coffe_integrals_t *integral = params->integral;
    double sep = params->sep;

    double z1 =
        interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, par->z_min) + sep/2.
        );
    double z2 =
        interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, par->z_max) - sep/2.
        );

    double z = (z2 - z1)*var[0] + z1;
    double mu = 2*var[1] - 1;
    double x1 = var[2], x2 = var[3];

#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_double_integrated(
            par, bg, integral, z, mu, sep, x1, x2
        )
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_double_integrated(
            par, bg, integral, z, mu, sep, x1, x2
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_double_integrated(
            par, bg, integral, z, mu, sep, x1, x2
        )
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
    else{
        return functions_double_integrated(
            par, bg, integral, z, mu, sep, x1, x2
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
#endif
}



/* computes the average multipole of all terms for given l and separation */

static double average_multipoles_compute(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t *integral,
    double sep,
    int l,
    enum coffe_types average_multipoles_flag
)
{
    const struct average_multipoles_params test = {
    .par = par,
    .bg = bg,
    .integral = integral,
    .sep = sep,
    .l = l
    };

    switch(average_multipoles_flag){
        case NONINTEGRATED:{
            const int dims = 2;
            int flag = 0;
            if (
                par->correlation_contrib.den ||
                par->correlation_contrib.rsd ||
                par->correlation_contrib.d1 ||
                par->correlation_contrib.d2 ||
                par->correlation_contrib.g1 ||
                par->correlation_contrib.g2 ||
                par->correlation_contrib.g3
            ) ++flag;
            if (flag == 0) return 0;

            #ifdef HAVE_CUBA
            int nregions, neval, fail;
            double result[1], error[1], prob[1];
            Cuhre(dims, 1,
                average_multipoles_nonintegrated_integrand,
                (void *)&test, 1,
                1e-3, 0, 0,
                1, par->integration_bins, 7,
                NULL, NULL,
                &nregions, &neval, &fail, result, error, prob
            );
            return (2*l + 1)*result[0]
            /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
            #else
            gsl_monte_function integrand;
            integrand.f = &average_multipoles_nonintegrated_integrand;
            integrand.dim = dims;
            integrand.params = &test;
            gsl_rng *random;
            gsl_rng_env_setup();
            const gsl_rng_type *T = gsl_rng_default;
            random = gsl_rng_alloc(T);
            double result, error;
            double lower[dims];
            double upper[dims];
            for (int i = 0; i<dims; ++i){
                lower[i] = 0.0;
                upper[i] = 1.0;
            }
            switch (par->integration_method){
                case 0:{
                    gsl_monte_plain_state *state =
                        gsl_monte_plain_alloc(dims);
                    gsl_monte_plain_integrate(
                        &integrand, lower, upper,
                        dims, par->integration_bins, random,
                        state,
                        &result, &error
                    );
                    gsl_monte_plain_free(state);
                    break;
                }
                case 1:{
                    gsl_monte_miser_state *state =
                        gsl_monte_miser_alloc(dims);
                    gsl_monte_miser_integrate(
                        &integrand, lower, upper,
                        dims, par->integration_bins, random,
                        state,
                        &result, &error
                    );
                    gsl_monte_miser_free(state);
                    break;
                }
                case 2:{
                    gsl_monte_vegas_state *state =
                        gsl_monte_vegas_alloc(dims);
                    gsl_monte_vegas_integrate(
                        &integrand, lower, upper,
                        dims, par->integration_bins, random,
                        state,
                        &result, &error
                    );
                    gsl_monte_vegas_free(state);
                    break;
                }
            }

            gsl_rng_free(random);

            return (2*l + 1)*result
            /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
            #endif
        }
        case SINGLE_INTEGRATED:{
            const int dims = 3;

            int flag = 0;
            if (
                (par->correlation_contrib.len  && par->correlation_contrib.den) ||
                (par->correlation_contrib.len  && par->correlation_contrib.rsd) ||
                (par->correlation_contrib.len  && par->correlation_contrib.d1) ||
                (par->correlation_contrib.len  && par->correlation_contrib.d2) ||
                (par->correlation_contrib.len  && par->correlation_contrib.g1) ||
                (par->correlation_contrib.len  && par->correlation_contrib.g2) ||
                (par->correlation_contrib.len  && par->correlation_contrib.g3) ||
                (par->correlation_contrib.g4  && par->correlation_contrib.den) ||
                (par->correlation_contrib.g4  && par->correlation_contrib.rsd) ||
                (par->correlation_contrib.g4  && par->correlation_contrib.d1) ||
                (par->correlation_contrib.g4  && par->correlation_contrib.d2) ||
                (par->correlation_contrib.g4  && par->correlation_contrib.g1) ||
                (par->correlation_contrib.g4  && par->correlation_contrib.g2) ||
                (par->correlation_contrib.g4  && par->correlation_contrib.g3) ||
                (par->correlation_contrib.g5  && par->correlation_contrib.den) ||
                (par->correlation_contrib.g5  && par->correlation_contrib.rsd) ||
                (par->correlation_contrib.g5  && par->correlation_contrib.d1) ||
                (par->correlation_contrib.g5  && par->correlation_contrib.d2) ||
                (par->correlation_contrib.g5  && par->correlation_contrib.g1) ||
                (par->correlation_contrib.g5  && par->correlation_contrib.g2) ||
                (par->correlation_contrib.g5  && par->correlation_contrib.g3)
            ) ++flag;
            if (flag == 0) return 0;

            #ifdef HAVE_CUBA
            int nregions, neval, fail;
            double result[1], error[1], prob[1];
            Cuhre(dims, 1,
                average_multipoles_single_integrated_integrand,
                (void *)&test, 1,
                5e-4, 0, 0,
                1, par->integration_bins, 7,
                NULL, NULL,
                &nregions, &neval, &fail, result, error, prob
            );
            return (2*l + 1)*result[0]
            /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
            #else
            gsl_monte_function integrand;
            integrand.f = &average_multipoles_single_integrated_integrand;
            integrand.dim = dims;
            integrand.params = &test;
            gsl_rng *random;
            gsl_rng_env_setup();
            const gsl_rng_type *T = gsl_rng_default;
            random = gsl_rng_alloc(T);
            double result, error;
            double lower[dims];
            double upper[dims];
            for (int i = 0; i<dims; ++i){
                lower[i] = 0.0;
                upper[i] = 1.0;
            }

            switch (par->integration_method){
                case 0:{
                    gsl_monte_plain_state *state =
                        gsl_monte_plain_alloc(dims);
                    gsl_monte_plain_integrate(
                        &integrand, lower, upper,
                        dims, par->integration_bins, random,
                        state,
                        &result, &error
                    );
                    gsl_monte_plain_free(state);
                    break;
                }
                case 1:{
                    gsl_monte_miser_state *state =
                        gsl_monte_miser_alloc(dims);
                    gsl_monte_miser_integrate(
                        &integrand, lower, upper,
                        dims, par->integration_bins, random,
                        state,
                        &result, &error
                    );
                    gsl_monte_miser_free(state);
                    break;
                }
                case 2:{
                    gsl_monte_vegas_state *state =
                        gsl_monte_vegas_alloc(dims);
                    gsl_monte_vegas_integrate(
                        &integrand, lower, upper,
                        dims, par->integration_bins, random,
                        state,
                        &result, &error
                    );
                    gsl_monte_vegas_free(state);
                    break;
                }
            }

            gsl_rng_free(random);

            return (2*l + 1)*result
            /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
            #endif
        }
        case DOUBLE_INTEGRATED:{
            const int dims = 4;

            int flag = 0;
            if (
                par->correlation_contrib.len ||
                par->correlation_contrib.g4 ||
                par->correlation_contrib.g5
            ) ++flag;
            if (flag == 0) return 0;

            #ifdef HAVE_CUBA
            int nregions, neval, fail;
            double result[1], error[1], prob[1];
            Cuhre(dims, 1,
                average_multipoles_double_integrated_integrand,
                (void *)&test, 1,
                5e-4, 0, 0,
                1, par->integration_bins, 7,
                NULL, NULL,
                &nregions, &neval, &fail, result, error, prob
            );
            return (2*l + 1)*result[0]
            /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
            #else
            gsl_monte_function integrand;
            integrand.f = &average_multipoles_double_integrated_integrand;
            integrand.dim = dims;
            integrand.params = &test;
            gsl_rng *random;
            gsl_rng_env_setup();
            const gsl_rng_type *T = gsl_rng_default;
            random = gsl_rng_alloc(T);
            double result, error;
            double lower[dims];
            double upper[dims];
            for (int i = 0; i<dims; ++i){
                lower[i] = 0.0;
                upper[i] = 1.0;
            }

            switch (par->integration_method){
                case 0:{
                    gsl_monte_plain_state *state =
                        gsl_monte_plain_alloc(dims);
                    gsl_monte_plain_integrate(
                        &integrand, lower, upper,
                        dims, par->integration_bins, random,
                        state,
                        &result, &error
                    );
                    gsl_monte_plain_free(state);
                    break;
                }
                case 1:{
                    gsl_monte_miser_state *state =
                        gsl_monte_miser_alloc(dims);
                    gsl_monte_miser_integrate(
                        &integrand, lower, upper,
                        dims, par->integration_bins, random,
                        state,
                        &result, &error
                    );
                    gsl_monte_miser_free(state);
                    break;
                }
                case 2:{
                    gsl_monte_vegas_state *state =
                        gsl_monte_vegas_alloc(dims);
                    gsl_monte_vegas_integrate(
                        &integrand, lower, upper,
                        dims, par->integration_bins, random,
                        state,
                        &result, &error
                    );
                    gsl_monte_vegas_free(state);
                    break;
                }
            }

            gsl_rng_free(random);

            return (2*l + 1)*result
            /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
            #endif
        }
        default:
            return 0.0;
    }
}


int coffe_average_multipoles_init(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t *integral,
    struct coffe_average_multipoles_t *ramp
)
{
#ifdef HAVE_CUBA
    cubacores(0, 10000);
#endif
    ramp->flag = 0;
    if (par->output_type == 3){
        ramp->flag = 1;
        clock_t start, end;
        start = clock();

        if (par->verbose)
            printf("Calculating the redshift averaged multipoles...\n");

        gsl_error_handler_t *default_handler =
            gsl_set_error_handler_off();

        alloc_double_matrix(
            &ramp->result,
            par->multipole_values_len,
            par->sep_len
        );
        ramp->l = (int *)coffe_malloc(sizeof(int)*par->multipole_values_len);
        for (int i = 0; i<par->multipole_values_len; ++i){
            ramp->l[i] = (int)par->multipole_values[i];
        }
        ramp->l_len = (size_t)par->multipole_values_len;

        ramp->sep = (double *)coffe_malloc(sizeof(double)*par->sep_len);
        for (size_t i = 0; i<par->sep_len; ++i){
            ramp->sep[i] = (double)par->sep[i];
        }
        ramp->sep_len = (size_t)par->sep_len;
        average_multipoles_check_range(
            &ramp->sep, &ramp->sep_len,
            par->z_min, par->z_max, bg
        );

        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<ramp->l_len; ++i){
            for (size_t j = 0; j<ramp->sep_len; ++j){
                ramp->result[i][j] =
                    average_multipoles_compute(
                        par, bg, integral,
                        ramp->sep[j]*COFFE_H0, ramp->l[i],
                        NONINTEGRATED
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<ramp->l_len; ++i){
            for (size_t j = 0; j<ramp->sep_len; ++j){
                ramp->result[i][j] +=
                    average_multipoles_compute(
                        par, bg, integral,
                        ramp->sep[j]*COFFE_H0, ramp->l[i],
                        SINGLE_INTEGRATED
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<ramp->l_len; ++i){
            for (size_t j = 0; j<ramp->sep_len; ++j){
                ramp->result[i][j] +=
                    average_multipoles_compute(
                        par, bg, integral,
                        ramp->sep[j]*COFFE_H0, ramp->l[i],
                        DOUBLE_INTEGRATED
                    );
            }
        }

        end = clock();

        if (par->verbose)
            printf("Redshift averaged multipoles calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);

        gsl_set_error_handler(default_handler);
    }

    return EXIT_SUCCESS;
}

int coffe_average_multipoles_free(
    struct coffe_average_multipoles_t *ramp
)
{
    if (ramp->flag){
        for (size_t i = 0; i<ramp->l_len; ++i){
            free(ramp->result[i]);
        }
        free(ramp->result);
        free(ramp->l);
        free(ramp->sep);
        ramp->flag = 0;
    }
    return EXIT_SUCCESS;
}
