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

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_errno.h>

#ifdef HAVE_DOUBLE_EXPONENTIAL
#include "tanhsinh.h"
#endif

#ifdef HAVE_CUBA
#include "cuba.h"
#else
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "common.h"
#include "background.h"
#include "integrals.h"
#include "multipoles.h"
#include "functions.h"


struct multipoles_params
{
    struct coffe_background_t *bg;
    struct coffe_parameters_t *par;
    struct coffe_integrals_t *integral;
    double sep;
    int l;
};

static int multipoles_check_range(
    double **separations,
    size_t *len,
    double z_mean,
    double deltaz,
    struct coffe_background_t *bg
)
{
    qsort(*separations, *len, sizeof(double), coffe_compare_descending);

    double min_separation = (*separations)[*len - 1];
    double lower_limit = 0.1; /* arbitrary limit */

    size_t counter_neg = 0;
    for (size_t i = 0; i<*len; ++i){
        if ((*separations)[i] >= lower_limit){
            ++counter_neg;
        }
    }
    *separations = (double *)realloc(*separations, sizeof(double)*counter_neg);
    *len = counter_neg;

    /* checking for max separation */
    qsort(*separations, *len, sizeof(double), coffe_compare_ascending);

    double max_separation = (*separations)[*len - 1];
    double upper_limit =
        2*(
            interp_spline(&bg->comoving_distance, z_mean + deltaz)
           -interp_spline(&bg->comoving_distance, z_mean)
        )/COFFE_H0;

    size_t counter = 0;
    for (size_t i = 0; i<*len; ++i){
        if ((*separations)[i] <= upper_limit){
            ++counter;
        }
    }
    *separations = (double *)realloc(*separations, sizeof(double)*counter);
    *len = counter;
    if (min_separation < lower_limit){
        fprintf(
            stderr,
            "WARNING: minimum separation too low; "
            "cutting off list at the value %.2f "
            "Mpc/h\n", (*separations)[0]
        );
    }
    if (max_separation > upper_limit){
        fprintf(
            stderr,
            "WARNING: maximum separation too high; "
            "cutting off list at the value %.2f "
            "Mpc/h\n", (*separations)[*len - 1]
        );
    }
    return EXIT_SUCCESS;
}


/**
    calculates all the nonintegrated terms
**/

static double multipoles_nonintegrated_integrand(
    double x,
#ifdef HAVE_DOUBLE_EXPONENTIAL
    const void *p
#else
    void *p
#endif
)
{
    struct multipoles_params *all_params =
        (struct multipoles_params *) p;
    struct coffe_parameters_t *par = all_params->par;
    struct coffe_background_t *bg = all_params->bg;
    struct coffe_integrals_t *integral = all_params->integral;
    double sep = all_params->sep;
    int l = all_params->l;

    double mu = 2*x - 1;
    if (l == 0){
        return functions_nonintegrated(
            par, bg, integral,
            par->z_mean, mu, sep
        );
    }
    else{
        return functions_nonintegrated(
            par, bg, integral,
            par->z_mean, mu, sep
            )
           *gsl_sf_legendre_Pl(l, mu);
    }
}

static double multipoles_nonintegrated(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t *integral,
    double sep,
    int l
)
{
    struct multipoles_params test;
    test.par = par;
    test.bg = bg;
    test.integral = integral;
    test.l = l;
    test.sep = sep;

    double result, error;

#ifdef HAVE_DOUBLE_EXPONENTIAL
    result = tanhsinh_quad(
        &multipoles_nonintegrated_integrand,
        &test,
        0., 1., 0.,
        &error, NULL
    );
#else
    double prec = 1E-5;
    gsl_function integrand;
    integrand.function = &multipoles_nonintegrated_integrand;
    integrand.params = &test;

    gsl_integration_workspace *wspace =
        gsl_integration_workspace_alloc(COFFE_MAX_INTSPACE);
    gsl_integration_qag(
        &integrand, 0., 1., 0,
        prec, COFFE_MAX_INTSPACE,
        GSL_INTEG_GAUSS61, wspace,
        &result, &error
    );
    gsl_integration_workspace_free(wspace);
#endif
    return (2*l + 1)*result
        /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
}

#ifdef HAVE_CUBA
static int multipoles_single_integrated_integrand(
    const int *ndim, const cubareal var[],
    const int *ncomp, cubareal value[],
    void *p
)
#else
static double multipoles_single_integrated_integrand(
    double *var, size_t dim, void *p
)
#endif
{
    struct multipoles_params *params = (struct multipoles_params *) p;
    struct coffe_parameters_t *par = params->par;
    struct coffe_background_t *bg = params->bg;
    struct coffe_integrals_t *integral = params->integral;
    double sep = params->sep;

    double mu = 2*var[0] - 1, x = var[1];
#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_single_integrated(
            par, bg, integral, par->z_mean, mu, sep, x
        );
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_single_integrated(
            par, bg, integral, par->z_mean, mu, sep, x
        )
       *gsl_sf_legendre_Pl(params->l, mu);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_single_integrated(
            par, bg, integral, par->z_mean, mu, sep, x
        );
    }
    else{
        return functions_single_integrated(
            par, bg, integral, par->z_mean, mu, sep, x
        )
       *gsl_sf_legendre_Pl(params->l, mu);
    }
#endif
}

static double multipoles_single_integrated(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t integral[],
    double sep,
    double l
)
{
    const int dims = 2;

    struct multipoles_params test;
    test.par = par;
    test.bg = bg;
    test.integral = integral;
    test.sep = sep;
    test.l = l;

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
        multipoles_single_integrated_integrand,
        (void *)&test, 1,
        5e-4, 0, 0,
        1, par->integration_bins, 7,
        NULL, NULL,
        &nregions, &neval, &fail, result, error, prob
    );
    return (2*l + 1)*result[0]/interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
#else
    gsl_monte_function integrand;
    integrand.dim = dims;
    integrand.params = &test;
    integrand.f = &multipoles_single_integrated_integrand;

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
        default:
            result = 0;
            break;
    }
    gsl_rng_free(random);
    return (2*l + 1)*result
        /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
#endif
}


#ifdef HAVE_CUBA
static int multipoles_double_integrated_integrand(
    const int *ndim, const cubareal var[],
    const int *ncomp, cubareal value[],
    void *p
)
#else
static double multipoles_double_integrated_integrand(
    double *var, size_t dim, void *p
)
#endif
{
    struct multipoles_params *params = (struct multipoles_params *) p;
    struct coffe_parameters_t *par = params->par;
    struct coffe_background_t *bg = params->bg;
    struct coffe_integrals_t *integral = params->integral;
    double sep = params->sep;

    double mu = 2*var[0] - 1, x1 = var[1], x2 = var[2];
#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_double_integrated(
            par, bg, integral, par->z_mean, mu, sep, x1, x2
        );
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_double_integrated(
            par, bg, integral, par->z_mean, mu, sep, x1, x2
        )
       *gsl_sf_legendre_Pl(params->l, mu);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_double_integrated(
            par, bg, integral, par->z_mean, mu, sep, x1, x2
        );
    }
    else{
        return functions_double_integrated(
            par, bg, integral, par->z_mean, mu, sep, x1, x2
        )
       *gsl_sf_legendre_Pl(params->l, mu);
    }
#endif
}

static double multipoles_double_integrated(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t integral[],
    double sep,
    double l
)
{
    const int dims = 3;
    struct multipoles_params test;
    test.par = par;
    test.bg = bg;
    test.integral = integral;
    test.sep = sep;
    test.l = l;

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
        multipoles_double_integrated_integrand,
        (void *)&test, 1,
        5e-4, 0, 0,
        1, par->integration_bins, 7,
        NULL, NULL,
        &nregions, &neval, &fail, result, error, prob
    );
    return (2*l + 1)*result[0]/interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
#else
    gsl_monte_function integrand;
    integrand.dim = dims;
    integrand.params = &test;
    integrand.f = &multipoles_double_integrated_integrand;

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
        default:
            result = 0;
            break;
    }
    gsl_rng_free(random);
    return (2*l + 1)*result
        /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
#endif
}


int coffe_multipoles_init(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t *integral,
    struct coffe_multipoles_t *mp
)
{
#ifdef HAVE_CUBA
    cubacores(0, 10000);
#endif
    if (par->output_type == 2){
        mp->flag = 1;
        clock_t start, end;
        start = clock();

        if (par->verbose)
            printf("Calculating multipoles...\n");

        gsl_error_handler_t *default_handler =
            gsl_set_error_handler_off();

        alloc_double_matrix(
            &mp->result,
            par->multipole_values_len,
            par->sep_len
        );
        mp->l = (int *)coffe_malloc(sizeof(int)*par->multipole_values_len);
        for (int i = 0; i<par->multipole_values_len; ++i){
            mp->l[i] = (int)par->multipole_values[i];
        }
        mp->l_len = (size_t)par->multipole_values_len;

        mp->sep = (double *)coffe_malloc(sizeof(double)*par->sep_len);
        for (size_t i = 0; i<par->sep_len; ++i){
            mp->sep[i] = (double)par->sep[i];
        }
        mp->sep_len = (size_t)par->sep_len;
        multipoles_check_range(
            &mp->sep,
            &mp->sep_len,
            par->z_mean,
            par->deltaz,
            bg
        );

        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<mp->l_len; ++i){
            for (size_t j = 0; j<mp->sep_len; ++j){
                mp->result[i][j] =
                    multipoles_nonintegrated(
                        par, bg, integral,
                        mp->sep[j]*COFFE_H0, mp->l[i]
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<mp->l_len; ++i){
            for (size_t j = 0; j<mp->sep_len; ++j){
                mp->result[i][j] +=
                    multipoles_single_integrated(
                        par, bg, integral,
                        mp->sep[j]*COFFE_H0, mp->l[i]
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<mp->l_len; ++i){
            for (size_t j = 0; j<mp->sep_len; ++j){
                mp->result[i][j] +=
                    multipoles_double_integrated(
                        par, bg, integral,
                        mp->sep[j]*COFFE_H0, mp->l[i]
                    );
            }
        }

        end = clock();

        if (par->verbose)
            printf("Multipoles calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);

        gsl_set_error_handler(default_handler);
    }

    return EXIT_SUCCESS;
}

int coffe_multipoles_free(
    struct coffe_multipoles_t *mp
)
{
    if (mp->flag){
        for (size_t i = 0; i<mp->l_len; ++i){
            free(mp->result[i]);
        }
        free(mp->result);
        free(mp->l);
        free(mp->sep);
        mp->flag = 0;
    }
    return EXIT_SUCCESS;
}
