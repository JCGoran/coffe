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


/* computes the average multipole of non-integrated terms for given l and separation */

static double average_multipoles_nonintegrated(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t *integral,
    double sep,
    int l
)
{
    const int dims = 2;

    struct average_multipoles_params test;
    test.par = par;
    test.bg = bg;
    test.integral = integral;
    test.sep = sep;
    test.l = l;

    int flag = 0;
    int len = par->correlation_sources_len*(par->correlation_sources_len + 1)/2;
    for (int i = 0; i<len; ++i){
        if (
            strcmp(par->corr_terms[i], "00") == 0 ||
            strcmp(par->corr_terms[i], "11") == 0 ||
            strcmp(par->corr_terms[i], "22") == 0 ||
            strcmp(par->corr_terms[i], "33") == 0 ||
            strcmp(par->corr_terms[i], "44") == 0 ||
            strcmp(par->corr_terms[i], "55") == 0 ||
            strcmp(par->corr_terms[i], "66") == 0 ||
            strcmp(par->corr_terms[i], "01") == 0 ||
            strcmp(par->corr_terms[i], "10") == 0 ||
            strcmp(par->corr_terms[i], "02") == 0 ||
            strcmp(par->corr_terms[i], "20") == 0 ||
            strcmp(par->corr_terms[i], "03") == 0 ||
            strcmp(par->corr_terms[i], "30") == 0 ||
            strcmp(par->corr_terms[i], "04") == 0 ||
            strcmp(par->corr_terms[i], "40") == 0 ||
            strcmp(par->corr_terms[i], "05") == 0 ||
            strcmp(par->corr_terms[i], "50") == 0 ||
            strcmp(par->corr_terms[i], "06") == 0 ||
            strcmp(par->corr_terms[i], "60") == 0 ||
            strcmp(par->corr_terms[i], "12") == 0 ||
            strcmp(par->corr_terms[i], "21") == 0 ||
            strcmp(par->corr_terms[i], "13") == 0 ||
            strcmp(par->corr_terms[i], "31") == 0 ||
            strcmp(par->corr_terms[i], "14") == 0 ||
            strcmp(par->corr_terms[i], "41") == 0 ||
            strcmp(par->corr_terms[i], "15") == 0 ||
            strcmp(par->corr_terms[i], "51") == 0 ||
            strcmp(par->corr_terms[i], "16") == 0 ||
            strcmp(par->corr_terms[i], "61") == 0 ||
            strcmp(par->corr_terms[i], "23") == 0 ||
            strcmp(par->corr_terms[i], "32") == 0 ||
            strcmp(par->corr_terms[i], "24") == 0 ||
            strcmp(par->corr_terms[i], "42") == 0 ||
            strcmp(par->corr_terms[i], "25") == 0 ||
            strcmp(par->corr_terms[i], "52") == 0 ||
            strcmp(par->corr_terms[i], "26") == 0 ||
            strcmp(par->corr_terms[i], "62") == 0 ||
            strcmp(par->corr_terms[i], "34") == 0 ||
            strcmp(par->corr_terms[i], "43") == 0 ||
            strcmp(par->corr_terms[i], "35") == 0 ||
            strcmp(par->corr_terms[i], "53") == 0 ||
            strcmp(par->corr_terms[i], "36") == 0 ||
            strcmp(par->corr_terms[i], "63") == 0 ||
            strcmp(par->corr_terms[i], "45") == 0 ||
            strcmp(par->corr_terms[i], "54") == 0 ||
            strcmp(par->corr_terms[i], "46") == 0 ||
            strcmp(par->corr_terms[i], "64") == 0 ||
            strcmp(par->corr_terms[i], "56") == 0 ||
            strcmp(par->corr_terms[i], "65") == 0
        ) ++flag;
    }
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


/* computes the average multipole of single integrated terms for given l and separation */

static double average_multipoles_single_integrated(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t *integral,
    double sep,
    int l
)
{
    const int dims = 3;

    struct average_multipoles_params test;
    test.par = par;
    test.bg = bg;
    test.integral = integral;
    test.sep = sep;
    test.l = l;

    int flag = 0;
    int len = par->correlation_sources_len*(par->correlation_sources_len + 1)/2;
    for (int i = 0; i<len; ++i){
        if (
            strcmp(par->corr_terms[i], "17") == 0 ||
            strcmp(par->corr_terms[i], "71") == 0 ||
            strcmp(par->corr_terms[i], "18") == 0 ||
            strcmp(par->corr_terms[i], "81") == 0 ||
            strcmp(par->corr_terms[i], "19") == 0 ||
            strcmp(par->corr_terms[i], "91") == 0 ||
            strcmp(par->corr_terms[i], "27") == 0 ||
            strcmp(par->corr_terms[i], "72") == 0 ||
            strcmp(par->corr_terms[i], "28") == 0 ||
            strcmp(par->corr_terms[i], "82") == 0 ||
            strcmp(par->corr_terms[i], "29") == 0 ||
            strcmp(par->corr_terms[i], "92") == 0 ||
            strcmp(par->corr_terms[i], "37") == 0 ||
            strcmp(par->corr_terms[i], "73") == 0 ||
            strcmp(par->corr_terms[i], "38") == 0 ||
            strcmp(par->corr_terms[i], "83") == 0 ||
            strcmp(par->corr_terms[i], "39") == 0 ||
            strcmp(par->corr_terms[i], "93") == 0 ||
            strcmp(par->corr_terms[i], "47") == 0 ||
            strcmp(par->corr_terms[i], "74") == 0 ||
            strcmp(par->corr_terms[i], "48") == 0 ||
            strcmp(par->corr_terms[i], "84") == 0 ||
            strcmp(par->corr_terms[i], "49") == 0 ||
            strcmp(par->corr_terms[i], "94") == 0 ||
            strcmp(par->corr_terms[i], "57") == 0 ||
            strcmp(par->corr_terms[i], "75") == 0 ||
            strcmp(par->corr_terms[i], "58") == 0 ||
            strcmp(par->corr_terms[i], "85") == 0 ||
            strcmp(par->corr_terms[i], "59") == 0 ||
            strcmp(par->corr_terms[i], "95") == 0 ||
            strcmp(par->corr_terms[i], "67") == 0 ||
            strcmp(par->corr_terms[i], "76") == 0 ||
            strcmp(par->corr_terms[i], "68") == 0 ||
            strcmp(par->corr_terms[i], "86") == 0 ||
            strcmp(par->corr_terms[i], "69") == 0 ||
            strcmp(par->corr_terms[i], "96") == 0
        ) ++flag;
    }
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


/* computes the average multipole of double integrated terms for given l and separation */

static double average_multipoles_double_integrated(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t *integral,
    double sep,
    int l
)
{
    const int dims = 4;

    struct average_multipoles_params test;
    test.par = par;
    test.bg = bg;
    test.integral = integral;
    test.sep = sep;
    test.l = l;

    int flag = 0;
    int len = par->correlation_sources_len*(par->correlation_sources_len + 1)/2;
    for (int i = 0; i<len; ++i){
        if (
            strcmp(par->corr_terms[i], "77") == 0 || // g4-g4
            strcmp(par->corr_terms[i], "88") == 0 || // g5-g5
            strcmp(par->corr_terms[i], "99") == 0 || // len-len
            strcmp(par->corr_terms[i], "78") == 0 || // g4-g5
            strcmp(par->corr_terms[i], "87") == 0 || // g5-g4
            strcmp(par->corr_terms[i], "79") == 0 || // g4-len
            strcmp(par->corr_terms[i], "97") == 0 || // len-g4
            strcmp(par->corr_terms[i], "89") == 0 || // g5-len
            strcmp(par->corr_terms[i], "98") == 0    // len-g5
        ) ++flag;
    }
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
    if (par->output_type == 3){
        ramp->flag = 1;
        clock_t start, end;
        printf("Calculating the redshift averaged multipoles...\n");
        start = clock();

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
                    average_multipoles_nonintegrated(
                        par, bg, integral,
                        ramp->sep[j]*COFFE_H0, ramp->l[i]
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<ramp->l_len; ++i){
            for (size_t j = 0; j<ramp->sep_len; ++j){
                ramp->result[i][j] +=
                    average_multipoles_single_integrated(
                        par, bg, integral,
                        ramp->sep[j]*COFFE_H0, ramp->l[i]
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<ramp->l_len; ++i){
            for (size_t j = 0; j<ramp->sep_len; ++j){
                ramp->result[i][j] +=
                    average_multipoles_double_integrated(
                        par, bg, integral,
                        ramp->sep[j]*COFFE_H0, ramp->l[i]
                    );
            }
        }

        end = clock();

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
