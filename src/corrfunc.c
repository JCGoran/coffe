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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline2d.h>
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
#include "functions.h"
#include "corrfunc.h"

const double r_parallel[] = {
0.1,0.2,0.4,0.8,1.,1.5,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.,50.,52.,54.,56.,58.,60.,62.,64.,66.,68.,70.,72.,74.,76.,78.,80.,81.,82.,83.,84.,85.,86.,87.,88.,89.,90.,91.,92.,93.,94.,95.,95.5,96.,96.5,97.,97.5,98.,98.5,99.,99.5,100.,100.5,101.,101.5,102.,102.5,103.,103.5,104.,104.5,105.,106.,107.,108.,109.,110.,112.,114.,116.,118.,120.,124.,128.,132.,136.,140.,144.,148.,152.,156.,160.,164.,168.,172.,176.,180.,185.,190.,195.,200.,205.,210.,215.,220.,225.,230.,235.,240.,250.,260.,270.,280.,290.,300.};

const double r_perpendicular[] = {
0.1,0.2,0.4,0.8,1.,1.5,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.,50.,52.,54.,56.,58.,60.,62.,64.,66.,68.,70.,72.,74.,76.,78.,80.,81.,82.,83.,84.,85.,86.,87.,88.,89.,90.,91.,92.,93.,94.,95.,95.5,96.,96.5,97.,97.5,98.,98.5,99.,99.5,100.,100.5,101.,101.5,102.,102.5,103.,103.5,104.,104.5,105.,106.,107.,108.,109.,110.,112.,114.,116.,118.,120.,124.,128.,132.,136.,140.,144.,148.,152.,156.,160.,164.,168.,172.,176.,180.,185.,190.,195.,200.,205.,210.,215.,220.,225.,230.,235.,240.,250.,260.,270.,280.,290.,300.};

const size_t r_p_len = sizeof(r_parallel)/sizeof(r_parallel[0]);

struct corrfunc_params
{
    struct coffe_parameters_t *par;
    struct coffe_background_t *bg;
    struct coffe_integrals_t *integral;
    double mu;
    double sep;
};

static int corrfunc_check_range(
    double **separations,
    size_t *len,
    double z_mean,
    double deltaz,
    struct coffe_background_t *bg
)
{
    /* checking for min separation */
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


static double corrfunc_single_integrated_integrand(
    double x,
#ifdef HAVE_DOUBLE_EXPONENTIAL
    const void *p
#else
    void *p
#endif
)
{
    struct corrfunc_params *test =
        (struct corrfunc_params *) p;
    struct coffe_background_t *bg = test->bg;
    struct coffe_parameters_t *par = test->par;
    struct coffe_integrals_t *integral = test->integral;
    double mu = test->mu;
    double sep = test->sep;
    return
        functions_single_integrated(
            par, bg, integral,
            par->z_mean, mu, sep, x
        );
}


#ifdef HAVE_CUBA
static int corrfunc_double_integrated_integrand(
    const int *ndim, const cubareal var[],
    const int *ncomp, cubareal value[],
    void *p
)
#else
static double corrfunc_double_integrated_integrand(
    double *var, size_t dim, void *p
)
#endif
{
    struct corrfunc_params *test =
        (struct corrfunc_params *) p;
    struct coffe_background_t *bg = test->bg;
    struct coffe_parameters_t *par = test->par;
    struct coffe_integrals_t *integral = test->integral;
    double mu = test->mu;
    double sep = test->sep;
    double x1 = var[0], x2 = var[1];
#ifdef HAVE_CUBA
    value[0] =
#else
    return
#endif
        functions_double_integrated(
            par, bg, integral,
            par->z_mean, mu, sep, x1, x2
        );
#ifdef HAVE_CUBA
    return EXIT_SUCCESS;
#endif
}


static double corrfunc_compute(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t integral[],
    double mu,
    double sep,
    enum coffe_types corrfunc_flag
)
{
    switch(corrfunc_flag){
        case NONINTEGRATED:{
            return functions_nonintegrated(
                par, bg, integral,
                par->z_mean, mu, sep
            )/interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
        }
        case SINGLE_INTEGRATED:{
            const struct corrfunc_params test = {
            .par = par,
            .bg = bg,
            .integral = integral,
            .mu = mu,
            .sep = sep
            };
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

            double result, error;

            #ifdef HAVE_DOUBLE_EXPONENTIAL
            result = tanhsinh_quad(
                &corrfunc_single_integrated_integrand,
                &test,
                0., 1., 0.,
                &error, NULL
            );
            #else
            double prec = 1E-5;
            gsl_function integrand;
            integrand.function = &corrfunc_single_integrated_integrand;
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

            return result/interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
        }
        case DOUBLE_INTEGRATED:{
            const int dims = 2;

            struct corrfunc_params test;
            test.par = par;
            test.bg = bg;
            test.integral = integral;
            test.mu = mu;
            test.sep = sep;
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
                corrfunc_double_integrated_integrand,
                (void *)&test, 1,
                5e-4, 0, 0,
                1, par->integration_bins, 7,
                NULL, NULL,
                &nregions, &neval, &fail, result, error, prob
            );
            return result[0]/interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
            #else
            gsl_monte_function integrand;
            integrand.dim = dims;
            integrand.params = &test;
            integrand.f = &corrfunc_single_integrated_integrand;

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
            return result/interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
            #endif
        }
        default:
            return 0.0;
    }
}

/**
    computes and stores the values of the correlation
    function
**/

int coffe_corrfunc_init(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t *integral,
    struct coffe_corrfunc_ang_t *cf_ang,
    struct coffe_corrfunc_t *corrfunc,
    struct coffe_corrfunc2d_t *cf2d
)
{
#ifdef HAVE_CUBA
    cubacores(0, 10000);
#endif
    cf_ang->flag = 0;
    corrfunc->flag = 0;
    cf2d->flag = 0;
    if (par->output_type == 0){
        cf_ang->flag = 1;
        clock_t start, end;
        start = clock();

        if (par->verbose)
            printf("Calculating the angular correlation function...\n");

        double chi_mean = interp_spline(&bg->comoving_distance, par->z_mean);
        gsl_error_handler_t *default_handler =
            gsl_set_error_handler_off();

        const size_t theta_len = (size_t)par->theta_len;
        cf_ang->theta_len = (size_t)theta_len;
        cf_ang->theta = (double *)coffe_malloc(sizeof(double)*theta_len);
        const double maxangle = M_PI/2.;

        cf_ang->result = (double *)coffe_malloc(sizeof(double)*theta_len);

        for (size_t i = 0; i<theta_len; ++i){
            cf_ang->theta[i] = maxangle*(i + 1)/theta_len;
        }

        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i<theta_len; ++i){
            (cf_ang->result)[i] =
                corrfunc_compute(
                    par, bg, integral, 0, chi_mean*sqrt(2*(1. - cos(cf_ang->theta[i]))),
                    NONINTEGRATED
                );
        }
        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i<theta_len; ++i){
            (cf_ang->result)[i] +=
                corrfunc_compute(
                    par, bg, integral, 0, chi_mean*sqrt(2*(1. - cos(cf_ang->theta[i]))),
                    SINGLE_INTEGRATED
                );
        }
        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i<theta_len; ++i){
            (cf_ang->result)[i] +=
                corrfunc_compute(
                    par, bg, integral, 0, chi_mean*sqrt(2*(1. - cos(cf_ang->theta[i]))),
                    DOUBLE_INTEGRATED
                );
        }

        gsl_set_error_handler(default_handler);

        end = clock();

        if (par->verbose)
            printf("Angular correlation function calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);
    }
    else if (par->output_type == 1){
        corrfunc->flag = 1;
        clock_t start, end;
        start = clock();

        if (par->verbose)
            printf("Calculating the correlation function...\n");

        /* first index mu, second separations */
        alloc_double_matrix(
            &corrfunc->result,
            par->mu_len, par->sep_len
        );

        corrfunc->sep_len = (size_t)par->sep_len;
        corrfunc->sep =
            (double *)coffe_malloc(sizeof(double)*corrfunc->sep_len);
        for (size_t i = 0; i<corrfunc->sep_len; ++i){
            corrfunc->sep[i] = (double)par->sep[i];
        }
        corrfunc_check_range(
            &corrfunc->sep,
            &corrfunc->sep_len,
            par->z_mean,
            par->deltaz,
            bg
        );

        corrfunc->mu_len = (size_t)par->mu_len;
        corrfunc->mu =
            (double *)coffe_malloc(sizeof(double)*corrfunc->mu_len);
        for (size_t i = 0; i<corrfunc->mu_len; ++i){
            corrfunc->mu[i] = par->mu[i];
        }

        gsl_error_handler_t *default_handler =
            gsl_set_error_handler_off();

        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<corrfunc->mu_len; ++i){
            for (size_t j = 0; j<corrfunc->sep_len; ++j){
                (corrfunc->result)[i][j] =
                    corrfunc_compute(
                        par, bg, integral,
                        corrfunc->mu[i], corrfunc->sep[j]*COFFE_H0,
                        NONINTEGRATED
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<corrfunc->mu_len; ++i){
            for (size_t j = 0; j<corrfunc->sep_len; ++j){
                (corrfunc->result)[i][j] +=
                    corrfunc_compute(
                        par, bg, integral,
                        corrfunc->mu[i], corrfunc->sep[j]*COFFE_H0,
                        SINGLE_INTEGRATED
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<corrfunc->mu_len; ++i){
            for (size_t j = 0; j<corrfunc->sep_len; ++j){
                (corrfunc->result)[i][j] +=
                    corrfunc_compute(
                        par, bg, integral,
                        corrfunc->mu[i], corrfunc->sep[j]*COFFE_H0,
                        DOUBLE_INTEGRATED
                    );
            }
        }
        gsl_set_error_handler(default_handler);

        end = clock();

        if (par->verbose)
            printf("Correlation function calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);
    }

    else if (par->output_type == 6){
        cf2d->flag = 1;
        clock_t start, end;
        start = clock();

        if (par->verbose)
            printf("Calculating the 2D correlation function...\n");

        const double chi_mean = interp_spline(&bg->comoving_distance, par->z_mean);
        if (chi_mean < 320.*COFFE_H0){
            fprintf(
                stderr,
                "ERROR: z_mean too small for 2D correlation function!\n"
            );
            exit(EXIT_FAILURE);
        }

        /* first index mu, second separations */
        alloc_double_matrix(
            &cf2d->result, r_p_len, r_p_len
        );

        cf2d->sep_len = (size_t)r_p_len;
        cf2d->sep_parallel =
            (double *)coffe_malloc(sizeof(double)*cf2d->sep_len);
        cf2d->sep_perpendicular =
            (double *)coffe_malloc(sizeof(double)*cf2d->sep_len);
        for (size_t i = 0; i<cf2d->sep_len; ++i){
            cf2d->sep_parallel[i] = r_parallel[i];
            cf2d->sep_perpendicular[i] = r_perpendicular[i];
        }

        gsl_error_handler_t *default_handler =
            gsl_set_error_handler_off();

        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<cf2d->sep_len; ++i){
            for (size_t j = 0; j<cf2d->sep_len; ++j){
                (cf2d->result)[i][j] =
                    corrfunc_compute(
                        par, bg, integral,
                        cf2d->sep_parallel[i]/sqrt(pow(cf2d->sep_parallel[i], 2) + pow(cf2d->sep_perpendicular[j], 2)),
                        sqrt(pow(cf2d->sep_parallel[i], 2) + pow(cf2d->sep_perpendicular[j], 2))*COFFE_H0,
                        NONINTEGRATED
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<cf2d->sep_len; ++i){
            for (size_t j = 0; j<cf2d->sep_len; ++j){
                (cf2d->result)[i][j] +=
                    corrfunc_compute(
                        par, bg, integral,
                        cf2d->sep_parallel[i]/sqrt(pow(cf2d->sep_parallel[i], 2) + pow(cf2d->sep_perpendicular[j], 2)),
                        sqrt(pow(cf2d->sep_parallel[i], 2) + pow(cf2d->sep_perpendicular[j], 2))*COFFE_H0,
                        SINGLE_INTEGRATED
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<cf2d->sep_len; ++i){
            for (size_t j = 0; j<cf2d->sep_len; ++j){
                (cf2d->result)[i][j] +=
                    corrfunc_compute(
                        par, bg, integral,
                        cf2d->sep_parallel[i]/sqrt(pow(cf2d->sep_parallel[i], 2) + pow(cf2d->sep_perpendicular[j], 2)),
                        sqrt(pow(cf2d->sep_parallel[i], 2) + pow(cf2d->sep_perpendicular[j], 2))*COFFE_H0,
                        DOUBLE_INTEGRATED
                    );
            }
        }
        gsl_set_error_handler(default_handler);

        end = clock();

        if (par->verbose)
            printf("2D correlation function calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);
    }

    return EXIT_SUCCESS;
}

int coffe_corrfunc_ang_free(
   struct coffe_corrfunc_ang_t *cf_ang
)
{
    if (cf_ang->flag){
        free(cf_ang->theta);
        free(cf_ang->result);
        cf_ang->flag = 0;
    }
    return EXIT_SUCCESS;
}

int coffe_corrfunc_free(
    struct coffe_corrfunc_t *cf
)
{
    if (cf->flag){
        for (size_t i = 0; i<cf->mu_len; ++i){
            free(cf->result[i]);
        }
        free(cf->result);
        free(cf->sep);
        free(cf->mu);
        cf->flag = 0;
    }
    return EXIT_SUCCESS;
}

int coffe_corrfunc2d_free(
    struct coffe_corrfunc2d_t *cf2d
)
{
    if (cf2d->flag){
        for (size_t i = 0; i<cf2d->sep_len; ++i){
            free(cf2d->result[i]);
        }
        free(cf2d->result);
        free(cf2d->sep_parallel);
        free(cf2d->sep_perpendicular);
        cf2d->flag = 0;
    }
    return EXIT_SUCCESS;
}
