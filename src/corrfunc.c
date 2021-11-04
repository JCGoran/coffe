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
#include "signal.h"

int coffe_corrfunc_coordinate_transform(
    const coffe_corrfunc_coordinate_array_t input,
    const enum coffe_corrfunc_coordinates_enum coord1,
    const enum coffe_corrfunc_coordinates_enum coord2,
    const enum coffe_corrfunc_coordinates_enum coord3,
    coffe_corrfunc_coordinate_array_t *output
)
{
    /* make sure they're all different */
    if (
        coord1 == coord2 ||
        coord2 == coord3 ||
        coord3 == coord1
    ){
        fprintf(
            stderr,
            "Duplicate coordinates encountered\n"
        );
        exit(EXIT_FAILURE);
    }

    /* placeholder: this is more complicated than I expected */
    return 0;
}

static coffe_corrfunc_t corrfunc_compute(
    coffe_parameters_t *par,
    coffe_background_t *bg,
    coffe_integral_array_t *integral,
    const double z_mean,
    const double separation,
    const double mu,
    const enum coffe_integral_type integral_type,
    const enum coffe_output_type output_type
)
{
    coffe_corrfunc_t result;

    result.coordinates.value[0].value = z_mean;
    result.coordinates.value[0].name = COFFE_COORDINATE_MEAN_REDSHIFT;

    result.coordinates.value[1].value = separation;
    result.coordinates.value[1].name = COFFE_COORDINATE_SEPARATION;

    result.coordinates.value[2].value = mu;
    result.coordinates.value[1].name = COFFE_COORDINATE_ANGLE_MU;

    result.value = coffe_integrate(
        par, bg, integral,
        z_mean,
        separation,
        mu,
        0,
        integral_type,
        output_type
    );

    return result;
}


const double r_parallel[] = {
0.1,0.2,0.4,0.8,1.,1.5,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.,50.,52.,54.,56.,58.,60.,62.,64.,66.,68.,70.,72.,74.,76.,78.,80.,81.,82.,83.,84.,85.,86.,87.,88.,89.,90.,91.,92.,93.,94.,95.,95.5,96.,96.5,97.,97.5,98.,98.5,99.,99.5,100.,100.5,101.,101.5,102.,102.5,103.,103.5,104.,104.5,105.,106.,107.,108.,109.,110.,112.,114.,116.,118.,120.,124.,128.,132.,136.,140.,144.,148.,152.,156.,160.,164.,168.,172.,176.,180.,185.,190.,195.,200.,205.,210.,215.,220.,225.,230.,235.,240.,250.,260.,270.,280.,290.,300.};

const double r_perpendicular[] = {
0.1,0.2,0.4,0.8,1.,1.5,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.,50.,52.,54.,56.,58.,60.,62.,64.,66.,68.,70.,72.,74.,76.,78.,80.,81.,82.,83.,84.,85.,86.,87.,88.,89.,90.,91.,92.,93.,94.,95.,95.5,96.,96.5,97.,97.5,98.,98.5,99.,99.5,100.,100.5,101.,101.5,102.,102.5,103.,103.5,104.,104.5,105.,106.,107.,108.,109.,110.,112.,114.,116.,118.,120.,124.,128.,132.,136.,140.,144.,148.,152.,156.,160.,164.,168.,172.,176.,180.,185.,190.,195.,200.,205.,210.,215.,220.,225.,230.,235.,240.,250.,260.,270.,280.,290.,300.};

const size_t r_p_len = sizeof(r_parallel)/sizeof(r_parallel[0]);

/**
    computes and stores the values of the correlation
    function
**/

int coffe_corrfunc_init(
    coffe_parameters_t *par,
    coffe_background_t *bg,
    coffe_integral_array_t *integral,
    coffe_corrfunc_array_t *corrfunc
)
{
#ifdef HAVE_CUBA
    {
        int n = 0, p = 10000;
        cubacores(&n, &p);
    }
#endif
    if (par->output_type == 1){

        clock_t start, end;
        start = clock();

        if (par->verbose)
            printf("Calculating the correlation function...\n");

        gsl_error_handler_t *default_handler =
            gsl_set_error_handler_off();

        corrfunc->size = par->corrfunc_output_coordinates.size;
        corrfunc->value =
            (coffe_corrfunc_t *)coffe_malloc(sizeof(coffe_corrfunc_t) * corrfunc->size);

        {
        size_t counter = 0;

        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i < corrfunc->size; ++i){
            if (
                coffe_check_range(
                    par->corrfunc_output_coordinates.value[i].separation,
                    par->corrfunc_output_coordinates.value[i].z_mean,
                    par->corrfunc_output_coordinates.value[i].deltaz,
                    bg
                )
            ){
                corrfunc->value[counter] = corrfunc_compute(
                    par, bg, integral,
                    par->corrfunc_output_coordinates.value[i].z_mean,
                    par->corrfunc_output_coordinates.value[i].separation,
                    par->corrfunc_output_coordinates.value[i].mu,
                    NONINTEGRATED, CORRFUNC
                );
            }
        }
        }

        {
        size_t counter = 0;
        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i < corrfunc->size; ++i){
            if (
                coffe_check_range(
                    par->corrfunc_output_coordinates.value[i].separation,
                    par->corrfunc_output_coordinates.value[i].z_mean,
                    par->corrfunc_output_coordinates.value[i].deltaz,
                    bg
                )
            ){
                corrfunc->value[counter].value += corrfunc_compute(
                    par, bg, integral,
                    par->corrfunc_output_coordinates.value[i].z_mean,
                    par->corrfunc_output_coordinates.value[i].separation,
                    par->corrfunc_output_coordinates.value[i].mu,
                    SINGLE_INTEGRATED, CORRFUNC
                ).value;
                ++counter;
            }
        }
        }

        {
        size_t counter = 0;
        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i < corrfunc->size; ++i){
            if (
                coffe_check_range(
                    par->corrfunc_output_coordinates.value[i].separation,
                    par->corrfunc_output_coordinates.value[i].z_mean,
                    par->corrfunc_output_coordinates.value[i].deltaz,
                    bg
                )
            ){
                corrfunc->value[counter].value += corrfunc_compute(
                    par, bg, integral,
                    par->corrfunc_output_coordinates.value[i].z_mean,
                    par->corrfunc_output_coordinates.value[i].separation,
                    par->corrfunc_output_coordinates.value[i].mu,
                    DOUBLE_INTEGRATED, CORRFUNC
                ).value;
                ++counter;
            }
        }
        }

        gsl_set_error_handler(default_handler);

        end = clock();

        if (par->verbose)
            printf("Correlation function calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);
    }

    return EXIT_SUCCESS;
}

int coffe_corrfunc_free(
    coffe_corrfunc_array_t *cf
)
{
    if (cf->size){
        free(cf->value);
        cf->size = 0;
    }
    return EXIT_SUCCESS;
}

