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
#include "signal.h"


static coffe_multipoles_t multipoles_compute(
    coffe_parameters_t *par,
    coffe_background_t *bg,
    coffe_integral_array_t *integral,
    const double z_mean,
    const double separation,
    const int l,
    const enum coffe_integral_type integral_type,
    const enum coffe_output_type output_type
)
{
    coffe_multipoles_t result;

    result.z_mean = z_mean;
    result.separation = separation;
    result.l = l;

    result.value = coffe_integrate(
        par, bg, integral,
        z_mean,
        separation,
        0,
        l,
        integral_type,
        output_type
    );

    return result;
}



int coffe_multipoles_init(
    coffe_parameters_t *par,
    coffe_background_t *bg,
    coffe_integral_array_t *integral,
    coffe_multipoles_array_t *mp
)
{
#ifdef HAVE_CUBA
    {
        int n = 0, p = 10000;
        cubacores(&n, &p);
    }
#endif
    if (par->output_type == 2){

        clock_t start, end;
        start = clock();

        if (par->verbose)
            printf("Calculating multipoles...\n");

        gsl_error_handler_t *default_handler =
            gsl_set_error_handler_off();

        mp->size = par->multipoles_output_coordinates.size;
        mp->value =
            (coffe_multipoles_t *)coffe_malloc(sizeof(coffe_multipoles_t) * mp->size);

        {
        size_t counter = 0;
        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i < mp->size; ++i){
            if (
                coffe_check_range(
                    par->multipoles_output_coordinates.value[i].separation,
                    par->multipoles_output_coordinates.value[i].z_mean,
                    par->multipoles_output_coordinates.value[i].deltaz,
                    bg
                )
            ){
                mp->value[counter] = multipoles_compute(
                    par, bg, integral,
                    par->multipoles_output_coordinates.value[i].z_mean,
                    par->multipoles_output_coordinates.value[i].separation,
                    par->multipoles_output_coordinates.value[i].l,
                    NONINTEGRATED, MULTIPOLES
                );
                ++counter;
            }
        }
        }
#if 0
        {
        size_t counter = 0;
        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i < mp->size; ++i){
            if (
                coffe_check_range(
                    par->multipoles_output_coordinates.value[i].separation,
                    par->multipoles_output_coordinates.value[i].z_mean,
                    par->multipoles_output_coordinates.value[i].deltaz,
                    bg
                )
            ){
                mp->value[counter].value += multipoles_compute(
                    par, bg, integral,
                    par->multipoles_output_coordinates.value[i].z_mean,
                    par->multipoles_output_coordinates.value[i].separation,
                    par->multipoles_output_coordinates.value[i].l,
                    SINGLE_INTEGRATED, MULTIPOLES
                ).value;
                ++counter;
            }
        }
        }

        {
        size_t counter = 0;
        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i < mp->size; ++i){
            if (
                coffe_check_range(
                    par->multipoles_output_coordinates.value[i].separation,
                    par->multipoles_output_coordinates.value[i].z_mean,
                    par->multipoles_output_coordinates.value[i].deltaz,
                    bg
                )
            ){
                mp->value[counter].value += multipoles_compute(
                    par, bg, integral,
                    par->multipoles_output_coordinates.value[i].z_mean,
                    par->multipoles_output_coordinates.value[i].separation,
                    par->multipoles_output_coordinates.value[i].l,
                    DOUBLE_INTEGRATED, MULTIPOLES
                ).value;
                ++counter;
            }
        }
        mp->size = counter;
        }
#endif

        mp->value = (coffe_multipoles_t *)
            realloc(mp->value, sizeof(coffe_multipoles_t) * mp->size);

        end = clock();

        if (par->verbose)
            printf("Multipoles calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);

        gsl_set_error_handler(default_handler);

    }

    return EXIT_SUCCESS;
}

int coffe_multipoles_free(
    coffe_multipoles_array_t *mp
)
{
    if (mp->size){
        free(mp->value);
        mp->size = 0;
    }
    return EXIT_SUCCESS;
}
