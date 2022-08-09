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


int coffe_multipoles_init(
    coffe_parameters_t *par,
    coffe_background_t *bg,
    coffe_integral_array_t *integral,
    coffe_multipoles_array_t *mp
)
{
    coffe_multipoles_free(mp);

    clock_t start, end;
    start = clock();

    if (par->verbose){
        printf("Calculating multipoles...\n");
        print_parameters(par);
    }

    gsl_error_handler_t *default_handler =
        gsl_set_error_handler_off();

    mp->size = par->multipole_values_len
        * par->z_mean_len
        * par->sep_len;

    mp->array = (coffe_multipoles_t *)coffe_malloc(
        sizeof(coffe_multipoles_t) * mp->size
    );

    {
    size_t counter = 0;
    for (size_t i = 0; i < par->z_mean_len; ++i){
    for (size_t j = 0; j < par->multipole_values_len; ++j){
    for (size_t k = 0; k < par->sep_len; ++k){
        mp->array[counter].coords.z_mean = par->z_mean[i];
        mp->array[counter].coords.l = par->multipole_values[j];
        mp->array[counter].coords.separation = par->sep[k];
        ++counter;
    }}}
    }

    #pragma omp parallel for num_threads(par->nthreads)
    for (size_t i = 0; i < mp->size; ++i){
        mp->array[i].value = coffe_integrate(
            par, bg, integral,
            mp->array[i].coords.z_mean,
            mp->array[i].coords.separation,
            0,
            mp->array[i].coords.l,
            NONINTEGRATED, MULTIPOLES
        );
    }

    #pragma omp parallel for num_threads(par->nthreads)
    for (size_t i = 0; i < mp->size; ++i){
        mp->array[i].value += coffe_integrate(
            par, bg, integral,
            mp->array[i].coords.z_mean,
            mp->array[i].coords.separation,
            0,
            mp->array[i].coords.l,
            SINGLE_INTEGRATED, MULTIPOLES
        );
    }

    #pragma omp parallel for num_threads(par->nthreads)
    for (size_t i = 0; i < mp->size; ++i){
        mp->array[i].value += coffe_integrate(
            par, bg, integral,
            mp->array[i].coords.z_mean,
            mp->array[i].coords.separation,
            0,
            mp->array[i].coords.l,
            DOUBLE_INTEGRATED, MULTIPOLES
        );
    }

    end = clock();

    if (par->verbose)
        printf("Multipoles calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);

    gsl_set_error_handler(default_handler);

    return EXIT_SUCCESS;
}

int coffe_multipoles_free(
    coffe_multipoles_array_t *mp
)
{
    if (mp->size)
        free(mp->array);
    mp->array = NULL;
    mp->size = 0;
    return EXIT_SUCCESS;
}
