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

#include "signal.h"

static double average_multipoles_compute(
    coffe_parameters_t *par,
    coffe_background_t *bg,
    coffe_integral_array_t *integral,
    const double z_min,
    const double z_max,
    const double separation,
    const int l,
    const enum coffe_integral_type integral_type,
    const enum coffe_output_type output_type
)
{
    return coffe_integrate(
        par, bg, integral,
        (z_min + z_max) / 2.,
        separation,
        0,
        l,
        integral_type,
        output_type
    );
}


int coffe_average_multipoles_init(
    coffe_parameters_t *par,
    coffe_background_t *bg,
    coffe_integral_array_t *integral,
    coffe_average_multipoles_array_t *ramp
)
{
    coffe_average_multipoles_free(ramp);
/*
    if (par->output_type == AVERAGE_MULTIPOLES){

        clock_t start, end;
        start = clock();

        if (par->verbose)
            printf("Calculating the redshift averaged multipoles...\n");

        gsl_error_handler_t *default_handler =
            gsl_set_error_handler_off();

        ramp->size = par->average_multipoles_coords.size;
        ramp->array = (coffe_average_multipoles_t *)coffe_malloc(
            sizeof(coffe_average_multipoles_t) * ramp->size
        );

        for (size_t i = 0; i < mp->size; ++i){
            ramp->array[i].coords.z_min = par->average_multipoles_coords.array[i].z_min;
            ramp->array[i].coords.z_max = par->average_multipoles_coords.array[i].z_max;
            ramp->array[i].coords.separation = par->average_multipoles_coords.array[i].separation;
            ramp->array[i].coords.l = par->multipoles_coords.array[i].l;
        }

        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i < ramp->size; ++i){
            ramp->array[i].value =
                average_multipoles_compute(
                    par, bg, integral,
                    ramp->array[i].coords.z_min,
                    ramp->array[i].coords.z_max,
                    ramp->array[i].coords.separation,
                    ramp->array[i].coords.l,
                    NONINTEGRATED, AVERAGE_MULTIPOLES
                );
        }

        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i < ramp->size; ++i){
            ramp->array[i].value +=
                average_multipoles_compute(
                    par, bg, integral,
                    ramp->array[i].coords.z_min,
                    ramp->array[i].coords.z_max,
                    ramp->array[i].coords.separation,
                    ramp->array[i].coords.l,
                    SINGLE_INTEGRATED, AVERAGE_MULTIPOLES
                );
        }

        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i < ramp->size; ++i){
            ramp->array[i].value +=
                average_multipoles_compute(
                    par, bg, integral,
                    ramp->array[i].coords.z_min,
                    ramp->array[i].coords.z_max,
                    ramp->array[i].coords.separation,
                    ramp->array[i].coords.l,
                    DOUBLE_INTEGRATED, AVERAGE_MULTIPOLES
                );
        }

        end = clock();

        if (par->verbose)
            printf("Redshift averaged multipoles calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);

        gsl_set_error_handler(default_handler);

    }

*/
    return EXIT_SUCCESS;
}


int coffe_average_multipoles_free(
    coffe_average_multipoles_array_t *ramp
)
{
    if (ramp->size)
        free(ramp->array);
    ramp->array = NULL;
    ramp->size = 0;
    return EXIT_SUCCESS;
}
