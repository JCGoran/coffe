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
    coffe_corrfunc_free(corrfunc);
    if (par->output_type == 1){

        clock_t start, end;
        start = clock();

        if (par->verbose)
            printf("Calculating the correlation function...\n");

        gsl_error_handler_t *default_handler =
            gsl_set_error_handler_off();

        corrfunc->size = par->corrfunc_coords.size;
        corrfunc->array = (coffe_corrfunc_t *)coffe_malloc(
            sizeof(coffe_corrfunc_t) * corrfunc->size
        );

        for (size_t i = 0; i < corrfunc->size; ++i){
            corrfunc->array[i].coords.z_mean = par->corrfunc_coords.array[i].z_mean;
            corrfunc->array[i].coords.separation = par->corrfunc_coords.array[i].separation;
            corrfunc->array[i].coords.mu = par->corrfunc_coords.array[i].mu;
        }

        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i < corrfunc->size; ++i){
            corrfunc->array[i].value = coffe_integrate(
                par, bg, integral,
                par->corrfunc_coords.array[i].z_mean,
                par->corrfunc_coords.array[i].separation,
                par->corrfunc_coords.array[i].mu,
                0,
                NONINTEGRATED, CORRFUNC
            );
        }

        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i < corrfunc->size; ++i){
            corrfunc->array[i].value += coffe_integrate(
                par, bg, integral,
                par->corrfunc_coords.array[i].z_mean,
                par->corrfunc_coords.array[i].separation,
                par->corrfunc_coords.array[i].mu,
                0,
                SINGLE_INTEGRATED, CORRFUNC
            );
        }

        #pragma omp parallel for num_threads(par->nthreads)
        for (size_t i = 0; i < corrfunc->size; ++i){
            corrfunc->array[i].value += coffe_integrate(
                par, bg, integral,
                par->corrfunc_coords.array[i].z_mean,
                par->corrfunc_coords.array[i].separation,
                par->corrfunc_coords.array[i].mu,
                0,
                DOUBLE_INTEGRATED, CORRFUNC
            );
        }

        end = clock();

        if (par->verbose)
            printf("Correlation function calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);

        gsl_set_error_handler(default_handler);

    }

    return EXIT_SUCCESS;
}

int coffe_corrfunc_free(
    coffe_corrfunc_array_t *cf
)
{
    if (cf->size)
        free(cf->array);
    cf->array = NULL;
    cf->size = 0;
    return EXIT_SUCCESS;
}

