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
#include "signal.h"


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
            coffe_interp_spline(&bg->comoving_distance, z_mean + deltaz)
           -coffe_interp_spline(&bg->comoving_distance, z_mean)
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
    mp->flag = 0;
    if (par->output_type == 2){

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

        for (size_t i = 0; i<mp->l_len; ++i){
            for (size_t j = 0; j<mp->sep_len; ++j){
                mp->result[i][j] = 0.0;
            }
        }

        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<mp->l_len; ++i){
            for (size_t j = 0; j<mp->sep_len; ++j){
                mp->result[i][j] +=
                    coffe_integrate(
                        par, bg, integral,
                        mp->sep[j]*COFFE_H0, 0, mp->l[i],
                        NONINTEGRATED, MULTIPOLES
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<mp->l_len; ++i){
            for (size_t j = 0; j<mp->sep_len; ++j){
                mp->result[i][j] +=
                    coffe_integrate(
                        par, bg, integral,
                        mp->sep[j]*COFFE_H0, 0, mp->l[i],
                        SINGLE_INTEGRATED, MULTIPOLES
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<mp->l_len; ++i){
            for (size_t j = 0; j<mp->sep_len; ++j){
                mp->result[i][j] +=
                    coffe_integrate(
                        par, bg, integral,
                        mp->sep[j]*COFFE_H0, 0, mp->l[i],
                        DOUBLE_INTEGRATED, MULTIPOLES
                    );
            }
        }

        end = clock();

        if (par->verbose)
            printf("Multipoles calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);

        gsl_set_error_handler(default_handler);

        mp->flag = 1;
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
