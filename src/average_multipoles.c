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

#include "signal.h"

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

    double lower_limit = 0.0; /* arbitrary limit */

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
        z1 = coffe_interp_spline(
            &bg->z_as_chi,
            coffe_interp_spline(&bg->comoving_distance, zmin) + temp_sep/2.
        );
        z2 = coffe_interp_spline(
            &bg->z_as_chi,
            coffe_interp_spline(&bg->comoving_distance, zmax) - temp_sep/2.
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
        for (size_t i = 0; i<par->multipole_values_len; ++i){
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
                    coffe_integrate(
                        par, bg, integral,
                        ramp->sep[j]*COFFE_H0, 0, ramp->l[i],
                        NONINTEGRATED, AVERAGE_MULTIPOLES
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<ramp->l_len; ++i){
            for (size_t j = 0; j<ramp->sep_len; ++j){
                ramp->result[i][j] +=
                    coffe_integrate(
                        par, bg, integral,
                        ramp->sep[j]*COFFE_H0, 0, ramp->l[i],
                        SINGLE_INTEGRATED, AVERAGE_MULTIPOLES
                    );
            }
        }
        #pragma omp parallel for num_threads(par->nthreads) collapse(2)
        for (size_t i = 0; i<ramp->l_len; ++i){
            for (size_t j = 0; j<ramp->sep_len; ++j){
                ramp->result[i][j] +=
                    coffe_integrate(
                        par, bg, integral,
                        ramp->sep[j]*COFFE_H0, 0, ramp->l[i],
                        DOUBLE_INTEGRATED, AVERAGE_MULTIPOLES
                    );
            }
        }

        end = clock();

        if (par->verbose)
            printf("Redshift averaged multipoles calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);

        gsl_set_error_handler(default_handler);

        ramp->flag = 1;
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
