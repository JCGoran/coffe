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

#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline2d.h>
#include <string.h>
#include <libconfig.h>
#include <unistd.h>
#include <getopt.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "common.h"
#include "errors.h"
#include "parser.h"
#include "background.h"
#include "covariance.h"
#include "integrals.h"
#include "corrfunc.h"
#include "multipoles.h"
#include "average_multipoles.h"
#include "output.h"


int main(int argc, char *argv[])
{
    coffe_parameters_t par;
    coffe_background_t bg;
    coffe_integral_array_t integral = {.array = NULL, .size = 0};
    coffe_corrfunc_array_t cf = {.array = NULL, .size = 0};
    coffe_multipoles_array_t mp = {.array = NULL, .size = 0};
    coffe_average_multipoles_array_t ramp = {.array = NULL, .size = 0};
    coffe_covariance_array_t cov_mp = {.array = NULL, .size = 0};
    coffe_covariance_array_t cov_ramp = {.array = NULL, .size = 0};

    char settings_file[COFFE_MAX_STRLEN];

    int command;
    char *nthreads_opt = 0, *settings_opt = 0;
    while ((command = getopt(argc, argv, "hvCs:n:")) != -1){
        switch (command){
            case 'h':
                printf("Usage: coffe [FLAGS] [FILE]\n");
                printf("Options:\n");
                printf("\t-s [FILE]\t Process settings file [FILE]\n");
                printf("\t-n [NUMTHREADS]\t Use [NUMTHREADS] threads for the computation\n");
                printf("\t-h\t\t Show this help page and exit\n");
                printf("\t-C\t\t Show copyright information and exit\n");
                printf("\t-v\t\t Show version information and exit\n");
                return EXIT_SUCCESS;
                break;
            case 'v':
                printf("%s\n", PACKAGE_STRING);
#ifdef HAVE_CUBA
                printf("CUBA: enabled\n");
#endif
#ifdef HAVE_CLASS
                printf("CLASS: enabled\n");
#endif
#ifdef HAVE_DOUBLE_EXPONENTIAL
                printf("Double exponential quadrature: enabled\n");
#endif
                return EXIT_SUCCESS;
                break;
            case 'C':
                printf("%s", COFFE_COPYRIGHT);
                return EXIT_SUCCESS;
                break;
            case 's':
                settings_opt = optarg;
                break;
            case 'n':
                nthreads_opt = optarg;
                break;
        }
    }
    if (settings_opt == NULL){
        printf("Usage: coffe [FLAGS] [FILE]\n");
        printf("Use \"coffe -h\" to see command-line options\n");
        return EXIT_SUCCESS;
    }
    else{
        strncpy(settings_file, settings_opt, COFFE_MAX_STRLEN);
    }
    clock_t start, end;
    start = clock();

    int n;
    if (nthreads_opt != NULL){
        n = atoi(nthreads_opt);
    }
    else{
        #ifdef _OPENMP
        n = omp_get_max_threads();
        #else
        n = 1;
        #endif
    }

    if (n <= 0){
        print_error_verbose(PROG_VALUE_ERROR, "NUMTHREADS");
        exit(EXIT_FAILURE);
    }

    /* the main sequence */

    coffe_parser_init(settings_file, &par);

    /* this goes after the parser since it sets it to the default (1) */
    par.nthreads = n;

    if (par.verbose)
        printf("Number of threads in use: %d\n", par.nthreads);

    coffe_background_init(&par, &bg);

    coffe_integrals_init(&par, &bg, &integral);

    if (par.output_type == CORRFUNC)
        coffe_corrfunc_init(&par, &bg, &integral, &cf);

    if (par.output_type == MULTIPOLES)
        coffe_multipoles_init(&par, &bg, &integral, &mp);

    if (
        par.output_type == COVARIANCE_MULTIPOLES ||
        par.output_type == COVARIANCE_AVERAGE_MULTIPOLES
    )
        coffe_covariance_init(&par, &bg, &cov_mp, &cov_ramp);

    coffe_output_init(
        &par, &bg,
#ifdef HAVE_INTEGRALS
        &integral,
#endif
        &cf,
        &mp, &ramp,
        &cov_mp, &cov_ramp
    );

    /* freeing the memory */

    coffe_background_free(&bg);

    coffe_integrals_free(&integral);

    coffe_corrfunc_free(&cf);

    coffe_multipoles_free(&mp);

    coffe_covariance_free(&cov_mp);

    coffe_covariance_free(&cov_ramp);

    coffe_parameters_free(&par);

    end = clock();
    if (par.verbose)
        printf("Total program runtime is: %.2f s\n",
            (double)(end - start) / CLOCKS_PER_SEC);

    return EXIT_SUCCESS;
}
