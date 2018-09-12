/*
 * This file is part of COFFE
 * Copyright (C) 2018 Goran Jelic-Cizmek
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
    struct coffe_parameters_t par;
    struct coffe_background_t bg;
#ifdef HAVE_NONLINEAR
    struct coffe_integrals_t integral[13];
#else
    struct coffe_integrals_t integral[9];
#endif
    struct coffe_corrfunc_ang_t cf_ang;
    struct coffe_corrfunc_t cf;
    struct coffe_multipoles_t mp;
    struct coffe_average_multipoles_t ramp;
    struct coffe_covariance_t cov_mp;
    struct coffe_covariance_t cov_ramp;
    struct coffe_corrfunc2d_t cf2d;

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
                printf("COFFE version %s\n", COFFE_VERSION_STRING);
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

    printf("   _____ ____  ______ ______ ______ \n");
    printf("  / ____/ __ \\|  ____|  ____|  ____|\n");
    printf(" | |   | |  | | |__  | |__  | |__   \n");
    printf(" | |   | |  | |  __| |  __| |  __|  \n");
    printf(" | |___| |__| | |    | |    | |____ \n");
    printf("  \\_____\\____/|_|    |_|    |______|\n");

    clock_t start, end;
    start = clock();

    int n;
    if (nthreads_opt != NULL){
        n = atoi(nthreads_opt);
    }
    else{
        n = 1;
    }

    if (n <= 0){
        print_error_verbose(PROG_VALUE_ERROR, "NUMTHREADS");
        exit(EXIT_FAILURE);
    }
    par.nthreads = n;
    printf("Number of threads in use: %d\n", par.nthreads);

    /* the main sequence */

    coffe_parser_init(settings_file, &par);

    coffe_background_init(&par, &bg);

    coffe_integrals_init(&par, &bg, integral);

    coffe_corrfunc_init(&par, &bg, integral, &cf_ang, &cf, &cf2d);

    coffe_multipoles_init(&par, &bg, integral, &mp);

    coffe_average_multipoles_init(&par, &bg, integral, &ramp);

    coffe_covariance_init(&par, &bg, &cov_mp, &cov_ramp);

    coffe_output_init(
        &par, &bg,
#ifdef HAVE_INTEGRALS
        integral,
#endif
        &cf_ang, &cf,
        &mp, &ramp,
        &cov_mp, &cov_ramp,
        &cf2d
    );

    /* freeing the memory */

    coffe_background_free(&bg);

    coffe_integrals_free(integral);

    coffe_corrfunc_ang_free(&cf_ang);

    coffe_corrfunc_free(&cf);

    coffe_corrfunc2d_free(&cf2d);

    coffe_multipoles_free(&mp);

    coffe_average_multipoles_free(&ramp);

    coffe_covariance_free(&cov_mp);

    coffe_covariance_free(&cov_ramp);

    end = clock();
    printf("Total program runtime is: %.2f s\n",
        (double)(end - start) / CLOCKS_PER_SEC);

    return EXIT_SUCCESS;
}
