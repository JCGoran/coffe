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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "common.h"
#include "background.h"
#include "covariance.h"
#include "integrals.h"
#include "corrfunc.h"
#include "multipoles.h"
#include "average_multipoles.h"
#include "output.h"
#include "errors.h"


/**
    makes the output directory if it doesn't exist;
    if it does, does nothing
**/

static int output_make_path(char *path)
{
    char temp_output[COFFE_MAX_STRLEN] = {0};
    struct stat temp_status = {0};
    char *current_dir;
    char temp_path[COFFE_MAX_STRLEN];
    strncpy(temp_path, path, COFFE_MAX_STRLEN);

    if (stat(path, &temp_status) == 0){
        if (path[strlen(path) - 1] != '/') strncat(path, "/", COFFE_MAX_STRLEN);
        fprintf(
                stderr,
                "WARNING: output directory exists, "
                "some files may be overwritten!\n"
        );
        return EXIT_SUCCESS;
    }
    else{
        if (path[0] == '/'){
            strncat(temp_output, "/", COFFE_MAX_STRLEN);
        }

        current_dir = strtok(temp_path, "/");
        while (current_dir != NULL){
            strncat(temp_output, current_dir, COFFE_MAX_STRLEN);
            if (stat(temp_output, &temp_status) == -1){
                if (mkdir(temp_output, 0700) != 0){
                    print_error(PROG_MKDIR_ERROR);
                    return EXIT_FAILURE;
                }
            }
            strncat(temp_output, "/", COFFE_MAX_STRLEN);
            current_dir = strtok(NULL, "/");
        }
    }
    strncpy(path, temp_output, COFFE_MAX_STRLEN);
    if (path[strlen(path) - 1] != '/') strncat(path, "/", COFFE_MAX_STRLEN);
    return EXIT_SUCCESS;
}

static int output_correlation_contributions_header(
    const struct coffe_correlation_contributions cc,
    char *header
)
{
    if (header != NULL){
        if (cc.den) strncat(header, "den ", COFFE_MAX_STRLEN);
        if (cc.rsd) strncat(header, "rsd ", COFFE_MAX_STRLEN);
        if (cc.d1) strncat(header, "d1 ", COFFE_MAX_STRLEN);
        if (cc.d2) strncat(header, "d2 ", COFFE_MAX_STRLEN);
        if (cc.g1) strncat(header, "g1 ", COFFE_MAX_STRLEN);
        if (cc.g2) strncat(header, "g2 ", COFFE_MAX_STRLEN);
        if (cc.g3) strncat(header, "g3 ", COFFE_MAX_STRLEN);
        if (cc.g4) strncat(header, "g4 ", COFFE_MAX_STRLEN);
        if (cc.g5) strncat(header, "g5 ", COFFE_MAX_STRLEN);
        if (cc.len) strncat(header, "len ", COFFE_MAX_STRLEN);
    }

    return EXIT_SUCCESS;
}


static int output_corrfunc(
    const char *filepath,
    const struct coffe_correlation_contributions cc,
    const coffe_corrfunc_array_t *cf
)
{
    char header[COFFE_MAX_STRLEN];
    FILE *data = fopen(filepath, "a");
    for (size_t i = 0; i < cf->size; ++i){
        snprintf(
            header, COFFE_MAX_STRLEN,
            "# correlation contributions: "
        );
        output_correlation_contributions_header(cc, header);
        strncat(header, "\n", COFFE_MAX_STRLEN);
        strncat(header, "# z_mean, sep, mu, value\n", COFFE_MAX_STRLEN);
        fprintf(data, header);
        fprintf(
            data,
            "%.10e %.10e %.10e %.10e\n",
            cf->array[i].coords.z_mean,
            cf->array[i].coords.separation / COFFE_H0,
            cf->array[i].coords.mu,
            cf->array[i].value
        );
    }
    fclose(data);
    return EXIT_SUCCESS;
}


static int output_multipoles(
    const char *filepath,
    const struct coffe_correlation_contributions cc,
    const coffe_multipoles_array_t *mp
)
{
    char header[COFFE_MAX_STRLEN];
    FILE *data = fopen(filepath, "a");
    snprintf(
        header, COFFE_MAX_STRLEN,
        "# correlation contributions: "
    );
    output_correlation_contributions_header(cc, header);
    strncat(header, "\n", COFFE_MAX_STRLEN);
    strncat(header, "# z_mean, l, sep, value\n", COFFE_MAX_STRLEN);
    fprintf(data, header);
    for (size_t i = 0; i < mp->size; ++i){
        fprintf(
            data,
            "%.10e %d %.10e %.10e\n",
            mp->array[i].coords.z_mean,
            mp->array[i].coords.l,
            mp->array[i].coords.separation / COFFE_H0,
            mp->array[i].value
        );
    }
    fclose(data);

    return EXIT_SUCCESS;
}


static int output_covariance(
    const char *filepath,
    const struct coffe_correlation_contributions cc,
    const coffe_covariance_array_t *cov
)
{
    char header[COFFE_MAX_STRLEN];
    FILE *output = fopen(filepath, "a");
    snprintf(
        header, COFFE_MAX_STRLEN,
        "# correlation contributions: "
    );
    output_correlation_contributions_header(cc, header);
    strncat(header, "\n", COFFE_MAX_STRLEN);
    strncat(header, "# z_mean, l1, l2, sep1, sep2, value\n", COFFE_MAX_STRLEN);
    fprintf(output, header);

    for (size_t i = 0; i < cov->size; ++i){
        fprintf(
            output,
            "%.10e %d %d %.10e %.10e %.10e\n",
            cov->array[i].coords.z_mean,
            cov->array[i].coords.l1,
            cov->array[i].coords.l2,
            cov->array[i].coords.separation1,
            cov->array[i].coords.separation2,
            cov->array[i].value
        );
    }
    fclose(output);

    return EXIT_SUCCESS;
}


/**
    specially to output the background functions
**/

static int output_background(
    char *filename,
    const char *sep,
    coffe_parameters_t *par,
    coffe_background_t *bg
)
{
    double *outputs[20] = {NULL}; // shouldn't be more than 20 anyway
    char header[COFFE_MAX_STRLEN] = {'\0'};
    char temp_header[COFFE_MAX_STRLEN];

    strncat(header, "# ", COFFE_MAX_STRLEN);
    for (size_t i = 0; i < par->type_bg_len; ++i){
        if (strcmp(par->type_bg[i], "z") == 0){
            outputs[i] = (double *)coffe_malloc(sizeof(double)*bg->a.spline->size);
            for (size_t n = 0; n<bg->a.spline->size; ++n){
                outputs[i][n] = bg->a.spline->x[n];
            }
            snprintf(temp_header, COFFE_MAX_STRLEN, "%zu)%s%s", i + 1, par->type_bg[i], sep);
            strncat(header, temp_header, COFFE_MAX_STRLEN);
        }
        if (strcmp(par->type_bg[i], "a") == 0){
            outputs[i] = (double *)coffe_malloc(sizeof(double)*bg->a.spline->size);
            for (size_t n = 0; n<bg->a.spline->size; ++n){
                outputs[i][n] = bg->a.spline->y[n];
            }
            snprintf(temp_header, COFFE_MAX_STRLEN, "%zu)%s%s", i + 1, par->type_bg[i], sep);
            strncat(header, temp_header, COFFE_MAX_STRLEN);
        }
        else if (strcmp(par->type_bg[i], "H") == 0){
            outputs[i] = (double *)coffe_malloc(sizeof(double)*bg->Hz.spline->size);
            for (size_t n = 0; n<bg->Hz.spline->size; ++n){
                outputs[i][n] = bg->Hz.spline->y[n]*COFFE_H0;
            }
            snprintf(temp_header, COFFE_MAX_STRLEN, "%zu)%s[h/Mpc]%s", i + 1, par->type_bg[i], sep);
            strncat(header, temp_header, COFFE_MAX_STRLEN);
        }
        else if (strcmp(par->type_bg[i], "conformal_H") == 0){
            outputs[i] = (double *)coffe_malloc(sizeof(double)*bg->conformal_Hz.spline->size);
            for (size_t n = 0; n<bg->conformal_Hz.spline->size; ++n){
                outputs[i][n] = bg->conformal_Hz.spline->y[n]*COFFE_H0;
            }
            snprintf(temp_header, COFFE_MAX_STRLEN, "%zu)%s[h/Mpc]%s", i + 1, par->type_bg[i], sep);
            strncat(header, temp_header, COFFE_MAX_STRLEN);
        }
        else if (strcmp(par->type_bg[i], "conformal_H_prime") == 0){
            outputs[i] = (double *)coffe_malloc(sizeof(double)*bg->conformal_Hz_prime.spline->size);
            for (size_t n = 0; n<bg->conformal_Hz_prime.spline->size; ++n){
                outputs[i][n] = bg->conformal_Hz_prime.spline->y[n]*COFFE_H0*COFFE_H0;
            }
            snprintf(temp_header, COFFE_MAX_STRLEN, "%zu)%s[h^2/Mpc^2]%s", i + 1, par->type_bg[i], sep);
            strncat(header, temp_header, COFFE_MAX_STRLEN);
        }
        else if (strcmp(par->type_bg[i], "D1") == 0){
            outputs[i] = (double *)coffe_malloc(sizeof(double)*bg->D1.spline->size);
            for (size_t n = 0; n<bg->D1.spline->size; ++n){
                outputs[i][n] = bg->D1.spline->y[n];
            }
            snprintf(temp_header, COFFE_MAX_STRLEN, "%zu)%s%s", i + 1, par->type_bg[i], sep);
            strncat(header, temp_header, COFFE_MAX_STRLEN);
        }
        else if (strcmp(par->type_bg[i], "f") == 0){
            outputs[i] = (double *)coffe_malloc(sizeof(double)*bg->f.spline->size);
            for (size_t n = 0; n<bg->f.spline->size; ++n){
                outputs[i][n] = bg->f.spline->y[n];
            }
            snprintf(temp_header, COFFE_MAX_STRLEN, "%zu)%s%s", i + 1, par->type_bg[i], sep);
            strncat(header, temp_header, COFFE_MAX_STRLEN);
        }
        else if (strcmp(par->type_bg[i], "comoving_distance") == 0){
            outputs[i] = (double *)coffe_malloc(sizeof(double)*bg->comoving_distance.spline->size);
            for (size_t n = 0; n<bg->comoving_distance.spline->size; ++n){
                outputs[i][n] = bg->comoving_distance.spline->y[n]/COFFE_H0;
            }
            snprintf(temp_header, COFFE_MAX_STRLEN, "%zu)%s[Mpc/h]%s", i + 1, par->type_bg[i], sep);
            strncat(header, temp_header, COFFE_MAX_STRLEN);
        }
    }
    strncat(header, "\n", COFFE_MAX_STRLEN);

    // I'm not proud of the way I handle this...
    write_ncol_null(
        filename,
        par->background_bins,
        header,
        sep,
        outputs[0],
        outputs[1],
        outputs[2],
        outputs[3],
        outputs[4],
        outputs[5],
        outputs[6],
        outputs[7],
        outputs[8],
        outputs[9],
        outputs[10],
        outputs[11],
        outputs[12],
        outputs[13],
        outputs[14],
        outputs[15],
        outputs[16],
        outputs[17],
        outputs[18],
        outputs[19]
    );
    for (size_t i = 0; i<20; ++i){
        free(outputs[i]);
    }
    return EXIT_SUCCESS;
}

int coffe_output_init(
    coffe_parameters_t *par,
    coffe_background_t *bg,
#ifdef HAVE_INTEGRALS
    coffe_integral_array_t *integral,
#endif
    coffe_corrfunc_array_t *cf,
    coffe_multipoles_array_t *mp,
    coffe_average_multipoles_array_t *ramp,
    coffe_covariance_array_t *cov_mp,
    coffe_covariance_array_t *cov_ramp
)
{
    clock_t start, end;
    start = clock();

    char filepath[COFFE_MAX_STRLEN];
    char prefix[COFFE_MAX_STRLEN];
    char header[COFFE_MAX_STRLEN];
    output_make_path(par->output_path);

    if (strcmp(par->output_prefix, "$TIME") == 0){
        snprintf(prefix, COFFE_MAX_STRLEN, "%s%s_", par->output_path, par->timestamp);
    }
    else{
        snprintf(prefix, COFFE_MAX_STRLEN, "%s%s", par->output_path, par->output_prefix);
    }

    if (par->verbose)
        printf("Writing output to \"%s\"...\n", prefix);

    /* settings file copy */
    snprintf(filepath, COFFE_MAX_STRLEN, "%ssettings.cfg", prefix);
    if (par->conf != NULL)
        config_write_file(par->conf, filepath);

    /* background */
    snprintf(filepath, COFFE_MAX_STRLEN, "%sbackground.dat", prefix);
    output_background(filepath, "\t", par, bg);

    /* correlation function (full) */
    if (par->output_type == 1){
        snprintf(filepath, COFFE_MAX_STRLEN, "%scorrfunc.dat", prefix);
        output_corrfunc(filepath, par->correlation_contrib, cf);
    }

    /* multipoles */
    else if (par->output_type == 2){
        snprintf(filepath, COFFE_MAX_STRLEN, "%smultipoles.dat", prefix);
        output_multipoles(filepath, par->correlation_contrib, mp);
    }

    /* covariance of multipoles */
    else if (par->output_type == 4){
        snprintf(filepath, COFFE_MAX_STRLEN, "%scovariance_multipoles.dat", prefix);
        output_covariance(filepath, par->correlation_contrib, cov_mp);
    }

    else if (par->output_type == 5){
        snprintf(
            filepath, COFFE_MAX_STRLEN,
            "%scovariance_average_multipoles.dat", prefix
        );

        FILE *output = fopen(filepath, "w");
        fprintf(output, "# sep[Mpc/h]\tsep[Mpc/h]\tresult\n");

        for (size_t i = 0; i < cov_ramp->size; ++i){
            fprintf(
                output,
                "%.10e %.10e %d %d %.10e %.10e\n",
                cov_ramp->array[i].coords.separation1,
                cov_ramp->array[i].coords.separation2,
                cov_ramp->array[i].coords.l1,
                cov_ramp->array[i].coords.l2,
                cov_ramp->array[i].coords.z_mean,
                cov_ramp->array[i].value
            );
        }
        fclose(output);
    }

#ifdef HAVE_INTEGRALS
    const struct nl_terms terms[] = {
        {.n = 0, .l = 0},
        {.n = 0, .l = 2},
        {.n = 0, .l = 4},
        {.n = 1, .l = 1},
        {.n = 1, .l = 3},
        {.n = 2, .l = 0},
        {.n = 2, .l = 2},
        {.n = 3, .l = 1}
    };

    for (size_t i = 0; i < sizeof(terms) / sizeof(*terms); ++i){
        snprintf(
            header, COFFE_MAX_STRLEN,
            "# n = %d, l = %d\n",
            terms[i].n,
            terms[i].l
        );
        snprintf(filepath, COFFE_MAX_STRLEN, "%sintegral%d.dat", prefix, i);
        write_ncol_null(
            filepath,
            coffe_find_integral(
                integral,
                terms[i].n,
                terms[i].l,
                COFFE_INTEGER,
                COFFE_INTEGER
            )->result.spline->size,
            header,
            " ",
            coffe_find_integral(
                integral,
                terms[i].n,
                terms[i].l,
                COFFE_INTEGER,
                COFFE_INTEGER
            )->result.spline->x,
            coffe_find_integral(
                integral,
                terms[i].n,
                terms[i].l,
                COFFE_INTEGER,
                COFFE_INTEGER
            )->result.spline->y,
            NULL
        );
    }

    if (par->divergent){

        /* the integral itself */
        snprintf(
            header, COFFE_MAX_STRLEN,
            "# n = 4, l = 0\n"
        );
        snprintf(filepath, COFFE_MAX_STRLEN, "%sintegral8.dat", prefix);
        write_ncol_null(
            filepath,
            coffe_find_integral(
                integral,
                4,
                0,
                COFFE_INTEGER,
                COFFE_INTEGER
            )->result.spline->size,
            header,
            " ",
            coffe_find_integral(
                integral,
                4,
                0,
                COFFE_INTEGER,
                COFFE_INTEGER
            )->result.spline->x,
            coffe_find_integral(
                integral,
                4,
                0,
                COFFE_INTEGER,
                COFFE_INTEGER
            )->result.spline->y,
            NULL
        );

        /* the renormalization */
        snprintf(
            header, COFFE_MAX_STRLEN,
            "# n = 4, l = 0\n"
        );
        snprintf(filepath, COFFE_MAX_STRLEN, "%sintegral8_renormalization.dat", prefix);
        FILE *file_renormalization = fopen(filepath, "w");
        for (size_t i = 0; i < coffe_find_integral(
            integral,
            1,
            -1,
            COFFE_HALF_INTEGER,
            COFFE_HALF_INTEGER
        )->renormalization.spline->interp_object.xsize; ++i){
            for (size_t j = 0; j < coffe_find_integral(
                integral,
                1,
                -1,
                COFFE_HALF_INTEGER,
                COFFE_HALF_INTEGER
            )->renormalization.spline->interp_object.ysize; ++j){
                fprintf(
                    file_renormalization,
                    "%e %e %e\n",
                    coffe_find_integral(
                        integral,
                        1,
                        -1,
                        COFFE_HALF_INTEGER,
                        COFFE_HALF_INTEGER
                    )->renormalization.spline->xarr[i],
                    coffe_find_integral(
                        integral,
                        1,
                        -1,
                        COFFE_HALF_INTEGER,
                        COFFE_HALF_INTEGER
                    )->renormalization.spline->yarr[j],
                    gsl_spline2d_get(
                        coffe_find_integral(
                            integral,
                            1,
                            -1,
                            COFFE_HALF_INTEGER,
                            COFFE_HALF_INTEGER
                        )->renormalization.spline,
                        coffe_find_integral(
                            integral,
                            1,
                            -1,
                            COFFE_HALF_INTEGER,
                            COFFE_HALF_INTEGER
                        )->renormalization.spline->zarr,
                        i, j
                    )
                );
            }
        }
        fclose(file_renormalization);
    }

    if (par->flatsky_local_nonlocal || par->flatsky_nonlocal){
        snprintf(
            header, COFFE_MAX_STRLEN,
            "# flatsky\n"
        );
        snprintf(filepath, COFFE_MAX_STRLEN, "%sintegral9.dat", prefix);
        write_ncol_null(
            filepath,
            coffe_find_integral(
                integral,
                1,
                -1,
                COFFE_HALF_INTEGER,
                COFFE_HALF_INTEGER
            )->result.spline->size,
            header,
            " ",
            coffe_find_integral(
                integral,
                1,
                -1,
                COFFE_HALF_INTEGER,
                COFFE_HALF_INTEGER
            )->result.spline->x,
            coffe_find_integral(
                integral,
                1,
                -1,
                COFFE_HALF_INTEGER,
                COFFE_HALF_INTEGER
            )->result.spline->y,
            NULL
        );
    }

#endif


    end = clock();

    if (par->verbose)
        printf("Output finished in %.2f s\n",
        (double)(end - start) / CLOCKS_PER_SEC);

    return EXIT_SUCCESS;
}
