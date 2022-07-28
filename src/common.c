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
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <math.h>
#include <gsl/gsl_version.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coupling.h>

#ifdef HAVE_CLASS
#include "class.h"
#endif

#ifdef HAVE_CUBA
#include "cuba.h"
#else
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#endif


#include "common.h"
#include "errors.h"

double coffe_max_array_double(
    const double *array,
    const size_t size
){
    if (size > 0){
        double current = array[0];
        for (size_t i = 1; i < size; ++i){
            if (current < array[i]){
                current = array[i];
            }
        }
        return current;
    }
    return NAN;
}

/**
    another malloc function, checks if malloc succeded
**/

void *coffe_malloc(size_t len)
{
    if (!len){
        print_error(PROG_ALLOC_ERROR);
    }
    void *values = malloc(len);
    if (values == NULL){
        print_error(PROG_ALLOC_ERROR);
    }
    return values;
}


/**
    gets the current time
**/

char *coffe_get_time(void)
{
    char *timestamp = (char *)malloc(sizeof(char) * COFFE_MAX_STRLEN);
    time_t ltime;
    ltime = time(NULL);
    struct tm *tm;
    tm = localtime(&ltime);

    snprintf(
        timestamp,
        COFFE_MAX_STRLEN,
        "%04d-%02d-%02d-%02d-%02d-%02d",
        /* see http://www.cplusplus.com/reference/ctime/tm/ for the offsets */
        tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday, tm->tm_hour, tm->tm_min, tm->tm_sec
    );
    return timestamp;
}


/**
    reads file <filename> into array <values>
    of length <len>
**/

int read_1col(
    const char *filename,
    double **values,
    size_t *len
)
{
    int error = 0;
    FILE *data = fopen(filename, "r");
    if (data == NULL){
        print_error_verbose(PROG_OPEN_ERROR, filename);
        return EXIT_FAILURE;
    }
    size_t n = 0;
    char c, temp_string[COFFE_MAX_STRLEN];

    while ((c=fgetc(data)) != EOF){
        if (c == '\n') ++n;
    }

    error = fseek(data, 0, SEEK_SET);
    if (error){
        print_error_verbose(PROG_POS_ERROR, filename);
        return EXIT_FAILURE;
    }

    size_t counter = 0;
    for (size_t i = 0; i<n; ++i){
        fgets(temp_string, COFFE_MAX_STRLEN, data);
        if (temp_string[0] != '#') ++counter;
    }
    *len = counter;

    error = fseek(data, 0, SEEK_SET);
    if (error){
        print_error_verbose(PROG_POS_ERROR, filename);
        return EXIT_FAILURE;
    }

    *values = (double *)coffe_malloc(sizeof(double)*counter);
    counter = 0;

    for (size_t i = 0; i<n; ++i){
        fgets(temp_string, COFFE_MAX_STRLEN, data);
        if (temp_string[0] != '#'){
            *(*values + counter) = atof(temp_string);
            ++counter;
        }
    }

    error = fclose(data);
    if (error){
        print_error_verbose(PROG_CLOSE_ERROR, filename);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


/**
    reads file <filename> into arrays <values1>
    and <values2> of length <len>
**/

int read_2col(
    const char *filename,
    double **values1,
    double **values2,
    size_t *len
)
{
    int error = 0;
    FILE *data = fopen(filename, "r");
    if (data == NULL){
        print_error_verbose(PROG_OPEN_ERROR, filename);
        return EXIT_FAILURE;
    }
    size_t n = 0;
    char c, temp_string[COFFE_MAX_STRLEN];
    char *temp_token;

    while ((c=fgetc(data)) != EOF){
        if (c == '\n') ++n;
    }

    error = fseek(data, 0, SEEK_SET);
    if (error){
        print_error_verbose(PROG_POS_ERROR, filename);
        return EXIT_FAILURE;
    }

    size_t counter = 0;
    for (size_t i = 0; i<n; ++i){
        fgets(temp_string, COFFE_MAX_STRLEN, data);
        if (temp_string[0] != '#') ++counter;
    }
    *len = counter;

    error = fseek(data, 0, SEEK_SET);
    if (error){
        print_error_verbose(PROG_POS_ERROR, filename);
        return EXIT_FAILURE;
    }

    *values1 = (double *)coffe_malloc(sizeof(double)*counter);
    *values2 = (double *)coffe_malloc(sizeof(double)*counter);
    counter = 0;

    for (size_t i = 0; i<n; ++i){
        fgets(temp_string, COFFE_MAX_STRLEN, data);
        if (temp_string[0] != '#'){
            temp_token = strtok(temp_string, ",\t: ");
            *(*values1 + counter) = atof(temp_token);
            temp_token = strtok(NULL, ",\t: ");
            *(*values2 + counter) = atof(temp_token);
            ++counter;
        }
    }

    error = fclose(data);
    if (error){
        print_error_verbose(PROG_CLOSE_ERROR, filename);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


/**
    writes array <values> of length <len>
    into file <filename> with header <header>
**/

int write_1col(
    const char *filename,
    double *values,
    size_t len,
    const char *header
)
{
    int error = 0;
    FILE *data = fopen(filename, "w");
    if (data == NULL){
        print_error_verbose(PROG_OPEN_ERROR, filename);
        return EXIT_FAILURE;
    }

    if (header != NULL)
        fprintf(data, "%s", header);

    for (size_t i = 0; i<len; ++i)
        fprintf(data, "%.15e\n", values[i]);

    error = fclose(data);
    if (error){
        print_error_verbose(PROG_CLOSE_ERROR, filename);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


/**
    writes arrays <values1> and <values2> of length <len>
    into file <filename> with header <header>
    using separator(s) <sep>
**/

int write_2col(
    const char *filename,
    const double *values1,
    const double *values2,
    size_t len,
    const char *header,
    const char *sep
)
{
    int error = 0;
    FILE *data = fopen(filename, "w");
    if (data == NULL){
        print_error_verbose(PROG_OPEN_ERROR, filename);
        return EXIT_FAILURE;
    }

    if (header != NULL)
        fprintf(data, "%s", header);

    for (size_t i = 0; i<len; ++i)
        fprintf(data, "%.15e%s%.15e\n", values1[i], sep, values2[i]);

    error = fclose(data);
    if (error){
        print_error_verbose(PROG_CLOSE_ERROR, filename);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


/**
    writes an arbitrary number of columns into file <filename>,
    each column having length <len>,
    with header <header> and using separator <sep>
    NOTE: C can't figure out by itself the number of
    arguments, so the first parameter is the number of columns
**/

int write_ncol(
    size_t ncolumns,
    const char *filename,
    size_t len,
    const char *header,
    const char *sep,
    const double *values,
    ...
)
{
    va_list args;
    FILE *data = fopen(filename, "w");
    if (data == NULL){
        print_error_verbose(PROG_OPEN_ERROR, filename);
        return EXIT_FAILURE;
    }
    if (header != NULL) fprintf(data, "%s", header);

    va_start(args, values);
    double **all_values = (double **)coffe_malloc(sizeof(double *)*ncolumns);
    all_values[0] = (double *)values;

    for (size_t i = 1; i<ncolumns; ++i){
        values = va_arg(args, double *);
        all_values[i] = (double *)values;
    }
    va_end(args);

    for (size_t i = 0; i<len; ++i){
        if (all_values[i] != NULL){
            for (size_t j = 0; j<ncolumns; ++j){
                if (j<ncolumns - 1)
                    fprintf(data, "%.10e%s", all_values[j][i], sep);
                else
                    fprintf(data, "%.10e%s\n", all_values[j][i], sep);
            }
        }
    }
    fclose(data);
    for (size_t i = 0; i<ncolumns; ++i){
        all_values[i] = NULL;
    }
    free(all_values);
    return EXIT_SUCCESS;
}


/**
    modification of the above, but without the counter;
    last argument MUST be NULL for the function to work properly!
**/

int write_ncol_null(
    const char *filename,
    size_t len,
    const char *header,
    const char *sep,
    const double *values,
    ...
)
{
    va_list args;
    va_start(args, values);
    size_t counter = 0;
    double **all_values = (double **)coffe_malloc(sizeof(double *)*COFFE_MAX_ALLOCABLE);
    all_values[0] = (double *)values;

    while (values != NULL){
        values = va_arg(args, double *);
        ++counter;
        all_values[counter] = (double *)values;
    }
    va_end(args);

    if (counter > 0){
        FILE *data = fopen(filename, "w");
        if (data == NULL){
            print_error_verbose(PROG_OPEN_ERROR, filename);
            return EXIT_FAILURE;
        }
        if (header != NULL) fprintf(data, "%s", header);

        for (size_t i = 0; i<len; ++i){
            for (size_t j = 0; j<counter; ++j){
                if (j<counter - 1)
                    fprintf(data, "%.10e%s", all_values[j][i], sep);
                else
                    fprintf(data, "%.10e%s\n", all_values[j][i], sep);
            }
        }
        fclose(data);
    }
    for (size_t i = 0; i<counter; ++i){
        all_values[i] = NULL;
    }
    free(all_values);
    return EXIT_SUCCESS;
}



/**
    writes a <len1>x<len2> matrix <values> into file <filename>
**/

int write_matrix(
    const char *filename,
    const double **values,
    size_t len1,
    size_t len2,
    const char *header,
    const char *sep
)
{
    int error = 0;
    FILE *data = fopen(filename, "w");
    if (data == NULL){
        print_error(PROG_OPEN_ERROR);
        return EXIT_FAILURE;
    }

    if (header != NULL){
        fprintf(data, "%s\n", header);
    }

    for (size_t i = 0; i<len1; ++i){
        for (size_t j = 0; j<len2; ++j){
            if (j != len2 - 1)
                fprintf(data, "%10e%s", values[i][j], sep);
            else
                fprintf(data, "\n");
        }
    }

    error = fclose(data);
    if (error){
        print_error(PROG_CLOSE_ERROR);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


/**
    allocates an <len1>x<len2> matrix
    and stores it into <values>
**/

int alloc_double_matrix(
    double ***values,
    size_t len1,
    size_t len2
)
{
    *values = (double **)coffe_malloc(sizeof(double *)*len1);
    for (size_t i = 0; i<len1; ++i){
        (*values)[i] = (double *)coffe_malloc(sizeof(double)*len2);
    }
    return EXIT_SUCCESS;
}


int free_double_matrix(
    double ***values,
    size_t len
)
{
    *values = (double **)coffe_malloc(sizeof(double *)*len);
    for (size_t i = 0; i<len; ++i){
        free((*values)[i]);
    }
    free(*values);
    return EXIT_SUCCESS;
}


int copy_matrix_array(
    double **destination,
    const double **source,
    size_t rows,
    size_t columns,
    size_t index,
    const char *type
)
{
    if (strcmp(type, "row") == 0){
        *destination = (double *)coffe_malloc(sizeof(double)*columns);
        if (index >= rows){
            print_error(PROG_FAIL);
            exit(EXIT_FAILURE);
        }
        for (size_t i = 0; i<columns; ++i)
            (*destination)[i] = source[index][i];
    }
    else if (strcmp(type, "column") == 0){
        *destination = (double *)coffe_malloc(sizeof(double)*rows);
        if (index >= columns){
            print_error(PROG_FAIL);
            exit(EXIT_FAILURE);
        }
        for (size_t i = 0; i<rows; ++i)
            (*destination)[i] = source[i][index];
    }

    return EXIT_SUCCESS;
}


int coffe_init_spline(
    coffe_interpolation *interp,
    const double *xi,
    const double *yi,
    const size_t bins,
    const enum coffe_interp1d_type interpolation_type
)
{
    if (bins <= 0){
        print_error(PROG_VALUE_ERROR);
        exit(EXIT_FAILURE);
    }
    const gsl_interp_type *T;
    switch(interpolation_type){
        case COFFE_INTERP_LINEAR:
            T = gsl_interp_linear;
            break;
        case COFFE_INTERP_POLYNOMIAL:
            T = gsl_interp_polynomial;
            break;
        case COFFE_INTERP_CSPLINE:
            T = gsl_interp_cspline;
            break;
        case COFFE_INTERP_CSPLINE_PERIODIC:
            T = gsl_interp_cspline_periodic;
            break;
        case COFFE_INTERP_AKIMA:
            T = gsl_interp_akima;
            break;
        case COFFE_INTERP_AKIMA_PERIODIC:
            T = gsl_interp_akima_periodic;
            break;
        case COFFE_INTERP_STEFFEN:{
#if GSL_MAJOR_VERSION > 1
                T = gsl_interp_steffen;
                break;
#else
                fprintf(stderr,
                    "Interpolation type \"Steffen\" requires GSL 2.1 or above!\n");
                exit(EXIT_SUCCESS);
#endif
        }
        default:
            T = gsl_interp_akima;
            break;
    }
    interp->spline
        = gsl_spline_alloc(T, bins);
    interp->accel
        = gsl_interp_accel_alloc();
    gsl_spline_init(interp->spline, xi, yi, bins);
    return EXIT_SUCCESS;
}


int coffe_init_spline2d(
    coffe_interpolation2d *interp,
    const double *xi,
    const double *yi,
    const double *zi,
    const size_t binsx,
    const size_t binsy,
    const enum coffe_interp2d_type interpolation_type
)
{
    if (binsx <= 0 || binsy <= 0){
        print_error(PROG_VALUE_ERROR);
        exit(EXIT_FAILURE);
    }
    const gsl_interp2d_type *T;
    switch(interpolation_type){
        case COFFE_INTERP2D_BILINEAR:
            T = gsl_interp2d_bilinear;
            break;
        case COFFE_INTERP2D_BICUBIC:
            T = gsl_interp2d_bicubic;
            break;
        default:
            T = gsl_interp2d_bicubic;
            break;
    }
    interp->spline
        = gsl_spline2d_alloc(T, binsx, binsy);
    interp->xaccel
        = gsl_interp_accel_alloc();
    interp->yaccel
        = gsl_interp_accel_alloc();
    gsl_spline2d_init(interp->spline, xi, yi, zi, binsx, binsy);
    return EXIT_SUCCESS;
}


int coffe_free_spline(
    coffe_interpolation *interp
)
{
    if (interp->spline != NULL)
        gsl_spline_free(interp->spline);
    if (interp->accel != NULL)
        gsl_interp_accel_free(interp->accel);

    interp->spline = NULL;
    interp->accel = NULL;

    return EXIT_SUCCESS;
}


int coffe_free_spline2d(
    coffe_interpolation2d *interp
)
{
    if (interp->spline != NULL)
        gsl_spline2d_free(interp->spline);
    if (interp->xaccel != NULL)
        gsl_interp_accel_free(interp->xaccel);
    if (interp->yaccel != NULL)
        gsl_interp_accel_free(interp->yaccel);

    interp->spline = NULL;
    interp->xaccel = NULL;
    interp->yaccel = NULL;

    return EXIT_SUCCESS;
}


int coffe_compare_ascending(
    const void *a,
    const void *b
)
{
    if (*(double *)a > *(double *)b) return 1;
    else if (*(double *)a < *(double *)b) return -1;
    else return 0;
}

int coffe_compare_descending(
    const void *a,
    const void *b
)
{
    if (*(double *)a > *(double *)b) return -1;
    else if (*(double *)a < *(double *)b) return 1;
    else return 0;
}


/**
    function describing the dark energy EOS
    for now only w(z) = w0 + wa*(1 - a),
    but others are straightforward to implement
**/

double coffe_dark_energy_eos(
    double z,
    const void *p
)
{
    const coffe_parameters_t *par = (coffe_parameters_t *)p;
    return par->w0 + par->wa * z / (1 + z);
}


/**
    function describing the galaxy bias analytically;
    by default, it is a fit of the form:
    b(z) = A * z^2 + B * z + C
    to the bias as described in GC of 1910.09273,
    table 3, but one can adapt it for their own needs
**/

double coffe_galaxy_bias(
    const double z
)
{
    const double A = -0.179886, B = 1.15485, C = 0.484564;
    return A * z * z + B * z + C;
}

double coffe_interp_spline(
    const coffe_interpolation *interp,
    double value
)
{
    return gsl_spline_eval(interp->spline, value, interp->accel);
}

double coffe_interp_spline2d(
    const coffe_interpolation2d *interp,
    const double value1,
    const double value2
)
{
    return gsl_spline2d_eval(interp->spline, value1, value2, interp->xaccel, interp->yaccel);
}

int coffe_parameters_free(
    coffe_parameters_t *par
)
{
    if (par->flag){
        #ifndef COFFE_CYTHON
        if (par->conf != NULL){
            config_destroy(par->conf);
            free(par->conf);
        }
        par->conf = NULL;
        #endif

        coffe_free_spline(&par->power_spectrum);
        coffe_free_spline(&par->power_spectrum_norm);
        coffe_free_spline2d(&par->power_spectrum2d);
        coffe_free_spline2d(&par->power_spectrum2d_norm);
        coffe_free_spline(&par->galaxy_bias1);
        coffe_free_spline(&par->galaxy_bias2);
        coffe_free_spline(&par->galaxy_bias3);
        coffe_free_spline(&par->galaxy_bias4);
        coffe_free_spline(&par->magnification_bias1);
        coffe_free_spline(&par->magnification_bias2);
        coffe_free_spline(&par->evolution_bias1);
        coffe_free_spline(&par->evolution_bias2);

        for (size_t i = 0; i<(size_t)par->type_bg_len; ++i){
            free(par->type_bg[i]);
        }
        free(par->type_bg);

        if (par->sep_len)
            free(par->sep);
        par->sep_len = 0;

        if (par->mu_len)
            free(par->mu);
        par->mu_len = 0;

        if (par->multipole_values_len)
            free(par->multipole_values);
        par->multipole_values_len = 0;

        if (par->z_mean_len)
            free(par->z_mean);
        par->z_mean_len = 0;

        if (par->deltaz_len)
            free(par->deltaz);
        par->deltaz_len = 0;

        if (par->fsky_len)
            free(par->fsky);
        par->fsky_len = 0;

        if (par->density1_len)
            free(par->density1);
        par->density1_len = 0;

        if (par->density2_len)
            free(par->density2);
        par->density2_len = 0;

        if (par->pixelsize_len)
            free(par->pixelsize);
        par->pixelsize_len = 0;

        if (par->zmin_len)
            free(par->zmin);

        if (par->zmax_len)
            free(par->zmax);
        par->zmax_len = 0;

        coffe_free_class_struct(&par->class_struct);

        coffe_free_fit_coefficients_array(&par->galaxy_bias1_coefficients);
        coffe_free_fit_coefficients_array(&par->galaxy_bias2_coefficients);
        coffe_free_fit_coefficients_array(&par->magnification_bias1_coefficients);
        coffe_free_fit_coefficients_array(&par->magnification_bias2_coefficients);

        par->flag = 0;
    }

    return EXIT_SUCCESS;
}

double coffe_resolution_window(
    double x
)
{
    return 3.0 * gsl_sf_bessel_j1(x) / x;
}


/**
    reads file <filename> into arrays <values1>..<valuesN>
    of length <len>
**/

int coffe_read_ncol(
    const char *filename,
    const size_t N,
    size_t *len,
    double **values,
    ...
)
{
    int error = 0;
    FILE *data = fopen(filename, "r");
    if (data == NULL){
        print_error_verbose(PROG_OPEN_ERROR, filename);
        return EXIT_FAILURE;
    }

    /* number of lines in the file */
    size_t n = 0;
    char c, temp_string[COFFE_MAX_STRLEN];
    char *temp_token;

    /* first we get how many lines are in the file, also including comments */
    while ((c=fgetc(data)) != EOF){
        if (c == '\n') ++n;
    }

    error = fseek(data, 0, SEEK_SET);
    if (error){
        print_error_verbose(PROG_POS_ERROR, filename);
        return EXIT_FAILURE;
    }

    /* get a line into temp_string */
    size_t counter = 0;
    for (size_t i = 0; i<n; ++i){
        fgets(temp_string, COFFE_MAX_STRLEN, data);
        if (temp_string[0] != '#') ++counter;
    }

    /* number of lines in the file */
    *len = counter;

    error = fseek(data, 0, SEEK_SET);
    if (error){
        print_error_verbose(PROG_POS_ERROR, filename);
        return EXIT_FAILURE;
    }

    /* we need to store <N> pointers to <values> */
    double **pointers = (double **)coffe_malloc(sizeof(double *) * N);

    for (size_t i = 0; i < N; ++i)
        /* allocating memory for <values> */
        pointers[i] = (double *)coffe_malloc(sizeof(double) * counter);

    /* kind of recycling <counter> here */
    counter = 0;

    for (size_t i = 0; i < n; ++i){
        fgets(temp_string, COFFE_MAX_STRLEN, data);
        if (temp_string[0] != '#'){
            temp_token = strtok(temp_string, ",\t: ");
            for (size_t j = 0; j < N; ++j){
                /* <pointers[j]> should hopefully be the address of <values> */
                if (temp_token != NULL)
                    *(pointers[j] + counter) = atof(temp_token);
                temp_token = strtok(NULL, ",\t: ");
            }
            ++counter;
        }
    }

    /* here we need to do the varg trick */
    va_list args;

    /* now we re-assign the pointers */
    va_start(args, values);
    for (size_t i = 0; i < N; ++i){
        *values = pointers[i];
        values = va_arg(args, double **);
    }
    va_end(args);

    /* memory cleanup of <pointers> */
    for (size_t i = 0; i < N; ++i)
        pointers[i] = NULL;
    free(pointers);

    error = fclose(data);
    if (error){
        print_error_verbose(PROG_CLOSE_ERROR, filename);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

/**
    multiplies every item of a double `array` of `size` by a double `factor`
    MUST be alloc'd beforehand
**/
void coffe_rescale_array(
    double *array,
    const size_t size,
    const double factor
)
{
    if (array != NULL){
        for (size_t i = 0; i < size; ++i) array[i] *= factor;
    }
    else{
        fprintf(
            stderr,
            "ERROR: NULL encountered in %s, exiting, fix your code!\n",
            __func__
        );
        exit(EXIT_FAILURE);
    }
}

/**
    multiplies two arrays, `array1` and `array2`, of same `size`,
    and stores the result into `array_out`
    MUST be alloc'd beforehand
**/
void coffe_multiply_array(
    double *array_out,
    const double *array1,
    const double *array2,
    const size_t size
)
{
    if (
        array_out != NULL &&
        array1 != NULL &&
        array2 != NULL
    ){
        for (size_t i = 0; i < size; ++i) array_out[i] = array1[i] * array2[i];
    }
    else{
        fprintf(
            stderr,
            "ERROR: NULL encountered in %s, exiting, fix your code!\n",
            __func__
        );
        exit(EXIT_FAILURE);
    }
}

/**
    computes `array1` * `array2`^`power`, of same `size`,
    and stores the result into `array_out`
    MUST be alloc'd beforehand
**/
void coffe_multiply_power_array(
    double *array_out,
    const double *array1,
    const double *array2,
    const size_t size,
    const double power
)
{
    if (
        array_out != NULL &&
        array1 != NULL &&
        array2 != NULL
    ){
        for (size_t i = 0; i < size; ++i) array_out[i] = array1[i] * pow(array2[i], power);
    }
    else{
        fprintf(
            stderr,
            "ERROR: NULL encountered in %s, exiting, fix your code!\n",
            __func__
        );
        exit(EXIT_FAILURE);
    }
}


/**
    integrates any 1D function `func` with arbitrary parameters `parameters`
    between `a` and `b`, with relative precision `prec`, and returns `result`
**/

double coffe_integrate_1d_prec(
    double (*func)(
        double,
#ifdef HAVE_DOUBLE_EXPONENTIAL
        const void*
#else
        void*
#endif
    ),
    const void *parameters,
    const double a,
    const double b,
    const double precision
)
{
    double result, error;

#ifdef HAVE_DOUBLE_EXPONENTIAL
    result = tanhsinh_quad(
        func,
        parameters,
        a,
        b,
        0.,
        &error,
        NULL
    );
#else

    gsl_function integrand;
    integrand.params = (void *)parameters;
    integrand.function = func;

    gsl_integration_workspace *wspace =
        gsl_integration_workspace_alloc(COFFE_MAX_INTSPACE);

    gsl_integration_qag(
        &integrand,
        a,
        b,
        0,
        precision,
        COFFE_MAX_INTSPACE,
        GSL_INTEG_GAUSS61,
        wspace,
        &result,
        &error
    );

    gsl_integration_workspace_free(wspace);
#endif

    return result;
}


double coffe_integrate_1d(
    double (*func)(
        double,
#ifdef HAVE_DOUBLE_EXPONENTIAL
        const void*
#else
        void*
#endif
    ),
    const void *parameters,
    const double a,
    const double b
)
{
    return coffe_integrate_1d_prec(
        func,
        parameters,
        a,
        b,
        1e-5
    );
}


/**
    integrate a multidimensional function (on the interval [0, 1] only!)
    using either GSL or Cuba (if available)
**/

double coffe_integrate_multidimensional(
#ifdef HAVE_CUBA
    int (*func)(
        const int *,
        const cubareal *,
        const int *,
        cubareal *,
        void *
    ),
#else
    double (*func)(
        double *,
        size_t,
        void *
    ),
#endif
    const void *parameters,
    const int integration_method,
    const int dims,
    const int integration_bins
)
{
#ifdef HAVE_CUBA

    int nregions, neval, fail;
    double result[1], error[1], prob[1];

    Cuhre(
        dims,
        1,
        func,
        (void *)parameters,
        1,
        5e-4,
        0,
        0,
        1,
        integration_bins,
        7,
        NULL,
        NULL,
        &nregions, &neval, &fail, result, error, prob
    );

    return result[0];

#else

    double result;
    gsl_monte_function integrand;
    integrand.f = func;
    integrand.dim = dims;
    integrand.params = (void *)parameters;
    gsl_rng_env_setup();
    const gsl_rng_type *rng = gsl_rng_default;
    gsl_rng *random = gsl_rng_alloc(rng);
    double lower[dims];
    double upper[dims];
    double error;
    for (int i = 0; i < dims; ++i){
        lower[i] = 0.0;
        upper[i] = 1.0;
    }
    switch (integration_method){
        case 0:{
            gsl_monte_plain_state *state =
                gsl_monte_plain_alloc(dims);
            gsl_monte_plain_integrate(
                &integrand, lower, upper,
                dims, integration_bins, random,
                state,
                &result, &error
            );
            gsl_monte_plain_free(state);
            break;
        }
        case 1:{
            gsl_monte_miser_state *state =
                gsl_monte_miser_alloc(dims);
            gsl_monte_miser_integrate(
                &integrand, lower, upper,
                dims, integration_bins, random,
                state,
                &result, &error
            );
            gsl_monte_miser_free(state);
            break;
        }
        case 2:{
            gsl_monte_vegas_state *state =
                gsl_monte_vegas_alloc(dims);
            gsl_monte_vegas_integrate(
                &integrand, lower, upper,
                dims, integration_bins, random,
                state,
                &result, &error
            );
            gsl_monte_vegas_free(state);
            break;
        }
    }
    gsl_rng_free(random);

    return result;
#endif
}

double *coffe_generate_range(
    const double xmin,
    const double xmax,
    const size_t steps
)
{
    /**
        Works the same as Python's range.
        User is responsible for memory cleanup.
    **/
    if (xmin >= xmax){
        fprintf(
            stderr, "ERROR: xmin (%e) larger than xmax (%e)\n",
            xmin, xmax
        );
        exit(EXIT_FAILURE);
    }
    double *result = (double *)coffe_malloc(sizeof(double) * steps);
    for (size_t i = 0; i < steps; ++i){
        result[i] = xmin + (xmax - xmin) * i / steps;
    }
    return result;
}

/* 1 if equal, else 0 */
int coffe_approx_equal(
    const double a,
    const double b,
    const double rel_epsilon,
    const double abs_epsilon
)
{
    if (a != 0 && b != 0){
        return (fabs(a - b) / fabs(a) < rel_epsilon) ? 1 : 0;
    }
    else{
        return (fabs(a - b) < abs_epsilon) ? 1 : 0;
    }
}

int coffe_new_spline(
    coffe_interpolation *input
)
{
    input->spline = NULL;
    input->accel = NULL;

    return EXIT_SUCCESS;
}

int coffe_new_spline2d(
    coffe_interpolation2d *input
)
{
    input->spline = NULL;
    input->xaccel = NULL;
    input->yaccel = NULL;

    return EXIT_SUCCESS;
}

int coffe_new_fit_coefficients(
    coffe_fit_coefficients_t *input
)
{
    input->coefficients = NULL;
    input->size = 0;

    return EXIT_SUCCESS;
}


int coffe_free_fit_coefficients(
    coffe_fit_coefficients_t *input
)
{
    if (input->size)
        free(input->coefficients);
    input->size = 0;
    input->coefficients = NULL;

    return EXIT_SUCCESS;
}


int coffe_new_fit_coefficients_array(
    coffe_fit_coefficients_array_t *input
)
{
    input->size = 0;
    input->array = NULL;

    return EXIT_SUCCESS;
}



int coffe_free_fit_coefficients_array(
    coffe_fit_coefficients_array_t *input
)
{
    if (input->size){
        for (size_t i = 0; i < input->size; ++i)
            coffe_free_fit_coefficients(&input->array[i]);
    }
    input->size = 0;

    return EXIT_SUCCESS;
}


int coffe_new_class_struct(
    coffe_class_struct_t *input
)
{
    input->file_content = NULL;
    input->background = NULL;
    input->thermodynamics = NULL;
    input->perturb = NULL;
    input->primordial = NULL;
    input->nonlinear = NULL;
    input->transfer = NULL;
    input->spectra = NULL;

    return EXIT_SUCCESS;
}

int coffe_free_class_struct(
    coffe_class_struct_t *input
)
{
#ifdef HAVE_CLASS
    if (input->file_content){
        parser_free((struct file_content *)input->file_content);
        free(input->file_content);
    }

    if (input->background){
        background_free((struct background *)input->background);
        free(input->background);
    }

    if (input->thermodynamics){
        thermodynamics_free((struct thermo *)input->thermodynamics);
        free(input->thermodynamics);
    }

    if (input->perturb){
        perturb_free((struct perturbs *)input->perturb);
        free(input->perturb);
    }

    if (input->primordial){
        primordial_free((struct primordial *)input->primordial);
        free(input->primordial);
    }

    if (input->nonlinear){
        nonlinear_free((struct nonlinear *)input->nonlinear);
        free(input->nonlinear);
    }

    if (input->transfer){
        transfer_free((struct transfers *)input->transfer);
        free(input->transfer);
    }

    if (input->spectra){
        spectra_free((struct spectra *)input->spectra);
        free(input->spectra);
    }
#endif
    input->background = NULL;
    input->thermodynamics = NULL;
    input->perturb = NULL;
    input->primordial = NULL;
    input->transfer = NULL;
    input->spectra = NULL;
    input->file_content = NULL;
    input->nonlinear = NULL;

    return EXIT_SUCCESS;
}



/**
    computes the integral over 4 Legendre polynomials of degrees n, m, a, and b
**/
double coffe_legendre_integral(
    int n,
    int m,
    int a,
    int b
)
{
    /* shortcut: if n + m + a + b is odd, the integral is identically zero */
    if ((n + m + a + b) % 2 != 0)
        return 0;

    const int lower_limit = (abs(n - m) > abs(a - b) ? abs(n - m) : abs(a - b));
    const int upper_limit = (n + m < a + b ? n + m : a + b);

    double result = 0;

    for (int l = lower_limit; l <= upper_limit; ++l){
        result +=
            (2. * l + 1)
           *pow(gsl_sf_coupling_3j(2 * n, 2 * m, 2 * l, 0, 0, 0), 2)
           *pow(gsl_sf_coupling_3j(2 * a, 2 * b, 2 * l, 0, 0, 0), 2);
    }

    return 2 * result;
}



/**
    calculates the Kronecker delta
**/
int coffe_kronecker_delta(
    int i,
    int j
)
{
    return (i == j ? 1 : 0);
}



/**
    calculates (-1)^m
**/
int coffe_sign(
    int m
)
{
    return (abs(m) % 2 == 0 ? 1 : -1);
}
