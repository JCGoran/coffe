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
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <gsl/gsl_integration.h>

#include "common.h"
#include "parser.h"
#include "errors.h"

#ifdef HAVE_CLASS
#include "class.h"
#endif

static int parse_double(
    config_t *conf,
    const char *setting,
    double *value,
    const int safe
)
{
    int error = config_lookup_float(conf, setting, value);
    if (safe == COFFE_TRUE){
        if (error != 1){
            fprintf(stderr,
                "ERROR: cannot parse value of setting %s\n",
            setting);
            exit(EXIT_FAILURE);
        }
    }
    else{
        if (error != 1){
            fprintf(stderr,
                "ERROR: cannot parse value of setting %s,"
                " continuing...\n", setting);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

static int parse_int(
    config_t *conf,
    const char *setting,
    int *value,
    const int safe
)
{
    int error = config_lookup_int(conf, setting, value);
    if (safe == COFFE_TRUE){
        if (error != 1){
            fprintf(stderr, "ERROR: cannot parse value of setting %s\n", setting);
            exit(EXIT_FAILURE);
        }
    }
    else{
        if (error != 1){
            fprintf(stderr,
                    "ERROR: cannot parse value of setting %s,"
                    " continuing...\n", setting);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

static int parse_string(
    config_t *conf,
    const char *setting,
    char value[],
    const int safe
)
{
    const char *temp_value;
    int error = config_lookup_string(conf, setting, &temp_value);
    strcpy(value, (char *)temp_value);
    if (safe == COFFE_TRUE){
        if (error != 1){
            fprintf(stderr,
                "ERROR: cannot parse value of setting %s\n", setting);
            exit(EXIT_FAILURE);
        }
    }
    else{
        if (error != 1){
            fprintf(stderr, "ERROR: cannot parse value of setting %s,"
                   " continuing...\n", setting);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}


/**
    parses an array of strings from the setting <setting>
    in config file <conf> into <values> with length <values_len>
**/

static int parse_string_array(
    config_t *conf,
    const char *setting,
    char ***values,
    int *values_len
)
{
    config_setting_t *type = config_lookup(conf, setting);
    if (type == NULL){
        print_error(PROG_NULL_PARAM);
        return EXIT_FAILURE;
    }
    *values_len = config_setting_length(type);

    *values = (char **)malloc(sizeof(char *) * *values_len);
    if (*values == NULL){
        print_error(PROG_ALLOC_ERROR);
        return EXIT_FAILURE;
    }

    const char *temp_setting;

    for (int i = 0; i<*values_len; ++i){
        temp_setting = config_setting_get_string_elem(type, i);
        (*values)[i] =
            (char *)coffe_malloc(sizeof(char)*(strlen(temp_setting) + 1));
        strcpy((char *)(*values)[i], temp_setting);
    }
    return EXIT_SUCCESS;
}


/**
    parses an array of integers from the setting <setting>
    in config file <conf> into <values> with length <values_len>
**/

static int parse_int_array(
    config_t *conf,
    const char *setting,
    int **values,
    int *values_len
)
{
    config_setting_t *type = config_lookup(conf, setting);
    if (type == NULL){
        print_error(PROG_NULL_PARAM);
        return EXIT_FAILURE;
    }
    *values_len = config_setting_length(type);

    *values = (int *)coffe_malloc(sizeof(int) * *values_len);

    int temp_value;

    for (int i = 0; i<*values_len; ++i){
        temp_value = config_setting_get_int_elem(type, i);
        *(*values + i) = temp_value;
    }
    return EXIT_SUCCESS;
}


/**
    parses an array of floats from the setting <setting>
    in config file <conf> into <values> with length <values_len>
**/

static int parse_double_array(
    config_t *conf,
    const char *setting,
    double **values,
    int *values_len
)
{
    config_setting_t *type = config_lookup(conf, setting);
    if (type == NULL){
        print_error(PROG_NULL_PARAM);
        return EXIT_FAILURE;
    }
    *values_len = config_setting_length(type);

    *values = (double *)coffe_malloc(sizeof(double) * *values_len);

    double temp_value;

    for (int i = 0; i<*values_len; ++i){
        temp_value = config_setting_get_float_elem(type, i);
        *(*values + i) = temp_value;
    }
    return EXIT_SUCCESS;
}


static int parse_bias_default(
    double value,
    struct coffe_interpolation *spline,
    int method
)
{
    double redshifts[] = {0, 25, 50, 75, 100};
    double values[] = {value, value, value, value, value};
    coffe_init_spline(
        spline, redshifts, values,
        sizeof(redshifts)/sizeof(*redshifts),
        method
    );

    return EXIT_SUCCESS;
}


/**
    parses the value of bias given by <setting> from config <conf>
    into coffe_interpolation <spline> depending on value of <flag>
**/

static int parse_bias(
    config_t *conf,
    char *setting,
    char *filename,
    char *setting_single,
    struct coffe_interpolation *spline,
    int method,
    int flag
)
{
    if (flag == COFFE_TRUE){
        double *redshifts, *values;
        size_t len;
        parse_string(conf, setting, filename, COFFE_TRUE);
        read_2col(
            filename,
            &redshifts,
            &values,
            &len
        );
        coffe_init_spline(spline, redshifts, values, len, method);
        free(redshifts);
        free(values);
    }
    else{
        double value;
        parse_double(conf, setting_single, &value, COFFE_TRUE);
        parse_bias_default(value, spline, method);
    }

    return EXIT_SUCCESS;
}


/**
    parsing the power spectrum from an external source (only CLASS for now)
**/
#ifdef HAVE_CLASS
static int parse_external_power_spectrum(
    struct coffe_parameters_t *par
)
{
    struct precision ppr;
    struct background pba;
    struct thermo pth;
    struct perturbs ppt;
    struct primordial ppm;
    struct spectra psp;
    struct nonlinear pnl;
    struct output pop;
    struct transfers ptr;
    struct lensing ple;
    ErrorMsg errmsg;

    clock_t class_start, class_end;

    if (par->verbose)
        printf("Launching CLASS...\n");

    class_start = clock();
    struct file_content fc;

    size_t class_parameters_len = 18, counter = 0;

    parser_init(&fc, class_parameters_len, "", errmsg);

    /*
        not sure which values are actually necessary
        (maybe give the user the ability to read all of them?)
    */

    sprintf(fc.name[counter], "h");
    sprintf(fc.value[counter], "%e", par->h);
    ++counter;

    sprintf(fc.name[counter], "T_cmb");
    sprintf(fc.value[counter], "%e", 2.7255);
    ++counter;

    sprintf(fc.name[counter], "Omega_b");
    sprintf(fc.value[counter], "%e", par->Omega0_baryon);
    ++counter;

    sprintf(fc.name[counter], "N_ur");
    sprintf(fc.value[counter], "%e", 3.046);
    ++counter;

    sprintf(fc.name[counter], "Omega_cdm");
    sprintf(fc.value[counter], "%e", par->Omega0_cdm);
    ++counter;

    sprintf(fc.name[counter], "Omega_k");
    sprintf(fc.value[counter], "%e", 0.0);
    ++counter;

    sprintf(fc.name[counter], "w0_fld");
    sprintf(fc.value[counter], "%e", par->w0);
    ++counter;

    sprintf(fc.name[counter], "wa_fld");
    sprintf(fc.value[counter], "%e", par->wa);
    ++counter;

    sprintf(fc.name[counter], "output");
    sprintf(fc.value[counter], "mPk");
    ++counter;

    sprintf(fc.name[counter], "gauge");
    sprintf(fc.value[counter], "synchronous");
    ++counter;

    sprintf(fc.name[counter], "P_k_ini type");
    sprintf(fc.value[counter], "analytic_Pk");
    ++counter;

    sprintf(fc.name[counter], "k_pivot");
    sprintf(fc.value[counter], "%e", par->k_pivot);
    ++counter;

    sprintf(fc.name[counter], "ln10^{10}A_s");
    sprintf(fc.value[counter], "%e", par->ln_10_pow_10_A_s);
    ++counter;

    sprintf(fc.name[counter], "n_s");
    sprintf(fc.value[counter], "%e", par->n_s);
    ++counter;

    sprintf(fc.name[counter], "alpha_s");
    sprintf(fc.value[counter], "%e", 0.0);
    ++counter;

    sprintf(fc.name[counter], "k_min_tau0");
    sprintf(fc.value[counter], "%e", 0.002);
    ++counter;

    sprintf(fc.name[counter], "P_k_max_h/Mpc");
    sprintf(fc.value[counter], "%e", par->k_max);
    ++counter;

    sprintf(fc.name[counter], "z_pk");
    sprintf(fc.value[counter], "%e", 0.0);
    ++counter;

    if (
        input_init(
            &fc, &ppr, &pba, &pth, &ppt, &ptr,
            &ppm, &psp, &pnl, &ple, &pop, errmsg
        ) == _FAILURE_
    ){
        fprintf(stderr, "\n\nError running input_init\n=>%s\n", errmsg);
        exit(EXIT_FAILURE);
    }

    /* the main CLASS sequence */

    background_init(&ppr, &pba);
    thermodynamics_init(&ppr, &pba, &pth);
    perturb_init(&ppr, &pba, &pth, &ppt);
    primordial_init(&ppr, &ppt, &ppm);
    nonlinear_init(&ppr, &pba, &pth, &ppt, &ppm, &pnl);
    transfer_init(&ppr, &pba, &pth, &ppt, &pnl, &ptr);
    spectra_init(&ppr, &pba, &ppt, &ppm, &pnl, &ptr, &psp);

    class_end = clock();
    if (par->verbose)
        printf(
            "CLASS finished in %.2f s\n",
            (double)(class_end - class_start) / CLOCKS_PER_SEC
        );

    double *k =(double *)coffe_malloc(sizeof(double)*psp.ln_k_size);
    double *pk =(double *)coffe_malloc(sizeof(double)*psp.ln_k_size);
    size_t pk_len = psp.ln_k_size;

    /* rescaling as CLASS internally uses units of 1/Mpc */
    for (size_t i = 0; i<(size_t)psp.ln_k_size; ++i){
        k[i] = exp(psp.ln_k[i]) / par->h;
        pk[i] = exp(psp.ln_pk_l[i]) * pow(par->h, 3);
    }
    coffe_init_spline(&par->power_spectrum, k, pk, pk_len, par->interp_method);

    free(k);
    free(pk);

    /* freeing the memory of CLASS structures */

    background_free(&pba);
    thermodynamics_free(&pth);
    perturb_free(&ppt);
    primordial_free(&ppm);
    nonlinear_free(&pnl);
    transfer_free(&ptr);
    spectra_free(&psp);
    parser_free(&fc);

    return EXIT_SUCCESS;
}
#endif

int coffe_parse_default_parameters(
    struct coffe_parameters_t *par
)
{
    par->Omega0_cdm = 0.25793;
    par->Omega0_baryon = 0.0484;
    par->Omega0_gamma = 9e-5;
    par->Omega0_de = 1 - (par->Omega0_cdm + par->Omega0_baryon + par->Omega0_gamma);
    par->w0 = -1.0;
    par->wa = 0.0;
    #ifdef HAVE_CLASS
    par->have_class = 0;
    par->h = 0.6781;
    par->k_pivot = 0.05;
    par->ln_10_pow_10_A_s = 3.062;
    par->n_s = 0.9677;
    #endif

    snprintf(par->file_sep, COFFE_MAX_STRLEN, "\0");
    const double separations[] = {10., 20., 40., 100., 150.};
    par->sep = coffe_malloc(
        sizeof(double) * sizeof(separations) / sizeof(*separations)
    );
    for (int i = 0; i < sizeof(separations) / sizeof(*separations); ++i)
        par->sep[i] = separations[i];
    par->sep_len = sizeof(separations) / sizeof(*separations);
    par->interp_method = 5;
    snprintf(par->file_power_spectrum, COFFE_MAX_STRLEN, "\0");
    /*
        turns out you can do this
        https://stackoverflow.com/a/1662202
    */
    const double k[] = {
        #include "WAVENUMBER_HEADER.dat"
    };
    const double pk[] = {
        #include "POWER_SPECTRUM_HEADER.dat"
    };
    coffe_init_spline(
        &par->power_spectrum,
        (double *)k, (double *)pk, sizeof(k) / sizeof(*k),
        par->interp_method
    );
    par->k_min = par->power_spectrum.spline->x[0];
    par->k_max = par->power_spectrum.spline->x[
        par->power_spectrum.spline->size - 1
    ];
    {
        size_t len = par->power_spectrum.spline->size;
        double *k_norm =
            (double *)coffe_malloc(sizeof(double)*len);
        double *pk_norm =
            (double *)coffe_malloc(sizeof(double)*len);
        for (size_t i = 0; i<len; ++i){
            k_norm[i] = par->power_spectrum.spline->x[i] / COFFE_H0;
            pk_norm[i] = par->power_spectrum.spline->y[i] * pow(COFFE_H0, 3);
        }
        coffe_init_spline(
            &par->power_spectrum_norm,
            k_norm, pk_norm, len,
            par->interp_method
        );
        par->k_min_norm = par->k_min / COFFE_H0;
        par->k_max_norm = par->k_max / COFFE_H0;
        free(k_norm);
        free(pk_norm);
    }

    parse_bias_default(
        1.0, &par->matter_bias1, par->interp_method
    );
    par->read_matter_bias1 = 0;
    snprintf(par->file_matter_bias1, COFFE_MAX_STRLEN, "\0");

    parse_bias_default(
        1.0, &par->matter_bias2, par->interp_method
    );
    par->read_matter_bias2 = 0;
    snprintf(par->file_matter_bias2, COFFE_MAX_STRLEN, "\0");

    parse_bias_default(
        0.0, &par->magnification_bias1, par->interp_method
    );
    par->read_magnification_bias1 = 0;
    snprintf(par->file_magnification_bias1, COFFE_MAX_STRLEN, "\0");

    parse_bias_default(
        0.0, &par->magnification_bias2, par->interp_method
    );
    par->read_magnification_bias2 = 0;
    snprintf(par->file_magnification_bias2, COFFE_MAX_STRLEN, "\0");

    parse_bias_default(
        0.0, &par->evolution_bias1, par->interp_method
    );
    par->read_evolution_bias1 = 0;
    snprintf(par->file_evolution_bias1, COFFE_MAX_STRLEN, "\0");

    parse_bias_default(
        0.0, &par->evolution_bias2, par->interp_method
    );
    par->read_evolution_bias2 = 0;
    snprintf(par->file_evolution_bias2, COFFE_MAX_STRLEN, "\0");

    par->output_type = 2;
    par->covariance_density = NULL;
    par->covariance_density_len = 0;
    par->covariance_z_mean = NULL;
    par->covariance_z_mean_len = 0;
    par->covariance_deltaz = NULL;
    par->covariance_deltaz_len = 0;
    par->covariance_pixelsize = 0.0;
    par->covariance_zmin = NULL;
    par->covariance_zmin_len = 0;
    par->covariance_zmax = NULL;
    par->covariance_zmax_len = 0;
    par->covariance_minimum_separation = 0.0;

    snprintf(par->output_path, COFFE_MAX_STRLEN, "./");
    snprintf(par->output_prefix, COFFE_MAX_STRLEN, "$TIME");
    snprintf(
        par->timestamp,
        COFFE_MAX_STRLEN,
        "%s",
        coffe_get_time()
    );

    par->correlation_contrib.den = 1;
    par->correlation_contrib.rsd = 1;
    par->correlation_contrib.d1 = 0;
    par->correlation_contrib.d2 = 0;
    par->correlation_contrib.g1 = 0;
    par->correlation_contrib.g2 = 0;
    par->correlation_contrib.g3 = 0;
    par->correlation_contrib.len = 0;
    par->correlation_contrib.g4 = 0;
    par->correlation_contrib.g5 = 0;

    par->mu = NULL;
    par->mu_len = 0;

    par->divergent = 0;

    par->nonzero_terms[0].n = 0, par->nonzero_terms[0].l = 0;
    par->nonzero_terms[1].n = 0, par->nonzero_terms[1].l = 2;
    par->nonzero_terms[2].n = 0, par->nonzero_terms[2].l = 4;
    par->nonzero_terms[3].n = 1, par->nonzero_terms[3].l = 1;
    par->nonzero_terms[4].n = 1, par->nonzero_terms[4].l = 3;
    par->nonzero_terms[5].n = 2, par->nonzero_terms[5].l = 0;
    par->nonzero_terms[6].n = 2, par->nonzero_terms[6].l = 2;
    par->nonzero_terms[7].n = 3, par->nonzero_terms[7].l = 1;

    par->type_bg = NULL;
    par->type_bg_len = 0;

    par->background_bins = 10000;
    par->bessel_bins = 10000;

    par->z_mean = 1.0;
    par->deltaz = 0.2;

    par->z_min = 0.9;
    par->z_max = 1.1;

    const double multipoles[] = {0, 2, 4};
    par->multipole_values = coffe_malloc(
        sizeof(double) * sizeof(multipoles) / sizeof(*multipoles)
    );
    for (int i = 0; i < sizeof(multipoles) / sizeof(*multipoles); ++i)
        par->multipole_values[i] = multipoles[i];
    par->multipole_values_len = sizeof(multipoles) / sizeof(*multipoles);

    par->flatsky = 0;

    par->theta_len = 0;

    par->verbose = 0;

    par->conf = NULL;

    return EXIT_SUCCESS;
}


/**
    parses all the settings from the input file
    (given by argv[1]) into the structure <par>
**/

int coffe_parser_init(
    char *filename,
    struct coffe_parameters_t *par
)
{
    coffe_parse_default_parameters(par);
    config_t *conf = (config_t *)coffe_malloc(sizeof(config_t));

    config_init(conf);

    if(!config_read_file(conf, filename))
    {
        fprintf(stderr, "ERROR: Cannot read the parameter file.\n");
        exit(EXIT_FAILURE);
    }

    clock_t start, end;

    start = clock();

    config_set_auto_convert(conf, CONFIG_TRUE);

    parse_int(conf, "verbose", &par->verbose, COFFE_TRUE);

    if (par->verbose){
        printf(COFFE_LOGO);
        printf("Parsing settings file \"%s\"...\n", filename);
    }

    parse_int(conf, "output_type", &par->output_type, COFFE_TRUE);
    parse_string_array(conf, "output_background", &par->type_bg, &par->type_bg_len);
    parse_int(conf, "background_sampling", &par->background_bins, COFFE_TRUE);

    /* cosmological parameters */
    parse_double(conf, "omega_cdm", &par->Omega0_cdm, COFFE_TRUE);
    parse_double(conf, "omega_baryon", &par->Omega0_baryon, COFFE_TRUE);
    parse_double(conf, "omega_gamma", &par->Omega0_gamma, COFFE_TRUE);

    par->Omega0_de = 1. - (par->Omega0_cdm + par->Omega0_baryon + par->Omega0_gamma);

    /* mean redshift */
    if (
        par->output_type == 0 ||
        par->output_type == 1 ||
        par->output_type == 2 ||
        par->output_type == 6
    ){
        parse_double(conf, "z_mean", &par->z_mean, COFFE_TRUE);
        if (par->z_mean <= 0){
            print_error_verbose(PROG_VALUE_ERROR, "z_mean");
            exit(EXIT_FAILURE);
        }
    }

    /* width of redshift bin */
    if (par->output_type == 1 || par->output_type == 2){
        parse_double(conf, "deltaz", &par->deltaz, COFFE_TRUE);
        if (par->deltaz <= 0){
            print_error_verbose(PROG_VALUE_ERROR, "deltaz");
            exit(EXIT_FAILURE);
        }
        /* safety check for the range of deltaz */
        if (par->deltaz > par->z_mean){
            fprintf(
                stderr,
                "ERROR: z_mean cannot be smaller than deltaz!\n"
            );
            exit(EXIT_FAILURE);
        }
    }

    /* range of integration for redshift averaged multipoles */
    if (par->output_type == 3){
        parse_double(conf, "z_min", &par->z_min, COFFE_TRUE);
        parse_double(conf, "z_max", &par->z_max, COFFE_TRUE);
    }

    /* the interpolation method for GSL */
    parse_int(conf, "interpolation", &par->interp_method, COFFE_FALSE);

    /* the cosine of the angle for the full sky correlation function */
    if (par->output_type == 1){
        parse_double_array(conf, "mu", &par->mu, &par->mu_len);
    }

    /* the multipoles of the correlation function */
    if (
        par->output_type == 2 ||
        par->output_type == 3 ||
        par->output_type == 4 ||
        par->output_type == 5
    ){
        parse_int_array(conf, "multipoles", &par->multipole_values, &par->multipole_values_len);
    }

    /* the custom separations for the ang/full correlation function or multipoles */
    if (par->output_type == 1 || par->output_type == 2 || par->output_type == 3){
        parse_string(conf, "input_separations", par->file_sep, COFFE_TRUE);
        read_1col(
            par->file_sep,
            &par->sep,
            &par->sep_len
        );
    }

    /* number of points to sample the integral of the Bessel function */
    parse_int(conf, "bessel_sampling", &par->bessel_bins, COFFE_TRUE);

#ifndef HAVE_CUBA
    /* parsing the integration method */
    parse_int(conf, "integration_method", &par->integration_method, COFFE_TRUE);
    if (par->integration_method < 0 || par->integration_method > 2){
        print_error_verbose(PROG_VALUE_ERROR, "integration_method");
        exit(EXIT_FAILURE);
    }
#endif

    /* number of points for the 2-3-4D integration */
    parse_int(conf, "integration_sampling", &par->integration_bins, COFFE_TRUE);

    /* parsing the w parameter */
    parse_double(conf, "w0", &par->w0, COFFE_TRUE);
    parse_double(conf, "wa", &par->wa, COFFE_TRUE);

    /* parsing the matter bias */
    parse_int(conf, "read_matter_bias1", &par->read_matter_bias1, COFFE_FALSE);
    parse_bias(
        conf,
        "input_matter_bias1",
        par->file_matter_bias1,
        "matter_bias1",
        &par->matter_bias1,
        par->interp_method,
        par->read_matter_bias1
    );

    parse_int(conf, "read_matter_bias2", &par->read_matter_bias2, COFFE_FALSE);
    parse_bias(
        conf,
        "input_matter_bias2",
        par->file_matter_bias2,
        "matter_bias2",
        &par->matter_bias2,
        par->interp_method,
        par->read_matter_bias2
    );

    /* parsing the magnification bias (s) */
    parse_int(conf, "read_magnification_bias1", &par->read_magnification_bias1, COFFE_FALSE);
    parse_bias(
        conf,
        "input_magnification_bias1",
        par->file_magnification_bias1,
        "magnification_bias1",
        &par->magnification_bias1,
        par->interp_method,
        par->read_magnification_bias1
    );

    parse_int(conf, "read_magnification_bias2", &par->read_magnification_bias2, COFFE_FALSE);
    parse_bias(
        conf,
        "input_magnification_bias2",
        par->file_magnification_bias2,
        "magnification_bias2",
        &par->magnification_bias2,
        par->interp_method,
        par->read_magnification_bias2
    );

    /* parsing the evolution_bias (f_evo) */
    parse_int(conf, "read_evolution_bias1", &par->read_evolution_bias1, COFFE_FALSE);
    parse_bias(
        conf,
        "input_evolution_bias1",
        par->file_evolution_bias1,
        "evolution_bias1",
        &par->evolution_bias1,
        par->interp_method,
        par->read_evolution_bias1
    );

    parse_int(conf, "read_evolution_bias2", &par->read_evolution_bias2, COFFE_FALSE);
    parse_bias(
        conf,
        "input_evolution_bias2",
        par->file_evolution_bias2,
        "evolution_bias2",
        &par->evolution_bias2,
        par->interp_method,
        par->read_evolution_bias2
    );

    /* parsing the covariance parameters */
    if (par->output_type == 4 || par->output_type == 5){
        parse_double_array(
            conf,
            "covariance_density",
            &par->covariance_density,
            &par->covariance_density_len
        );
        parse_double_array(
            conf,
            "covariance_fsky",
            &par->covariance_fsky,
            &par->covariance_fsky_len
        );
        parse_double(
            conf,
            "covariance_pixelsize",
            &par->covariance_pixelsize,
            COFFE_TRUE
        );
        parse_double(
            conf,
            "covariance_minimum_separation",
            &par->covariance_minimum_separation,
            COFFE_TRUE
        );
    }

    if (par->output_type == 4){
        parse_double_array(
            conf,
            "covariance_z_mean",
            &par->covariance_z_mean,
            &par->covariance_z_mean_len
        );
        parse_double_array(
            conf,
            "covariance_deltaz",
            &par->covariance_deltaz,
            &par->covariance_deltaz_len
        );
        if (
            par->covariance_density_len != par->covariance_fsky_len ||
            par->covariance_density_len != par->covariance_z_mean_len ||
            par->covariance_density_len != par->covariance_deltaz_len
        ){
            fprintf(
                stderr,
                "ERROR: covariance parameters have mismatching lengths; "
                "please ensure they have the same length!\n");
            exit(EXIT_FAILURE);
        }
    }

    if (par->output_type == 5){
        parse_double_array(
            conf,
            "covariance_zmin",
            &par->covariance_zmin,
            &par->covariance_zmin_len
        );
        parse_double_array(
            conf,
            "covariance_zmax",
            &par->covariance_zmax,
            &par->covariance_zmax_len
        );
        if (
            par->covariance_density_len != par->covariance_fsky_len ||
            par->covariance_density_len != par->covariance_zmin_len ||
            par->covariance_density_len != par->covariance_zmax_len
        ){
            fprintf(
                stderr,
                "ERROR: covariance parameters have mismatching lengths; "
                "please ensure they have the same length!\n");
            exit(EXIT_FAILURE);
        }
    }

    /* if we just want the angular correlation function */
    if (par->output_type == 0){
        parse_int(
            conf,
            "theta_sampling",
            &par->theta_len,
            COFFE_TRUE
        );
    }

    for (int i = 0; i<9; ++i){
        par->nonzero_terms[i].l = -1, par->nonzero_terms[i].n = -1;
    }

    /* parsing the contributions to the correlation function and covariance */
    char **correlation_contributions = NULL;
    int correlation_contributions_len = 0;
    parse_string_array(
        conf, "correlation_contributions",
        &correlation_contributions,
        &correlation_contributions_len
    );

    par->correlation_contrib.den = 0;
    par->correlation_contrib.rsd = 0;
    par->correlation_contrib.d1 = 0;
    par->correlation_contrib.d2 = 0;
    par->correlation_contrib.g1 = 0;
    par->correlation_contrib.g2 = 0;
    par->correlation_contrib.g3 = 0;
    par->correlation_contrib.len = 0;
    par->correlation_contrib.g4 = 0;
    par->correlation_contrib.g5 = 0;

    for (int i = 0; i<correlation_contributions_len; ++i){
        if (strcmp(correlation_contributions[i], "den") == 0)
            par->correlation_contrib.den = 1;
        else if (strcmp(correlation_contributions[i], "rsd") == 0)
            par->correlation_contrib.rsd = 1;
        else if (strcmp(correlation_contributions[i], "len") == 0)
            par->correlation_contrib.len = 1;
        else if (strcmp(correlation_contributions[i], "d1") == 0)
            par->correlation_contrib.d1 = 1;
        else if (strcmp(correlation_contributions[i], "d2") == 0)
            par->correlation_contrib.d1 = 1;
        else if (strcmp(correlation_contributions[i], "g1") == 0)
            par->correlation_contrib.g1 = 1;
        else if (strcmp(correlation_contributions[i], "g2") == 0)
            par->correlation_contrib.g2 = 1;
        else if (strcmp(correlation_contributions[i], "g3") == 0)
            par->correlation_contrib.g3 = 1;
        else if (strcmp(correlation_contributions[i], "g4") == 0)
            par->correlation_contrib.g4 = 1;
        else if (strcmp(correlation_contributions[i], "g5") == 0)
            par->correlation_contrib.g5 = 1;
    }

    for (int i = 0; i < correlation_contributions_len; ++i)
        free(correlation_contributions[i]);
    free(correlation_contributions);

    par->divergent = 0;

    par->nonzero_terms[0].n = 0, par->nonzero_terms[0].l = 0;
    par->nonzero_terms[1].n = 0, par->nonzero_terms[1].l = 2;
    par->nonzero_terms[2].n = 0, par->nonzero_terms[2].l = 4;
    par->nonzero_terms[3].n = 1, par->nonzero_terms[3].l = 1;
    par->nonzero_terms[4].n = 1, par->nonzero_terms[4].l = 3;
    par->nonzero_terms[5].n = 2, par->nonzero_terms[5].l = 0;
    par->nonzero_terms[6].n = 2, par->nonzero_terms[6].l = 2;
    par->nonzero_terms[7].n = 3, par->nonzero_terms[7].l = 1;

    /* isolating the term requiring renormalization */
    if (
        par->correlation_contrib.d2 ||
        par->correlation_contrib.g1 ||
        par->correlation_contrib.g2 ||
        par->correlation_contrib.g3 ||
        par->correlation_contrib.g4 ||
        par->correlation_contrib.g5
    ){
        par->nonzero_terms[8].n = 4, par->nonzero_terms[8].l = 0;
        par->divergent = 1;
    }

    /* the output path */
    parse_string(conf, "output_path", par->output_path, COFFE_TRUE);

    /* the prefix for the output files */
    parse_string(conf, "output_prefix", par->output_prefix, COFFE_TRUE);

    /* flatsky parameter */
    parse_int(conf, "flatsky", &par->flatsky, COFFE_TRUE);

    if (par->flatsky){
        par->nonzero_terms[9].n = -10, par->nonzero_terms[9].l = -10;
    }

    /* parsing the k range */
    parse_double(conf, "k_min", &par->k_min, COFFE_TRUE);
    parse_double(conf, "k_max", &par->k_max, COFFE_TRUE);

    /* parsing the power spectrum */
#ifdef HAVE_CLASS

    parse_int(conf, "have_class", &par->have_class, COFFE_TRUE);
    if (par->have_class){
        parse_double(conf, "h", &par->h, COFFE_TRUE);
        parse_double(conf, "ln_10_pow_10_A_s", &par->ln_10_pow_10_A_s, COFFE_TRUE);
        parse_double(conf, "n_s", &par->n_s, COFFE_TRUE);
        parse_double(conf, "k_pivot", &par->k_pivot, COFFE_TRUE);
        parse_external_power_spectrum(par);
    }
    else{
#endif
        /* the power spectrum */
        parse_string(
            conf,
            "input_power_spectrum",
            par->file_power_spectrum,
            COFFE_TRUE
        );
        double *k, *pk;
        size_t pk_len;

        read_2col(
            par->file_power_spectrum,
            &k,
            &pk,
            &pk_len
        );

        coffe_init_spline(&par->power_spectrum, k, pk, pk_len, par->interp_method);

        free(k);
        free(pk);
#ifdef HAVE_CLASS
    }
#endif

    if (par->k_min < par->power_spectrum.spline->x[0]){
        fprintf(
            stderr,
            "WARNING: k_min smaller than minimum value in the file; "
            "using minimum value instead\n"
        );
        par->k_min = par->power_spectrum.spline->x[0];
    }

    if (par->k_max > par->power_spectrum.spline->x[
            par->power_spectrum.spline->size - 1
        ]){
        fprintf(
            stderr,
            "WARNING: k_max larger than maximum value in the file; "
            "using maximum value instead\n"
        );
        par->k_max =
            par->power_spectrum.spline->x[
                par->power_spectrum.spline->size - 1
            ];
    }


    /* normalizing the power spectrum to be dimensionless in our units */
    {
        size_t len = par->power_spectrum.spline->size;
        double *k_norm =
            (double *)coffe_malloc(sizeof(double)*len);
        double *pk_norm =
            (double *)coffe_malloc(sizeof(double)*len);
        for (size_t i = 0; i<len; ++i){
            k_norm[i] = par->power_spectrum.spline->x[i]/COFFE_H0;
            pk_norm[i] = par->power_spectrum.spline->y[i]*pow(COFFE_H0, 3);
        }
        coffe_init_spline(
            &par->power_spectrum_norm,
            k_norm, pk_norm, len,
            par->interp_method
        );
        par->k_min_norm = par->k_min/COFFE_H0;
        par->k_max_norm = par->k_max/COFFE_H0;
        free(k_norm);
        free(pk_norm);
    }


    /* saving the timestamp */
    snprintf(
        par->timestamp,
        COFFE_MAX_STRLEN,
        "%s",
        coffe_get_time()
    );

    par->conf = conf;

    end = clock();

    if (par->verbose)
        printf(
            "Settings file \"%s\" parsed in %.2f s\n",
            filename,
            (double)(end - start) / CLOCKS_PER_SEC
        );
    return EXIT_SUCCESS;
}
