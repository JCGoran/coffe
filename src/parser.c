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
    size_t *values_len
)
{
    config_setting_t *type = config_lookup(conf, setting);
    if (type == NULL){
        print_error(PROG_NULL_PARAM);
        return EXIT_FAILURE;
    }
    *values_len = (size_t)config_setting_length(type);

    *values = (char **)malloc(sizeof(char *) * *values_len);
    if (*values == NULL){
        print_error(PROG_ALLOC_ERROR);
        return EXIT_FAILURE;
    }

    const char *temp_setting;

    for (size_t i = 0; i < *values_len; ++i){
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
    size_t *values_len
)
{
    config_setting_t *type = config_lookup(conf, setting);
    if (type == NULL){
        print_error(PROG_NULL_PARAM);
        return EXIT_FAILURE;
    }
    *values_len = (size_t)config_setting_length(type);

    *values = (int *)coffe_malloc(sizeof(int) * *values_len);

    int temp_value;

    for (size_t i = 0; i < *values_len; ++i){
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
    size_t *values_len
)
{
    config_setting_t *type = config_lookup(conf, setting);
    if (type == NULL){
        print_error(PROG_NULL_PARAM);
        return EXIT_FAILURE;
    }
    *values_len = (size_t)config_setting_length(type);

    *values = (double *)coffe_malloc(sizeof(double) * *values_len);

    double temp_value;

    for (size_t i = 0; i < *values_len; ++i){
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
        COFFE_ARRAY_SIZE(redshifts),
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
    struct precision *ppr = coffe_malloc(sizeof(struct precision));
    struct background *pba = coffe_malloc(sizeof(struct background));
    struct thermo *pth = coffe_malloc(sizeof(struct thermo));
    struct perturbs *ppt = coffe_malloc(sizeof(struct perturbs));
    struct primordial *ppm = coffe_malloc(sizeof(struct primordial));
    struct spectra *psp = coffe_malloc(sizeof(struct spectra));
    struct nonlinear *pnl = coffe_malloc(sizeof(struct nonlinear));
    struct output *pop = coffe_malloc(sizeof(struct output));
    struct transfers *ptr = coffe_malloc(sizeof(struct transfers));
    struct lensing *ple = coffe_malloc(sizeof(struct lensing));
    ErrorMsg errmsg;

    clock_t class_start, class_end;

    if (par->verbose)
        printf("Launching CLASS...\n");

    class_start = clock();
    struct file_content *fc = coffe_malloc(sizeof(struct file_content));

    size_t class_parameters_len = 24, counter = 0;

    parser_init(fc, class_parameters_len, "", errmsg);

    /*
        not sure which values are actually necessary
        (maybe give the user the ability to read all of them?)
    */

    sprintf(fc->name[counter], "h");
    sprintf(fc->value[counter], "%e", par->h);
    ++counter;

    sprintf(fc->name[counter], "T_cmb");
    sprintf(fc->value[counter], "%e", par->T_cmb);
    ++counter;

    sprintf(fc->name[counter], "Omega_b");
    sprintf(fc->value[counter], "%e", par->Omega0_baryon);
    ++counter;

    sprintf(fc->name[counter], "N_ur");
    sprintf(fc->value[counter], "%e", par->N_ur);
    ++counter;

    sprintf(fc->name[counter], "Omega_cdm");
    sprintf(fc->value[counter], "%e", par->Omega0_cdm);
    ++counter;

    sprintf(fc->name[counter], "Omega_k");
    sprintf(fc->value[counter], "%e", 0.0);
    ++counter;

    sprintf(fc->name[counter], "w0_fld");
    sprintf(fc->value[counter], "%e", par->w0);
    ++counter;

    sprintf(fc->name[counter], "wa_fld");
    sprintf(fc->value[counter], "%e", par->wa);
    ++counter;

    sprintf(fc->name[counter], "output");
    sprintf(fc->value[counter], "mPk");
    ++counter;

    sprintf(fc->name[counter], "gauge");
    sprintf(fc->value[counter], "synchronous");
    ++counter;

    sprintf(fc->name[counter], "P_k_ini type");
    sprintf(fc->value[counter], "analytic_Pk");
    ++counter;

    sprintf(fc->name[counter], "N_ncdm");
    sprintf(fc->value[counter], "%d", par->N_ncdm);
    ++counter;

    sprintf(fc->name[counter], "m_ncdm");
    sprintf(fc->value[counter], "%e", par->m_ncdm);
    ++counter;

    sprintf(fc->name[counter], "k_pivot");
    sprintf(fc->value[counter], "%e", par->k_pivot);
    ++counter;

    sprintf(fc->name[counter], "n_s");
    sprintf(fc->value[counter], "%e", par->n_s);
    ++counter;

    sprintf(fc->name[counter], "sigma8");
    sprintf(fc->value[counter], "%e", par->sigma8);
    ++counter;

    sprintf(fc->name[counter], "alpha_s");
    sprintf(fc->value[counter], "%e", 0.0);
    ++counter;

    sprintf(fc->name[counter], "k_min_tau0");
    sprintf(fc->value[counter], "%e", 0.002);
    ++counter;

    sprintf(fc->name[counter], "P_k_max_h/Mpc");
    sprintf(fc->value[counter], "%e", par->k_max);
    ++counter;

    sprintf(fc->name[counter], "tau_reio");
    sprintf(fc->value[counter], "%e", 0.0925);
    ++counter;

    sprintf(fc->name[counter], "z_pk");
    sprintf(fc->value[counter], "%e", 0.0);
    ++counter;

    sprintf(fc->name[counter], "z_max_pk");
    sprintf(fc->value[counter], "%e", 3.0);
    ++counter;

    sprintf(fc->name[counter], "non linear");

    if (par->pk_type == 2)
        sprintf(fc->value[counter], "%s", "halofit");
    else if (par->pk_type == 3)
        sprintf(fc->value[counter], "%s", "hmcode");
    ++counter;

    if (
        input_init(
            fc,
            ppr,
            pba,
            pth,
            ppt,
            ptr,
            ppm,
            psp,
            pnl,
            ple,
            pop,
            errmsg
        ) == _FAILURE_
    ){
        fprintf(stderr, "\n\nError running input_init\n=>%s\n", errmsg);
        exit(EXIT_FAILURE);
    }

    /* the main CLASS sequence */

    background_init(ppr, pba);
    thermodynamics_init(ppr, pba, pth);
    perturb_init(ppr, pba, pth, ppt);
    primordial_init(ppr, ppt, ppm);
    nonlinear_init(ppr, pba, pth, ppt, ppm, pnl);
    transfer_init(ppr, pba, pth, ppt, pnl, ptr);
    spectra_init(ppr, pba, ppt, ppm, pnl, ptr, psp);

    class_end = clock();
    if (par->verbose)
        printf(
            "CLASS finished in %.2f s\n",
            (double)(class_end - class_start) / CLOCKS_PER_SEC
        );

    double *k = (double *)coffe_malloc(sizeof(double) * pnl->k_size);
    double *pk = (double *)coffe_malloc(sizeof(double) * pnl->k_size);
    size_t pk_len = pnl->k_size;

    nonlinear_pk_at_z(
        pba,
        pnl,
        logarithmic,
        pk_linear,
        0,
        pnl->index_pk_total,
        pk,
        NULL
    );

    /* rescaling as CLASS internally uses units of 1/Mpc */
    for (size_t i = 0; i < pk_len; ++i){
        k[i] = pnl->k[i] / par->h;
        pk[i] = exp(pk[i]) * pow(par->h, 3);
    }
    if (par->have_window){
        for (size_t i = 0; i < pk_len; ++i)
            pk[i] *= pow(
                coffe_resolution_window(par->window_size * k[i]),
                2
            );
    }

    double *k_norm = (double *)coffe_malloc(sizeof(double) * pnl->k_size);
    double *pk_norm = (double *)coffe_malloc(sizeof(double) * pnl->k_size);

    /* rescaling to the normalized one */
    for (size_t i = 0; i < pk_len; ++i){
        k_norm[i] = k[i] / COFFE_H0;
        pk_norm[i] = pk[i] * pow(COFFE_H0, 3);
    }

    coffe_init_spline(&par->power_spectrum, k, pk, pk_len, par->interp_method);
    coffe_init_spline(&par->power_spectrum_norm, k_norm, pk_norm, pk_len, par->interp_method);

    free(k);
    free(pk);
    free(k_norm);
    free(pk_norm);

    /* keeping it for later reuse */

    par->class_file_content = fc;
    par->class_background = pba;
    par->class_thermodynamics = pth;
    par->class_perturb = ppt;
    par->class_primordial = ppm;
    par->class_nonlinear = pnl;
    par->class_transfer = ptr;
    par->class_spectra = psp;

    return EXIT_SUCCESS;
}
#endif

int coffe_parse_default_parameters(
    struct coffe_parameters_t *par
)
{
    par->nthreads = 1;
    /* cosmological parameters */
    par->Omega0_m = 0.3;
    par->Omega0_baryon = 0.05;
    par->Omega0_cdm = par->Omega0_m - par->Omega0_baryon;
    par->Omega0_gamma = 9e-5;
    par->Omega0_de = 1 - (par->Omega0_cdm + par->Omega0_baryon + par->Omega0_gamma);
    par->w0 = -1.0;
    par->wa = 0.0;
    par->have_class = 0;
    par->N_ur = 2.0328;
    par->T_cmb = 2.726;
    par->N_ncdm = 1;
    par->m_ncdm = 0.06;
    par->h = 0.67;
    par->k_pivot = 0.05;
    par->sigma8 = 0.8156;
    par->n_s = 0.96;
    par->b_derivative = 0;
    par->f_derivative = 0;
    par->b_tilde_derivative = 0;
    par->f_tilde_derivative = 0;

    par->only_cross_correlations = 0;

    par->file_sep[0] = 0;
    const double separations[] = {10., 20., 40., 100., 120., 150., 200., 250., 300., 350.};
    par->sep = coffe_malloc(
        sizeof(double) * COFFE_ARRAY_SIZE(separations)
    );
    for (size_t i = 0; i < COFFE_ARRAY_SIZE(separations); ++i)
        par->sep[i] = separations[i];
    par->sep_len = COFFE_ARRAY_SIZE(separations);
    par->interp_method = 5;
    par->covariance_integration_method = 1;
    par->covariance_integration_bins = 8000;
    par->covariance_interpolation_method = 2;
    par->file_power_spectrum[0] = 0;
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
        k, pk, COFFE_ARRAY_SIZE(k),
        par->interp_method
    );

    par->k_min = 1e-5;
    par->k_max = 300.;
    par->pk_type = 0;
    par->zeldovich_approximation = 0;
    par->midpoint_approximation = 0;

    {
        size_t len = par->power_spectrum.spline->size;
        double *k_norm =
            (double *)coffe_malloc(sizeof(double)*len);
        double *pk_norm =
            (double *)coffe_malloc(sizeof(double)*len);
        for (size_t i = 0; i < len; ++i){
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

    par->galaxy_bias_analytic = 0;
    parse_bias_default(
        1.0, &par->galaxy_bias1, par->interp_method
    );
    par->read_galaxy_bias1 = 0;
    par->file_galaxy_bias1[0] = 0;

    parse_bias_default(
        1.0, &par->galaxy_bias2, par->interp_method
    );
    par->read_galaxy_bias2 = 0;
    par->file_galaxy_bias2[0] = 0;

    parse_bias_default(
        0.0, &par->magnification_bias1, par->interp_method
    );
    par->read_magnification_bias1 = 0;
    par->file_magnification_bias1[0] = 0;

    parse_bias_default(
        0.0, &par->magnification_bias2, par->interp_method
    );
    par->read_magnification_bias2 = 0;
    par->file_magnification_bias2[0] = 0;

    parse_bias_default(
        0.0, &par->evolution_bias1, par->interp_method
    );
    par->read_evolution_bias1 = 0;
    par->file_evolution_bias1[0] = 0;

    parse_bias_default(
        0.0, &par->evolution_bias2, par->interp_method
    );
    par->read_evolution_bias2 = 0;
    par->file_evolution_bias2[0] = 0;

    par->output_type = 2;
    par->covariance_density = NULL;
    par->covariance_density_len = 0;
    par->covariance_z_mean = NULL;
    par->covariance_z_mean_len = 0;
    par->covariance_deltaz = NULL;
    par->covariance_deltaz_len = 0;
    par->covariance_fsky = NULL;
    par->covariance_fsky_len = 0;
    par->covariance_pixelsize = NULL;
    par->covariance_pixelsize_len = 0;
    par->covariance_step_size = 10.0;
    par->covariance_zmin = NULL;
    par->covariance_zmin_len = 0;
    par->covariance_zmax = NULL;
    par->covariance_zmax_len = 0;
    par->covariance_minimum_separation = 0.0;

    par->have_window = 0;
    par->window_size = 0;

    snprintf(par->output_path, COFFE_MAX_STRLEN, "./");
    snprintf(par->output_prefix, COFFE_MAX_STRLEN, "$TIME");
    {
        char *temp_time = coffe_get_time();
        snprintf(
            par->timestamp,
            COFFE_MAX_STRLEN,
            "%s",
            temp_time
        );
        free(temp_time);
    }

    par->correlation_contrib.den = 1;
    par->correlation_contrib.rsd = 0;
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

    par->type_bg = NULL;
    par->type_bg_len = 0;

    par->background_bins = 10000;
    par->bessel_bins = 10000;
    #ifndef HAVE_CUBA
    par->integration_method = 2;
    #endif
    par->integration_bins = 750000;

    par->z_mean = 1.0;
    par->deltaz = 0.2;

    par->z_min = 0.9;
    par->z_max = 1.1;

    const int multipoles[] = {0, 2, 4};
    par->multipole_values = coffe_malloc(
        sizeof(double) * COFFE_ARRAY_SIZE(multipoles)
    );
    for (size_t i = 0; i < COFFE_ARRAY_SIZE(multipoles); ++i)
        par->multipole_values[i] = multipoles[i];
    par->multipole_values_len = COFFE_ARRAY_SIZE(multipoles);

    par->flatsky_local = 0;

    par->flatsky_local_nonlocal = 0;

    par->flatsky_nonlocal = 0;

    par->theta_len = 0;

    par->verbose = 0;

    par->conf = NULL;

    par->flag = 1;

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
    parse_double(conf, "omega_m", &par->Omega0_m, COFFE_TRUE);
    parse_double(conf, "omega_baryon", &par->Omega0_baryon, COFFE_TRUE);
    par->Omega0_cdm = par->Omega0_m - par->Omega0_baryon;
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
        const int error = read_1col(
            par->file_sep,
            &par->sep,
            &par->sep_len
        );
        /* hack: use default separations if the file is not available */
        if (error == EXIT_FAILURE){
            fprintf(
                stderr,
                "WARNING: cannot process file %s; using default separations instead!\n",
                par->file_sep
            );
        }
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

    /* integration method for covariance */
    parse_int(conf, "covariance_integration_method", &par->covariance_integration_method, COFFE_TRUE);
    if (
        par->covariance_integration_method != 1 &&
        par->covariance_integration_method != 2
    ){
        print_error_verbose(PROG_VALUE_ERROR, "covariance_integration_method");
        exit(EXIT_FAILURE);
    }

    if (par->covariance_integration_method == 2){
        /* integration method for covariance (only for FFT log) */
        parse_int(conf, "covariance_integration_bins", &par->covariance_integration_bins, COFFE_TRUE);
        if (par->covariance_integration_bins <= 0){
            print_error_verbose(PROG_VALUE_ERROR, "covariance_integration_bins");
            exit(EXIT_FAILURE);
        }
        parse_int(conf, "covariance_interpolation_method", &par->covariance_interpolation_method, COFFE_TRUE);
        if (
            par->covariance_interpolation_method != 1 &&
            par->covariance_interpolation_method != 2
        ){
            print_error_verbose(PROG_VALUE_ERROR, "covariance_interpolation_method");
            exit(EXIT_FAILURE);
        }
    }

    /* parsing the w parameter */
    parse_double(conf, "w0", &par->w0, COFFE_TRUE);
    parse_double(conf, "wa", &par->wa, COFFE_TRUE);

    /* the analytic galaxy bias has been disabled for now */
    //parse_int(conf, "galaxy_bias_analytic", &par->galaxy_bias_analytic, COFFE_FALSE);

    /* parsing the galaxy bias */
    if (!par->galaxy_bias_analytic){
        parse_int(conf, "read_galaxy_bias1", &par->read_galaxy_bias1, COFFE_FALSE);
        parse_bias(
            conf,
            "input_galaxy_bias1",
            par->file_galaxy_bias1,
            "galaxy_bias1",
            &par->galaxy_bias1,
            par->interp_method,
            par->read_galaxy_bias1
        );

        parse_int(conf, "read_galaxy_bias2", &par->read_galaxy_bias2, COFFE_FALSE);
        parse_bias(
            conf,
            "input_galaxy_bias2",
            par->file_galaxy_bias2,
            "galaxy_bias2",
            &par->galaxy_bias2,
            par->interp_method,
            par->read_galaxy_bias2
        );
    }
    else{
        /* need to hack this up */
        const size_t bins = 1024;
        const double redshift_max = 10.;
        double *x = (double *)coffe_malloc(sizeof(double) * bins);
        double *y = (double *)coffe_malloc(sizeof(double) * bins);
        for (size_t i = 0; i < bins; ++i){
            x[i] = (double)i * redshift_max / bins;
            y[i] = coffe_galaxy_bias(x[i]);
        }
        coffe_init_spline(
            &par->galaxy_bias1,
            x, y, bins,
            par->interp_method
        );

        coffe_init_spline(
            &par->galaxy_bias2,
            x, y, bins,
            par->interp_method
        );

        free(x);
        free(y);
    }

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
        parse_double_array(
            conf,
            "covariance_pixelsize",
            &par->covariance_pixelsize,
            &par->covariance_pixelsize_len
        );
        parse_double(
            conf,
            "covariance_step_size",
            &par->covariance_step_size,
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
            par->covariance_density_len != par->covariance_deltaz_len ||
            par->covariance_density_len != par->covariance_pixelsize_len
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

    /* parsing the contributions to the correlation function and covariance */
    char **correlation_contributions = NULL;
    size_t correlation_contributions_len = 0;
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

    for (size_t i = 0; i < correlation_contributions_len; ++i){
        if (strcmp(correlation_contributions[i], "den") == 0)
            par->correlation_contrib.den = 1;
        else if (strcmp(correlation_contributions[i], "rsd") == 0)
            par->correlation_contrib.rsd = 1;
        else if (strcmp(correlation_contributions[i], "len") == 0)
            par->correlation_contrib.len = 1;
        else if (strcmp(correlation_contributions[i], "d1") == 0)
            par->correlation_contrib.d1 = 1;
        else if (strcmp(correlation_contributions[i], "d2") == 0)
            par->correlation_contrib.d2 = 1;
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

    for (size_t i = 0; i < correlation_contributions_len; ++i)
        free(correlation_contributions[i]);
    free(correlation_contributions);

    /* isolating the term requiring renormalization */
    if (
        par->correlation_contrib.d2 ||
        par->correlation_contrib.g1 ||
        par->correlation_contrib.g2 ||
        par->correlation_contrib.g3 ||
        par->correlation_contrib.g4 ||
        par->correlation_contrib.g5
    ){
        par->divergent = 1;
    }

    /* the output path */
    parse_string(conf, "output_path", par->output_path, COFFE_TRUE);

    /* the prefix for the output files */
    parse_string(conf, "output_prefix", par->output_prefix, COFFE_TRUE);

    /* flatsky parameters */
    parse_int(conf, "flatsky_local", &par->flatsky_local, COFFE_TRUE);
    parse_int(conf, "flatsky_local_nonlocal", &par->flatsky_local_nonlocal, COFFE_TRUE);
    parse_int(conf, "flatsky_nonlocal", &par->flatsky_nonlocal, COFFE_TRUE);

    /* parameters for the Fisher forcat */
    parse_int(conf, "b_derivative", &par->b_derivative, COFFE_TRUE);
    parse_int(conf, "f_derivative", &par->f_derivative, COFFE_TRUE);
    parse_int(conf, "b_tilde_derivative", &par->b_tilde_derivative, COFFE_TRUE);
    parse_int(conf, "f_tilde_derivative", &par->f_tilde_derivative, COFFE_TRUE);

    /* cant have two or more at the same time */
    if (
        par->b_derivative
        +
        par->f_derivative
        +
        par->b_tilde_derivative
        +
        par->f_tilde_derivative > 1
    ){
        print_error_verbose(PROG_VALUE_ERROR, "derivatives");
        exit(EXIT_FAILURE);
    }

    /* parsing the k range */
    parse_double(conf, "k_min", &par->k_min, COFFE_TRUE);
    parse_double(conf, "k_max", &par->k_max, COFFE_TRUE);

    /* parsing the window */
    parse_int(conf, "have_window", &par->have_window, COFFE_TRUE);

    /* parsing whether we want just cross-correlations or not */
    parse_int(conf, "only_cross_correlations", &par->only_cross_correlations, COFFE_TRUE);

    /* parsing the size of the window (in Mpc/h) */
    if (par->have_window){
        parse_double(conf, "window_size", &par->window_size, COFFE_TRUE);
        /* boundary checking */
        if (par->window_size <= 0.0 || par->window_size >= 1000.){
            print_error_verbose(PROG_VALUE_ERROR, "window_size");
            exit(EXIT_FAILURE);
        }
    }

    /* parsing the power spectrum */
#ifdef HAVE_CLASS

    parse_int(conf, "have_class", &par->have_class, COFFE_TRUE);
    if (par->have_class){
        parse_int(conf, "midpoint_approximation", &par->midpoint_approximation, COFFE_TRUE);
        parse_int(conf, "pk_type", &par->pk_type, COFFE_TRUE);
        if (par->pk_type){
            parse_int(conf, "zeldovich_approximation", &par->zeldovich_approximation, COFFE_TRUE);
        }
        parse_double(conf, "h", &par->h, COFFE_TRUE);
        parse_double(conf, "sigma8", &par->sigma8, COFFE_TRUE);
        parse_double(conf, "n_s", &par->n_s, COFFE_TRUE);
        parse_double(conf, "k_pivot", &par->k_pivot, COFFE_TRUE);
        // new stuff for forecast
        parse_double(conf, "T_cmb", &par->T_cmb, COFFE_FALSE);
        parse_double(conf, "N_ur", &par->N_ur, COFFE_FALSE);
        parse_double(conf, "m_ncdm", &par->m_ncdm, COFFE_FALSE);
        parse_int(conf, "N_ncdm", &par->N_ncdm, COFFE_FALSE);
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

        const int error = read_2col(
            par->file_power_spectrum,
            &k,
            &pk,
            &pk_len
        );

        if (error == EXIT_SUCCESS){
            if (par->have_window){
                for (size_t i = 0; i < pk_len; ++i)
                    pk[i] *= pow(
                        coffe_resolution_window(par->window_size * k[i]),
                        2
                    );
            }

            coffe_init_spline(&par->power_spectrum, k, pk, pk_len, par->interp_method);

            free(k);
            free(pk);
        }
        /* hack: use default power spectrum */
        else{
            fprintf(
                stderr,
                "WARNING: cannot process file %s; using default P(k) instead!\n",
                par->file_power_spectrum
            );
#ifdef HAVE_CLASS
            fprintf(
                stderr,
                "Hint: consider setting the value of `have_class` to 1 to generate P(k) on the fly\n"
            );
#endif
        }
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
        for (size_t i = 0; i < len; ++i){
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
    {
        char *temp_time = coffe_get_time();
        snprintf(
            par->timestamp,
            COFFE_MAX_STRLEN,
            "%s",
            temp_time
        );
        free(temp_time);
    }

    par->conf = conf;

    end = clock();

    if (par->verbose)
        printf(
            "Settings file \"%s\" parsed in %.2f s\n",
            filename,
            (double)(end - start) / CLOCKS_PER_SEC
        );

    par->flag = 1;

    return EXIT_SUCCESS;
}
