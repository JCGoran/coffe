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


/**
    parses all the settings from the input file
    (given by argv[1]) into the structure <par>
**/

int coffe_parser_init(
    char *filename,
    struct coffe_parameters_t *par
)
{
    config_t *conf = (config_t *)coffe_malloc(sizeof(config_t));

    config_init(conf);

    if(!config_read_file(conf, filename))
    {
        fprintf(stderr, "ERROR: Cannot read the parameter file.\n");
        exit(EXIT_FAILURE);
    }

    clock_t start, end;

    printf("Parsing settings file \"%s\"...\n", filename);

    start = clock();

    config_set_auto_convert(conf, CONFIG_TRUE);

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
    if (par->read_matter_bias1 == COFFE_TRUE){
        double *bias_z, *bias_value;
        size_t bias_len;
        parse_string(conf, "input_matter_bias1", par->file_matter_bias1, COFFE_TRUE);
        read_2col(
            par->file_matter_bias1,
            &bias_z,
            &bias_value,
            &bias_len
        );
        init_spline(&par->matter_bias1, bias_z, bias_value, bias_len, par->interp_method);
        free(bias_z);
        free(bias_value);
    }
    else{
        double bias;
        parse_double(conf, "matter_bias1", &bias, COFFE_TRUE);
        /* a hacky way to init; if you need more range, increase redshift */
        double bias_redshift[] = {0, 25, 50, 75, 100};
        double bias_value[] = {bias, bias, bias, bias, bias};
        init_spline(&par->matter_bias1, bias_redshift, bias_value, sizeof(bias_redshift)/sizeof(bias_redshift[0]), par->interp_method);
    }

    parse_int(conf, "read_matter_bias2", &par->read_matter_bias2, COFFE_FALSE);
    if (par->read_matter_bias2 == COFFE_TRUE){
        double *bias_z, *bias_value;
        size_t bias_len;
        parse_string(conf, "input_matter_bias2", par->file_matter_bias2, COFFE_TRUE);
        read_2col(
            par->file_matter_bias2,
            &bias_z,
            &bias_value,
            &bias_len
        );
        init_spline(&par->matter_bias2, bias_z, bias_value, bias_len, par->interp_method);
        free(bias_z);
        free(bias_value);
    }
    else{
        double bias;
        parse_double(conf, "matter_bias2", &bias, COFFE_TRUE);
        /* a hacky way to init; if you need more range, increase redshift */
        double bias_redshift[] = {0, 25, 50, 75, 100};
        double bias_value[] = {bias, bias, bias, bias, bias};
        init_spline(&par->matter_bias2, bias_redshift, bias_value, sizeof(bias_redshift)/sizeof(bias_redshift[0]), par->interp_method);
    }


    /* parsing the magnification bias (s) */
    parse_int(conf, "read_magnification_bias1", &par->read_magnification_bias1, COFFE_FALSE);
    if (par->read_magnification_bias1 == COFFE_TRUE){
        double *s_z, *s_value;
        size_t s_len;
        parse_string(
            conf,
            "input_magnification_bias1",
            par->file_magnification_bias1,
            COFFE_TRUE
        );
        read_2col(
            par->file_magnification_bias1,
            &s_z,
            &s_value,
            &s_len
        );
        init_spline(
            &par->magnification_bias1,
            s_z, s_value, s_len,
            par->interp_method
        );
        free(s_z);
        free(s_value);
    }
    else{
        double s;
        parse_double(conf, "magnification_bias1", &s, COFFE_TRUE);
        /* a hacky way to init; if you need more range, increase redshift */
        double s_redshift[] = {0, 25, 50, 75, 100};
        double s_value[] = {s, s, s, s, s};
        init_spline(&par->magnification_bias1, s_redshift, s_value, sizeof(s_redshift)/sizeof(s_redshift[0]), par->interp_method);
    }

    parse_int(conf, "read_magnification_bias2", &par->read_magnification_bias2, COFFE_FALSE);
    if (par->read_magnification_bias2 == COFFE_TRUE){
        double *s_z, *s_value;
        size_t s_len;
        parse_string(conf, "input_magnification_bias2", par->file_magnification_bias2, COFFE_TRUE);
        read_2col(
            par->file_magnification_bias2,
            &s_z,
            &s_value,
            &s_len
        );
        init_spline(&par->magnification_bias2, s_z, s_value, s_len, par->interp_method);
        free(s_z);
        free(s_value);
    }
    else{
        double s;
        parse_double(conf, "magnification_bias2", &s, COFFE_TRUE);
        /* a hacky way to init; if you need more range, increase redshift */
        double s_redshift[] = {0, 25, 50, 75, 100};
        double s_value[] = {s, s, s, s, s};
        init_spline(&par->magnification_bias2, s_redshift, s_value, sizeof(s_redshift)/sizeof(s_redshift[0]), par->interp_method);
    }

    /* parsing the evolution_bias (f_evo) */
    parse_int(conf, "read_evolution_bias1", &par->read_evolution_bias1, COFFE_FALSE);
    if (par->read_evolution_bias1 == COFFE_TRUE){
        double *f_evo_z, *f_evo_value;
        size_t f_evo_len;
        parse_string(conf, "input_evolution_bias1", par->file_evolution_bias1, COFFE_TRUE);
        read_2col(
            par->file_evolution_bias1,
            &f_evo_z,
            &f_evo_value,
            &f_evo_len
        );
        init_spline(&par->evolution_bias1, f_evo_z, f_evo_value, f_evo_len, par->interp_method);
        free(f_evo_z);
        free(f_evo_value);
    }
    else{
        double f_evo;
        parse_double(conf, "evolution_bias1", &f_evo, COFFE_TRUE);
        /* a hacky way to init; if you need more range, increase redshift */
        double f_evo_redshift[] = {0, 25, 50, 75, 200};
        double f_evo_value[] = {f_evo, f_evo, f_evo, f_evo, f_evo};
        init_spline(&par->evolution_bias1, f_evo_redshift, f_evo_value, sizeof(f_evo_redshift)/sizeof(f_evo_redshift[0]), par->interp_method);
    }

    parse_int(conf, "read_evolution_bias2", &par->read_evolution_bias2, COFFE_FALSE);
    if (par->read_evolution_bias2 == COFFE_TRUE){
        double *f_evo_z, *f_evo_value;
        size_t f_evo_len;
        parse_string(conf, "input_evolution_bias2", par->file_evolution_bias2, COFFE_TRUE);
        read_2col(
            par->file_evolution_bias2,
            &f_evo_z,
            &f_evo_value,
            &f_evo_len
        );
        init_spline(&par->evolution_bias2, f_evo_z, f_evo_value, f_evo_len, par->interp_method);
        free(f_evo_z);
        free(f_evo_value);
    }
    else{
        double f_evo;
        parse_double(conf, "evolution_bias2", &f_evo, COFFE_TRUE);
        /* a hacky way to init; if you need more range, increase redshift */
        double f_evo_redshift[] = {0, 25, 50, 75, 200};
        double f_evo_value[] = {f_evo, f_evo, f_evo, f_evo, f_evo};
        init_spline(&par->evolution_bias2, f_evo_redshift, f_evo_value, sizeof(f_evo_redshift)/sizeof(f_evo_redshift[0]), par->interp_method);
    }

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

    /* parsing the power spectrum */
#ifdef HAVE_CLASS

    parse_int(conf, "have_class", &par->have_class, COFFE_TRUE);
    if (par->have_class){
        parse_double(conf, "h", &par->h, COFFE_TRUE);
        parse_double(conf, "ln_10_pow_10_A_s", &par->ln_10_pow_10_A_s, COFFE_TRUE);
        parse_double(conf, "n_s", &par->n_s, COFFE_TRUE);
        parse_double(conf, "k_pivot", &par->k_pivot, COFFE_TRUE);

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

        parse_double(conf, "k_min", &par->k_min, COFFE_TRUE);
        parse_double(conf, "k_max", &par->k_max, COFFE_TRUE);

        printf("Launching CLASS...\n");

        class_start = clock();
        struct file_content fc;

        parser_init(&fc, 18, "", errmsg);

        /* not sure which values I actually need (maybe give the user the ability to read all of them?) */

        sprintf(fc.name[0], "h");
        sprintf(fc.value[0], "%e", par->h);

        sprintf(fc.name[1], "T_cmb");
        sprintf(fc.value[1], "%e", 2.7255);

        sprintf(fc.name[2], "Omega_b");
        sprintf(fc.value[2], "%e", par->Omega0_baryon);

        sprintf(fc.name[3], "N_ur");
        sprintf(fc.value[3], "%e", 3.046);

        sprintf(fc.name[4], "Omega_cdm");
        sprintf(fc.value[4], "%e", par->Omega0_cdm);

        sprintf(fc.name[5], "Omega_k");
        sprintf(fc.value[5], "%e", 0.0);

        sprintf(fc.name[6], "w0_fld");
        sprintf(fc.value[6], "%e", par->w0);

        sprintf(fc.name[7], "wa_fld");
        sprintf(fc.value[7], "%e", par->wa);

        sprintf(fc.name[8], "output");
        sprintf(fc.value[8], "mPk");

        sprintf(fc.name[9], "gauge");
        sprintf(fc.value[9], "synchronous");

        sprintf(fc.name[10], "P_k_ini type");
        sprintf(fc.value[10], "analytic_Pk");

        sprintf(fc.name[11], "k_pivot");
        sprintf(fc.value[11], "%e", par->k_pivot);

        sprintf(fc.name[12], "ln10^{10}A_s");
        sprintf(fc.value[12], "%e", par->ln_10_pow_10_A_s);

        sprintf(fc.name[13], "n_s");
        sprintf(fc.value[13], "%e", par->n_s);

        sprintf(fc.name[14], "alpha_s");
        sprintf(fc.value[14], "%e", 0.0);

        sprintf(fc.name[15], "k_min_tau0");
        sprintf(fc.value[15], "%e", 0.002);

        sprintf(fc.name[16], "P_k_max_h/Mpc");
        sprintf(fc.value[16], "%e", par->k_max);

        sprintf(fc.name[17], "z_pk");
        sprintf(fc.value[17], "%e", 0.0);

        if (
            input_init(
                &fc, &ppr, &pba, &pth, &ppt, &ptr,
                &ppm, &psp, &pnl, &ple, &pop, errmsg
            ) == _FAILURE_
        ){
            printf("\n\nError running input_init\n=>%s\n", errmsg);
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
        printf(
            "CLASS finished in %.2f s\n",
            (double)(class_end - class_start) / CLOCKS_PER_SEC
        );

        double *k =(double *)coffe_malloc(sizeof(double)*psp.ln_k_size);
        double *pk =(double *)coffe_malloc(sizeof(double)*psp.ln_k_size);
        size_t pk_len = psp.ln_k_size;

        for (size_t i = 0; i<(size_t)psp.ln_k_size; ++i){
            k[i] = exp(psp.ln_k[i])/par->h;
            pk[i] = exp(psp.ln_pk_l[i])*pow(par->h, 3);
        }
        init_spline(&par->power_spectrum, k, pk, pk_len, par->interp_method);

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
    }
    else{
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

        init_spline(&par->power_spectrum, k, pk, pk_len, par->interp_method);

        free(k);
        free(pk);

        parse_double(conf, "k_min", &par->k_min, COFFE_FALSE);

        if (par->k_min < par->power_spectrum.spline->x[0]){
            fprintf(
                stderr,
                "WARNING: k_min smaller than minimum value in the file; "
                "using minimum value instead\n"
            );
            par->k_min = par->power_spectrum.spline->x[0];
        }

        parse_double(conf, "k_max", &par->k_max, COFFE_FALSE);

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
    }
#else

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
        &k, &pk, &pk_len
    );

    init_spline(&par->power_spectrum, k, pk, pk_len, par->interp_method);

    free(k);
    free(pk);

    /* the lower limit of integration for P(k) */
    parse_double(conf, "k_min", &par->k_min, COFFE_FALSE);

    if (par->k_min < par->power_spectrum.spline->x[0]){
        fprintf(
            stderr,
            "WARNING: k_min smaller than minimum value in the file; "
            "using minimum value instead\n"
        );
        par->k_min = par->power_spectrum.spline->x[0];
    }

    /* the upper limit of integration for P(k) */
    parse_double(conf, "k_max", &par->k_max, COFFE_FALSE);

    if (par->k_max > par->power_spectrum.spline->x[par->power_spectrum.spline->size - 1]){
        fprintf(
            stderr,
            "WARNING: k_max larger than maximum value in the file; "
            "using maximum value instead\n"
        );
        par->k_max =
            par->power_spectrum.spline->x[par->power_spectrum.spline->size - 1];
    }
#endif

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
        init_spline(
            &par->power_spectrum_norm,
            k_norm, pk_norm, len,
            par->interp_method
        );
        par->k_min_norm = par->k_min/COFFE_H0;
        par->k_max_norm = par->k_max/COFFE_H0;
        free(k_norm);
        free(pk_norm);
    }

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

    /* parsing the contributions to the correlation function */
    if (
        par->output_type == 0 ||
        par->output_type == 1 ||
        par->output_type == 2 ||
        par->output_type == 3 ||
        par->output_type == 6
    ){
        parse_string_array(
            conf, "correlation_contributions",
            &par->correlation_sources,
            &par->correlation_sources_len
        );

        char possible_inputs[10][10]
            = {"den", "rsd", "d1", "d2", "g1", "g2", "g3", "g4", "g5", "len"};

        int counter = 0;
        char temp_input[10];
        for (int i = 0; i<par->correlation_sources_len*(par->correlation_sources_len + 1)/2; ++i){
            strcpy(par->corr_terms[i], "\0");
        }

        /* input of the possible correlation terms */
        for (int i = 0; i<par->correlation_sources_len; ++i){
            for (int j = i; j<par->correlation_sources_len; ++j){
                for (int l = 0; l<10; ++l){
                    sprintf(temp_input, "%d", l);
                    if (strcmp(par->correlation_sources[i], possible_inputs[l]) == 0)
                        strcat(par->corr_terms[counter], temp_input);
                    if (strcmp(par->correlation_sources[j], possible_inputs[l]) == 0)
                        strcat(par->corr_terms[counter], temp_input);
                }
                ++counter;
            }
        }

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
        for (int i = 0; i<counter; ++i){
            if (
                strcmp(par->corr_terms[i], "33") == 0 || // d2-d2 term
                strcmp(par->corr_terms[i], "44") == 0 || // g1-g1 term
                strcmp(par->corr_terms[i], "55") == 0 || // g2-g2 term
                strcmp(par->corr_terms[i], "66") == 0 || // g3-g3 term
                strcmp(par->corr_terms[i], "77") == 0 || // g4-g4 term
                strcmp(par->corr_terms[i], "88") == 0 || // g5-g5 term
                /* I don't think these are necessary anymore */
                strcmp(par->corr_terms[i], "34") == 0 ||
                strcmp(par->corr_terms[i], "43") == 0 ||
                strcmp(par->corr_terms[i], "35") == 0 ||
                strcmp(par->corr_terms[i], "53") == 0 ||
                strcmp(par->corr_terms[i], "36") == 0 ||
                strcmp(par->corr_terms[i], "63") == 0 ||
                strcmp(par->corr_terms[i], "45") == 0 ||
                strcmp(par->corr_terms[i], "54") == 0 ||
                strcmp(par->corr_terms[i], "46") == 0 ||
                strcmp(par->corr_terms[i], "64") == 0 ||
                strcmp(par->corr_terms[i], "56") == 0 ||
                strcmp(par->corr_terms[i], "65") == 0
            ){
                par->nonzero_terms[8].n = 4, par->nonzero_terms[8].l = 0;
                par->divergent = 1;
            }
        }
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

    /* saving the timestamp */
    sprintf(par->timestamp, "%s", coffe_get_time());

    par->conf = conf;

    end = clock();

    printf(
        "Settings file \"%s\" parsed in %.2f s\n",
        filename,
        (double)(end - start) / CLOCKS_PER_SEC
    );
    return EXIT_SUCCESS;
}
