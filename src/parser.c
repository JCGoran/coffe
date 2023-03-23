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

#include "common.h"
#include "parser.h"
#include "errors.h"

#ifdef HAVE_CLASS
#include "class.h"
#endif



static int parse_bias_default(
    double value,
    coffe_interpolation *spline,
    const enum coffe_interp1d_type method
)
{
    double redshifts[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double values[] = {value, value, value, value, value, value, value, value, value, value, value};
    coffe_init_spline(
        spline, redshifts, values,
        COFFE_ARRAY_SIZE(redshifts),
        method
    );

    return EXIT_SUCCESS;
}



#ifndef COFFE_CYTHON

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


/**
    parses the value of bias given by <setting> from config <conf>
    into coffe_interpolation <spline> depending on value of <flag>
**/

static int parse_bias(
    config_t *conf,
    char *setting,
    char *filename,
    char *setting_single,
    coffe_interpolation *spline,
    const enum coffe_interp1d_type method,
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

#endif




/**
    parsing the power spectrum from an external source (only CLASS for now)
**/
#ifdef HAVE_CLASS
int parse_external_power_spectrum(
    coffe_parameters_t *par
)
{
    /* threshold for comparing with zero */
    const double EPSILON_ABS = 1e-10;

    /* so we don't leak memory */
    coffe_free_spline(&par->power_spectrum);
    coffe_free_spline(&par->power_spectrum_norm);

    struct precision *ppr = (struct precision *)coffe_malloc(sizeof(struct precision));
    struct background *pba = (struct background *)coffe_malloc(sizeof(struct background));
    struct thermo *pth = (struct thermo *)coffe_malloc(sizeof(struct thermo));
    struct perturbs *ppt = (struct perturbs *)coffe_malloc(sizeof(struct perturbs));
    struct primordial *ppm = (struct primordial *)coffe_malloc(sizeof(struct primordial));
    struct spectra *psp = (struct spectra *)coffe_malloc(sizeof(struct spectra));
    struct nonlinear *pnl = (struct nonlinear *)coffe_malloc(sizeof(struct nonlinear));
    struct output *pop = (struct output *)coffe_malloc(sizeof(struct output));
    struct transfers *ptr = (struct transfers *)coffe_malloc(sizeof(struct transfers));
    struct lensing *ple = (struct lensing *)coffe_malloc(sizeof(struct lensing));
    ErrorMsg errmsg;

    clock_t class_start, class_end;

    if (par->verbose)
        printf("Launching CLASS...\n");

    class_start = clock();
    struct file_content *fc = (struct file_content *)coffe_malloc(sizeof(struct file_content));

    size_t class_parameters_len = 35, counter = 0;

    parser_init(fc, class_parameters_len, "", errmsg);

    /*
        not sure which values are actually necessary
        (maybe give the user the ability to read all of them?)
    */

    sprintf(fc->name[counter], "h");
    sprintf(fc->value[counter], "%e", par->h);
    if (par->verbose)
        fprintf(stderr, "h = %e\n", par->h);
    ++counter;

    sprintf(fc->name[counter], "T_cmb");
    sprintf(fc->value[counter], "%e", par->T_cmb);
    if (par->verbose)
        fprintf(stderr, "T_cmb = %e\n", par->T_cmb);
    ++counter;

    sprintf(fc->name[counter], "Omega_b");
    sprintf(fc->value[counter], "%e", par->Omega0_baryon);
    if (par->verbose)
        fprintf(stderr, "Omega_b = %e\n", par->Omega0_baryon);
    ++counter;

    sprintf(fc->name[counter], "N_ur");
    sprintf(fc->value[counter], "%e", par->N_ur);
    if (par->verbose)
        fprintf(stderr, "N_ur = %e\n", par->N_ur);
    ++counter;

    /* set the helium fraction to a constant value (must be different from zero) */
    if (fabs(par->YHe) > EPSILON_ABS){
        sprintf(fc->name[counter], "YHe");
        sprintf(fc->value[counter], "%e", par->YHe);
        if (par->verbose)
            fprintf(stderr, "YHe = %e\n", par->YHe);
        ++counter;
    }

    sprintf(fc->name[counter], "Omega_cdm");
    sprintf(fc->value[counter], "%e", par->Omega0_cdm);
    if (par->verbose)
        fprintf(stderr, "Omega_cdm = %e\n", par->Omega0_cdm);
    ++counter;

    sprintf(fc->name[counter], "Omega_k");
    sprintf(fc->value[counter], "%e", 0.0);
    if (par->verbose)
        fprintf(stderr, "Omega_k = %e\n", 0.0);
    ++counter;

    sprintf(fc->name[counter], "w0_fld");
    sprintf(fc->value[counter], "%e", par->w0);
    if (par->verbose)
        fprintf(stderr, "w0_fld = %e\n", par->w0);
    ++counter;

    sprintf(fc->name[counter], "wa_fld");
    sprintf(fc->value[counter], "%e", par->wa);
    if (par->verbose)
        fprintf(stderr, "wa_fld = %e\n", par->wa);
    ++counter;

    sprintf(fc->name[counter], "output");
    sprintf(fc->value[counter], "mPk");
    if (par->verbose)
        fprintf(stderr, "output = mPk\n");
    ++counter;

    sprintf(fc->name[counter], "gauge");
    sprintf(fc->value[counter], "synchronous");
    if (par->verbose)
        fprintf(stderr, "gauge = synchronous\n");
    ++counter;

    sprintf(fc->name[counter], "P_k_ini type");
    sprintf(fc->value[counter], "analytic_Pk");
    if (par->verbose)
        fprintf(stderr, "P_k_ini type = analytic_Pk\n");
    ++counter;

    sprintf(fc->name[counter], "N_ncdm");
    sprintf(fc->value[counter], "%d", par->N_ncdm);
    if (par->verbose)
        fprintf(stderr, "N_ncdm = %d\n", par->N_ncdm);
    ++counter;

    sprintf(fc->name[counter], "m_ncdm");
    sprintf(fc->value[counter], "%e", par->m_ncdm);
    if (par->verbose)
        fprintf(stderr, "m_ncdm = %e\n", par->m_ncdm);
    ++counter;

    sprintf(fc->name[counter], "k_pivot");
    sprintf(fc->value[counter], "%e", par->k_pivot);
    if (par->verbose)
        fprintf(stderr, "k_pivot = %e\n", par->k_pivot);
    ++counter;

    sprintf(fc->name[counter], "n_s");
    sprintf(fc->value[counter], "%e", par->n_s);
    if (par->verbose)
        fprintf(stderr, "n_s = %e\n", par->n_s);
    ++counter;

    if (fabs(par->A_s) > EPSILON_ABS){
        sprintf(fc->name[counter], "A_s");
        sprintf(fc->value[counter], "%e", par->A_s);
        if (par->verbose)
            fprintf(stderr, "A_s = %e\n", par->A_s);
        ++counter;
    }
    else{
        sprintf(fc->name[counter], "sigma8");
        sprintf(fc->value[counter], "%e", par->sigma8);
        if (par->verbose)
            fprintf(stderr, "sigma8 = %e\n", par->sigma8);
        ++counter;
    }

    sprintf(fc->name[counter], "alpha_s");
    sprintf(fc->value[counter], "%e", 0.0);
    if (par->verbose)
        fprintf(stderr, "alpha_s = %e\n", 0.0);
    ++counter;

    sprintf(fc->name[counter], "k_min_tau0");
    sprintf(fc->value[counter], "%e", 0.002);
    if (par->verbose)
        fprintf(stderr, "k_min_tau0 = %e\n", 0.002);
    ++counter;

    sprintf(fc->name[counter], "P_k_max_h/Mpc");
    sprintf(fc->value[counter], "%e", par->k_max);
    if (par->verbose)
        fprintf(stderr, "P_k_max_h/Mpc = %e\n", par->k_max);
    ++counter;

    sprintf(fc->name[counter], "tau_reio");
    sprintf(fc->value[counter], "%e", 0.0925);
    if (par->verbose)
        fprintf(stderr, "tau_reio = %e\n", 0.0925);
    ++counter;

    sprintf(fc->name[counter], "k_per_decade_for_bao");
    sprintf(fc->value[counter], "%e", par->class_precision.k_per_decade_for_bao);
    if (par->verbose)
        fprintf(stderr, "k_per_decade_for_bao = %e\n", par->class_precision.k_per_decade_for_bao);
    ++counter;

    sprintf(fc->name[counter], "k_per_decade_for_pk");
    sprintf(fc->value[counter], "%e", par->class_precision.k_per_decade_for_pk);
    if (par->verbose)
        fprintf(stderr, "k_per_decade_for_pk = %e\n", par->class_precision.k_per_decade_for_pk);
    ++counter;

    sprintf(fc->name[counter], "start_large_k_at_tau_h_over_tau_k");
    sprintf(fc->value[counter], "%e", par->class_precision.start_large_k_at_tau_h_over_tau_k);
    if (par->verbose)
        fprintf(stderr, "start_large_k_at_tau_h_over_tau_k = %e\n", par->class_precision.start_large_k_at_tau_h_over_tau_k);
    ++counter;

    sprintf(fc->name[counter], "l_max_g");
    sprintf(fc->value[counter], "%d", par->class_precision.l_max_g);
    if (par->verbose)
        fprintf(stderr, "l_max_g = %d\n", par->class_precision.l_max_g);
    ++counter;

    sprintf(fc->name[counter], "l_max_ur");
    sprintf(fc->value[counter], "%d", par->class_precision.l_max_ur);
    if (par->verbose)
        fprintf(stderr, "l_max_ur = %d\n", par->class_precision.l_max_ur);
    ++counter;

    sprintf(fc->name[counter], "tol_perturb_integration");
    sprintf(fc->value[counter], "%e", par->class_precision.tol_perturb_integration);
    if (par->verbose)
        fprintf(stderr, "tol_perturb_integration = %e\n", par->class_precision.tol_perturb_integration);
    ++counter;

    sprintf(fc->name[counter], "radiation_streaming_trigger_tau_over_tau_k");
    sprintf(fc->value[counter], "%e", par->class_precision.radiation_streaming_trigger_tau_over_tau_k);
    if (par->verbose)
        fprintf(stderr, "radiation_streaming_trigger_tau_over_tau_k = %e\n", par->class_precision.radiation_streaming_trigger_tau_over_tau_k);
    ++counter;

    sprintf(fc->name[counter], "ur_fluid_trigger_tau_over_tau_k");
    sprintf(fc->value[counter], "%e", par->class_precision.ur_fluid_trigger_tau_over_tau_k);
    if (par->verbose)
        fprintf(stderr, "ur_fluid_trigger_tau_over_tau_k = %e\n", par->class_precision.ur_fluid_trigger_tau_over_tau_k);
    ++counter;

    sprintf(fc->name[counter], "z_pk");
    sprintf(fc->value[counter], "%e", 0.0);
    if (par->verbose)
        fprintf(stderr, "z_pk = %e\n", 0.0);
    ++counter;

    sprintf(fc->name[counter], "z_max_pk");
    sprintf(fc->value[counter], "%e", 3.0);
    if (par->verbose)
        fprintf(stderr, "z_max_pk = %e\n", 3.0);
    ++counter;

    sprintf(fc->name[counter], "non linear");

    if (par->pk_type == COFFE_PK_NONLINEAR_HALOFIT){
        sprintf(fc->value[counter], "%s", "halofit");
        if (par->verbose)
            fprintf(stderr, "non linear = %s\n", "halofit");
    }
    else if (par->pk_type == COFFE_PK_NONLINEAR_HMCODE){
        sprintf(fc->value[counter], "%s", "hmcode");
        if (par->verbose)
            fprintf(stderr, "non linear = %s\n", "hmcode");
    }
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

    for (size_t i = 0; i < pk_len; ++i){
        k[i] = pnl->k[i];
        pk[i] = exp(pk[i]) * exp(-pnl->k[i] * par->inv_k_window);
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
        k_norm[i] = k[i];
        pk_norm[i] = pk[i];
    }

    if (par->k_min < k[0]){
        par->k_min = k[0];
        par->k_min_norm = k_norm[0];
    }
    if (par->k_max > k[pk_len - 1]){
        par->k_max = k[pk_len - 1];
        par->k_max_norm = k_norm[pk_len - 1];
    }

    coffe_init_spline(&par->power_spectrum, k, pk, pk_len, par->interp_method);
    coffe_init_spline(&par->power_spectrum_norm, k_norm, pk_norm, pk_len, par->interp_method);

    free(k);
    free(pk);
    free(k_norm);
    free(pk_norm);

    /* keeping it for later reuse */

    par->class_struct.file_content = fc;
    par->class_struct.background = pba;
    par->class_struct.thermodynamics = pth;
    par->class_struct.perturb = ppt;
    par->class_struct.primordial = ppm;
    par->class_struct.nonlinear = pnl;
    par->class_struct.transfer = ptr;
    par->class_struct.spectra = psp;

    lensing_free(ple);
    free(ple);
    free(ppr);
    free(pop);

    if (par->pk_type != COFFE_PK_LINEAR){

        /* so we don't leak memory */
        coffe_free_spline2d(&par->power_spectrum2d);
        coffe_free_spline2d(&par->power_spectrum2d_norm);

        enum pk_outputs pk_type = pk_linear;
        if (
            par->pk_type == COFFE_PK_NONLINEAR_HALOFIT ||
            par->pk_type == COFFE_PK_NONLINEAR_HMCODE
        )
            pk_type = pk_nonlinear;

        /* TODO make all of this modular */
        const double z_min = 0.0;
        const double z_max = 3.0;

        /* list of all the redshifts */
        const size_t z_size = 100;
        double *z = (double *)coffe_malloc(sizeof(double) * z_size);
        for (size_t i = 0; i < z_size; ++i)
            z[i] = z_min + z_max * (double)i / z_size;

        /* size of the wavevector array k */
        const size_t k_size = ((struct nonlinear *)par->class_struct.nonlinear)->k_size;
        double *k = (double *)coffe_malloc(sizeof(double) * k_size);
        double *k_norm = (double *)coffe_malloc(sizeof(double) * k_size);
        /* alloc memory for P(k) that CLASS outputs */
        double *pk = (double *)coffe_malloc(sizeof(double) * k_size);

        for (size_t j = 0; j < k_size; ++j){
            k[j] = ((struct nonlinear *)par->class_struct.nonlinear)->k[j];
            k_norm[j] = ((struct nonlinear *)par->class_struct.nonlinear)->k[j];
        }

        /* alloc memory for 2D interpolation */
        double *pk_at_z2d = (double *)coffe_malloc(sizeof(double) * z_size * k_size);
        double *pk_at_z2d_norm = (double *)coffe_malloc(sizeof(double) * z_size * k_size);

        /* go over each redshift */
        for (size_t i = 0; i < z_size; ++i){
            /* get the (non)linear power spectrum at redshift z (and store it in pk) */
            nonlinear_pk_at_z(
                (struct background *)par->class_struct.background,
                (struct nonlinear *)par->class_struct.nonlinear,
                logarithmic,
                pk_type,
                z[i],
                ((struct nonlinear *)par->class_struct.nonlinear)->index_pk_total,
                pk,
                NULL
            );

            for (size_t j = 0; j < k_size; ++j){
                pk_at_z2d[j * z_size + i] = exp(pk[j]);
                pk_at_z2d_norm[j * z_size + i] = exp(pk[j]);
            }
        }

        coffe_init_spline2d(
            &par->power_spectrum2d,
            z,
            k,
            pk_at_z2d,
            z_size,
            k_size,
            COFFE_INTERP2D_BICUBIC
        );

        coffe_init_spline2d(
            &par->power_spectrum2d_norm,
            z,
            k,
            pk_at_z2d_norm,
            z_size,
            k_size,
            COFFE_INTERP2D_BICUBIC
        );

        free(k);
        free(k_norm);
        free(pk);
        free(pk_at_z2d);
        free(pk_at_z2d_norm);
        free(z);
    }

    return EXIT_SUCCESS;
}
#endif


int coffe_parse_default_parameters(
    coffe_parameters_t *par
)
{
    #ifdef _OPENMP
    par->nthreads = omp_get_max_threads();
    #else
    par->nthreads = 1;
    #endif
    /* cosmological parameters */
    par->Omega0_m = 0.3;
    par->Omega0_baryon = 0.05;
    par->Omega0_gamma = 9e-5;
    par->h = 0.67;
    par->w0 = -1.0;
    par->wa = 0.0;
    par->have_class = 0;
    par->N_ur = 2.0328;
    par->T_cmb = 2.726;
    par->N_ncdm = 1;
    par->m_ncdm = 0.00;
    par->YHe = 0; /* this indicates the BBN table should be used instead */
    /* see eq. (19) of https://arxiv.org/abs/1212.6154 */
    par->Omega0_nu = par->m_ncdm / 93.14 / par->h / par->h;
    par->Omega0_cdm = par->Omega0_m - par->Omega0_baryon - par->Omega0_nu;
    par->Omega0_de = 1 - (par->Omega0_m + par->Omega0_gamma);
    par->k_pivot = 0.05;
    par->sigma8 = 0.8156;
    par->A_s = 0;
    par->n_s = 0.96;
    par->b_derivative = 0;
    par->f_derivative = 0;
    par->b_tilde_derivative = 0;
    par->f_tilde_derivative = 0;
    par->inv_k_window = 0;
    par->integration_1d_type = COFFE_INTEGRATION_GSL;
    par->integration_1d_prec = 1e-5;

    coffe_new_fit_coefficients_array(&par->galaxy_bias1_coefficients);
    coffe_new_fit_coefficients_array(&par->galaxy_bias2_coefficients);
    coffe_new_fit_coefficients_array(&par->magnification_bias1_coefficients);
    coffe_new_fit_coefficients_array(&par->magnification_bias2_coefficients);
    coffe_new_fit_coefficients_array(&par->evolution_bias1_coefficients);
    coffe_new_fit_coefficients_array(&par->evolution_bias2_coefficients);

    #ifdef HAVE_CLASS
    par->has_class = 1;
    #else
    par->has_class = 0;
    #endif

    #ifdef HAVE_CUBA
    par->has_cuba = 1;
    #else
    par->has_cuba = 0;
    #endif

    coffe_new_class_struct(&par->class_struct);

    par->class_precision.k_per_decade_for_bao = 70.;
    par->class_precision.k_per_decade_for_pk = 10.;
    par->class_precision.start_large_k_at_tau_h_over_tau_k = 0.07;
    par->class_precision.l_max_g = 12;
    par->class_precision.l_max_ur = 17;
    par->class_precision.tol_perturb_integration = 1e-5;
    par->class_precision.radiation_streaming_trigger_tau_over_tau_k = 45.;
    par->class_precision.ur_fluid_trigger_tau_over_tau_k = 30.;

    par->only_cross_correlations = 0;

    const size_t separations_size = 100;
    par->sep = coffe_generate_range(5, 305, separations_size);
    par->sep_len = separations_size;

    par->interp_method = COFFE_INTERP_AKIMA;
    par->covariance_integration_method = COFFE_COVARIANCE_INTEGRATION_STANDARD;
    par->covariance_integration_bins = 8000;
    par->covariance_interpolation_method = COFFE_INTERP2D_BICUBIC;
    par->file_power_spectrum[0] = 0;
    /*
        turns out you can do this
        https://stackoverflow.com/a/1662202
    */
    double k[] = {
        #include "WAVENUMBER_HEADER.dat"
    };
    double pk[] = {
        #include "POWER_SPECTRUM_HEADER.dat"
    };

    for (size_t i = 0; i < COFFE_ARRAY_SIZE(k); ++i){
        k[i] *= par->h;
        pk[i] /= pow(par->h, 3);
    }

    coffe_new_spline(&par->power_spectrum);
    coffe_new_spline(&par->power_spectrum_norm);

    coffe_init_spline(
        &par->power_spectrum,
        k, pk, COFFE_ARRAY_SIZE(k),
        par->interp_method
    );

    coffe_new_spline2d(&par->power_spectrum2d);
    coffe_new_spline2d(&par->power_spectrum2d_norm);

    par->k_min = 1e-5;
    par->k_max = 300.;
    par->pk_type = COFFE_PK_LINEAR;
    par->midpoint_approximation = 1;

    {
        size_t len = par->power_spectrum.spline->size;
        double *k_norm =
            (double *)coffe_malloc(sizeof(double)*len);
        double *pk_norm =
            (double *)coffe_malloc(sizeof(double)*len);
        for (size_t i = 0; i < len; ++i){
            k_norm[i] = par->power_spectrum.spline->x[i];
            pk_norm[i] = par->power_spectrum.spline->y[i];
        }
        coffe_init_spline(
            &par->power_spectrum_norm,
            k_norm, pk_norm, len,
            par->interp_method
        );
        par->k_min_norm = par->k_min;
        par->k_max_norm = par->k_max;
        free(k_norm);
        free(pk_norm);
    }

    par->galaxy_bias_analytic = 0;

    coffe_new_spline(&par->galaxy_bias1);
    coffe_new_spline(&par->galaxy_bias2);
    coffe_new_spline(&par->galaxy_bias3);
    coffe_new_spline(&par->galaxy_bias4);
    coffe_new_spline(&par->magnification_bias1);
    coffe_new_spline(&par->magnification_bias2);
    coffe_new_spline(&par->evolution_bias1);
    coffe_new_spline(&par->evolution_bias2);

    parse_bias_default(
        1.0, &par->galaxy_bias1, par->interp_method
    );
    par->read_galaxy_bias1 = 0;
    par->file_galaxy_bias1[0] = 0;
    par->degree_galaxy_bias1 = 1;

    parse_bias_default(
        1.0, &par->galaxy_bias2, par->interp_method
    );
    par->read_galaxy_bias2 = 0;
    par->file_galaxy_bias2[0] = 0;
    par->degree_galaxy_bias2 = 1;

    parse_bias_default(
        1.0, &par->galaxy_bias3, par->interp_method
    );

    parse_bias_default(
        1.0, &par->galaxy_bias4, par->interp_method
    );

    parse_bias_default(
        0.0, &par->magnification_bias1, par->interp_method
    );
    par->read_magnification_bias1 = 0;
    par->file_magnification_bias1[0] = 0;
    par->degree_magnification_bias1 = 1;

    parse_bias_default(
        0.0, &par->magnification_bias2, par->interp_method
    );
    par->read_magnification_bias2 = 0;
    par->file_magnification_bias2[0] = 0;
    par->degree_magnification_bias2 = 1;

    parse_bias_default(
        0.0, &par->evolution_bias1, par->interp_method
    );
    par->read_evolution_bias1 = 0;
    par->file_evolution_bias1[0] = 0;
    par->degree_evolution_bias1 = 1;

    parse_bias_default(
        0.0, &par->evolution_bias2, par->interp_method
    );
    par->read_evolution_bias2 = 0;
    par->file_evolution_bias2[0] = 0;
    par->degree_evolution_bias2 = 1;

    par->output_type = MULTIPOLES;
    par->density1 = (double *)coffe_malloc(sizeof(double));
    par->density1[0] = 1e-3;
    par->density1_len = 1;
    par->density2 = (double *)coffe_malloc(sizeof(double));
    par->density2[0] = 1e-3;
    par->density2_len = 1;
    par->fsky = (double *)coffe_malloc(sizeof(double));
    par->fsky[0] = 1.0;
    par->fsky_len = 1;
    par->pixelsize = (double *)coffe_malloc(sizeof(double));
    par->pixelsize[0] = 10;
    par->pixelsize_len = 1;
    par->zmin = NULL;
    par->zmin_len = 0;
    par->zmax = NULL;
    par->zmax_len = 0;
    par->covariance_window = 0;
    par->covariance_pop1 = 1;
    par->covariance_pop2 = 1;
    par->covariance_pop3 = 1;
    par->covariance_pop4 = 1;
    par->covariance_cosmic = 1;
    par->covariance_mixed = 1;
    par->covariance_poisson = 1;

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

    par->mu = (double *)coffe_malloc(sizeof(double));
    par->mu[0] = 0.5;
    par->mu_len = 1;

    par->divergent = 0;

    par->type_bg = NULL;
    par->type_bg_len = 0;

    par->background_bins = 10000;
    par->bessel_bins = 10000;
    #ifndef HAVE_CUBA
    par->integration_method = 2;
    #endif
    par->integration_bins = 750000;

    par->z_mean = (double *)coffe_malloc(sizeof(double));
    par->deltaz = (double *)coffe_malloc(sizeof(double));
    par->z_mean[0] = 1.0;
    par->deltaz[0] = 0.2;
    par->z_mean_len = 1;
    par->deltaz_len = 1;

    par->zmin = (double *)coffe_malloc(sizeof(double));
    par->zmax = (double *)coffe_malloc(sizeof(double));
    par->zmin[0] = 0.9;
    par->zmax[0] = 1.1;
    par->zmin_len = 1;
    par->zmax_len = 1;

    const int multipoles[] = {0, 2, 4};
    par->multipole_values = (int *)coffe_malloc(
        sizeof(int) * COFFE_ARRAY_SIZE(multipoles)
    );
    for (size_t i = 0; i < COFFE_ARRAY_SIZE(multipoles); ++i)
        par->multipole_values[i] = multipoles[i];
    par->multipole_values_len = COFFE_ARRAY_SIZE(multipoles);

    par->flatsky_local = 0;

    par->flatsky_local_nonlocal = 0;

    par->flatsky_nonlocal = 0;

    par->verbose = 0;

    #ifndef COFFE_CYTHON
    par->conf = NULL;
    #endif

    par->flag = 1;

    return EXIT_SUCCESS;
}


/**
    parses all the settings from the input file
    (given by argv[1]) into the structure <par>
**/
#ifndef COFFE_CYTHON
int coffe_parser_init(
    char *filename,
    coffe_parameters_t *par
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
    parse_double(conf, "omega_gamma", &par->Omega0_gamma, COFFE_TRUE);

    par->Omega0_de = 1. - (par->Omega0_m + par->Omega0_gamma);

    /* mean redshift */
    if (
        par->output_type == CORRFUNC ||
        par->output_type == MULTIPOLES ||
        par->output_type == COVARIANCE_MULTIPOLES
    ){
        if (par->z_mean_len)
            free(par->z_mean);
        par->z_mean_len = 0;
        parse_double_array(conf, "z_mean", &par->z_mean, &par->z_mean_len);
        for (size_t i = 0; i < par->z_mean_len; ++i){
            if (par->z_mean[i] <= 0){
                print_error_verbose(PROG_VALUE_ERROR, "z_mean");
                exit(EXIT_FAILURE);
            }
        }
    }

    /* range of integration for redshift averaged multipoles */
    if (par->output_type == AVERAGE_MULTIPOLES){
        parse_double_array(conf, "z_min", &par->zmin, &par->zmin_len);
        parse_double_array(conf, "z_max", &par->zmax, &par->zmax_len);
    }

    /* the interpolation method for GSL */
    {
    int temp;
    parse_int(conf, "interpolation", &temp, COFFE_FALSE);
    par->interp_method = (enum coffe_interp1d_type)temp;
    }

    /* the cosine of the angle for the full sky correlation function */
    if (par->output_type == CORRFUNC){
        parse_double_array(conf, "mu", &par->mu, &par->mu_len);
    }

    /* the multipoles of the correlation function */
    if (
        par->output_type == MULTIPOLES ||
        par->output_type == AVERAGE_MULTIPOLES ||
        par->output_type == COVARIANCE_MULTIPOLES ||
        par->output_type == COVARIANCE_AVERAGE_MULTIPOLES
    ){
        if (par->multipole_values_len)
            free(par->multipole_values);
        par->multipole_values_len = 0;
        parse_int_array(conf, "multipoles", &par->multipole_values, &par->multipole_values_len);
    }

    /* the custom separations for the ang/full correlation function or multipoles */
    if (
        par->output_type == CORRFUNC ||
        par->output_type == MULTIPOLES ||
        par->output_type == AVERAGE_MULTIPOLES ||
        par->output_type == COVARIANCE_MULTIPOLES
    ){
        double xmin, xmax;
        parse_double(conf, "separations_min", &xmin, COFFE_TRUE);
        parse_double(conf, "separations_max", &xmax, COFFE_TRUE);
        double step = 0;
        int sampling = 0;
        const int error1 = parse_int(conf, "separations_sampling", &sampling, COFFE_FALSE);
        const int error2 = parse_double(conf, "separations_step", &step, COFFE_FALSE);
        if (
            (error1 == EXIT_SUCCESS) && (error2 == EXIT_SUCCESS)
        ){
            fprintf(stderr, "You need to specify one of `separations_sampling` and `separations_step`, but not both");
            exit(EXIT_FAILURE);
        }
        if (
            (error1 == EXIT_FAILURE) && (error2 == EXIT_FAILURE)
        ){
            fprintf(stderr, "You need to specify one of `separations_sampling` and `separations_step`");
            exit(EXIT_FAILURE);
        }
        if (par->sep_len)
            free(par->sep);
        par->sep = coffe_generate_range(xmin, xmax, sampling);
        par->sep_len = (size_t)sampling;
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
        par->covariance_integration_method != COFFE_COVARIANCE_INTEGRATION_STANDARD &&
        par->covariance_integration_method != COFFE_COVARIANCE_INTEGRATION_FFTLOG
    ){
        print_error_verbose(PROG_VALUE_ERROR, "covariance_integration_method");
        exit(EXIT_FAILURE);
    }

    if (par->covariance_integration_method == COFFE_COVARIANCE_INTEGRATION_FFTLOG){
        /* integration method for covariance (only for FFT log) */
        parse_int(conf, "covariance_integration_bins", &par->covariance_integration_bins, COFFE_TRUE);
        if (par->covariance_integration_bins <= 0){
            print_error_verbose(PROG_VALUE_ERROR, "covariance_integration_bins");
            exit(EXIT_FAILURE);
        }
        {
        int temp;
        parse_int(conf, "covariance_interpolation_method", &temp, COFFE_TRUE);
        par->covariance_interpolation_method = (enum coffe_interp2d_type)temp;
        }
        if (
            par->covariance_interpolation_method != COFFE_INTERP2D_BILINEAR &&
            par->covariance_interpolation_method != COFFE_INTERP2D_BICUBIC
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
        parse_int(conf, "degree_galaxy_bias1", &par->degree_galaxy_bias1, COFFE_FALSE);
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
        parse_int(conf, "degree_galaxy_bias2", &par->degree_galaxy_bias2, COFFE_FALSE);
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

        coffe_init_spline(
            &par->galaxy_bias3,
            x, y, bins,
            par->interp_method
        );

        coffe_init_spline(
            &par->galaxy_bias4,
            x, y, bins,
            par->interp_method
        );

        free(x);
        free(y);
    }

    /* parsing the magnification bias (s) */
    parse_int(conf, "read_magnification_bias1", &par->read_magnification_bias1, COFFE_FALSE);
    parse_int(conf, "degree_magnification_bias1", &par->degree_magnification_bias1, COFFE_FALSE);
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
    parse_int(conf, "degree_magnification_bias2", &par->degree_magnification_bias2, COFFE_FALSE);
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
    parse_int(conf, "degree_evolution_bias1", &par->degree_evolution_bias1, COFFE_FALSE);
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
    parse_int(conf, "degree_evolution_bias2", &par->degree_evolution_bias2, COFFE_FALSE);
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
    if (
        par->output_type == COVARIANCE_MULTIPOLES ||
        par->output_type == COVARIANCE_AVERAGE_MULTIPOLES
    ){
        parse_double_array(
            conf,
            "density1",
            &par->density1,
            &par->density1_len
        );
        parse_double_array(
            conf,
            "density2",
            &par->density2,
            &par->density2_len
        );
        parse_double_array(
            conf,
            "fsky",
            &par->fsky,
            &par->fsky_len
        );
        parse_double_array(
            conf,
            "pixelsize",
            &par->pixelsize,
            &par->pixelsize_len
        );
        parse_int(
            conf,
            "covariance_window",
            &par->covariance_window,
            COFFE_TRUE
        );
    }

    if (par->output_type == COVARIANCE_MULTIPOLES){
        parse_double_array(
            conf,
            "z_mean",
            &par->z_mean,
            &par->z_mean_len
        );
        parse_double_array(
            conf,
            "deltaz",
            &par->deltaz,
            &par->deltaz_len
        );
        parse_double_array(
            conf,
            "fsky",
            &par->fsky,
            &par->fsky_len
        );

        if (
            par->density1_len != par->fsky_len ||
            par->density1_len != par->density2_len ||
            par->density1_len != par->z_mean_len ||
            par->density1_len != par->deltaz_len ||
            par->density1_len != par->pixelsize_len
        ){
            fprintf(
                stderr,
                "ERROR: covariance parameters have mismatching lengths; "
                "please ensure they have the same length!\n");
            exit(EXIT_FAILURE);
        }
    }

    if (par->output_type == COVARIANCE_AVERAGE_MULTIPOLES){
        parse_double_array(
            conf,
            "zmin",
            &par->zmin,
            &par->zmin_len
        );
        parse_double_array(
            conf,
            "zmax",
            &par->zmax,
            &par->zmax_len
        );
        if (
            par->density1_len != par->fsky_len ||
            par->density1_len != par->density2_len ||
            par->density1_len != par->zmin_len ||
            par->density1_len != par->zmax_len
        ){
            fprintf(
                stderr,
                "ERROR: covariance parameters have mismatching lengths; "
                "please ensure they have the same length!\n");
            exit(EXIT_FAILURE);
        }
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
        {
        int temp;
        parse_int(conf, "pk_type", &temp, COFFE_TRUE);
        par->pk_type = (enum coffe_pk_type)temp;
        }
        parse_double(conf, "h", &par->h, COFFE_TRUE);
        const int error_sigma8 = parse_double(conf, "sigma8", &par->sigma8, COFFE_FALSE);
        const int error_A_s = parse_double(conf, "A_s", &par->A_s, COFFE_FALSE);

        if (error_sigma8 == EXIT_SUCCESS && error_A_s == EXIT_SUCCESS){
            fprintf(
                stderr,
                "ERROR: unable to set both sigma8 and A_s simultaneously"
            );
            exit(EXIT_FAILURE);
        }

        parse_double(conf, "n_s", &par->n_s, COFFE_TRUE);
        parse_double(conf, "k_pivot", &par->k_pivot, COFFE_TRUE);
        // new stuff for forecast
        parse_double(conf, "T_cmb", &par->T_cmb, COFFE_FALSE);
        parse_double(conf, "N_ur", &par->N_ur, COFFE_FALSE);
        parse_double(conf, "m_ncdm", &par->m_ncdm, COFFE_FALSE);
        parse_int(conf, "N_ncdm", &par->N_ncdm, COFFE_FALSE);
        parse_double(conf, "YHe", &par->YHe, COFFE_FALSE);
        par->Omega0_nu = par->m_ncdm / 93.14 / par->h / par->h;
        par->Omega0_cdm = par->Omega0_m - par->Omega0_baryon - par->Omega0_nu;
        parse_external_power_spectrum(par);
    }
    else{
#endif
        par->Omega0_cdm = par->Omega0_m - par->Omega0_baryon - par->Omega0_nu;

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
            k_norm[i] = par->power_spectrum.spline->x[i];
            pk_norm[i] = par->power_spectrum.spline->y[i];
        }
        coffe_init_spline(
            &par->power_spectrum_norm,
            k_norm, pk_norm, len,
            par->interp_method
        );
        par->k_min_norm = par->k_min;
        par->k_max_norm = par->k_max;
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
#endif
