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
#include <time.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>

#include "common.h"
#include "background.h"



/**
    fits some data to a polynomial of some degree
**/

int coffe_fit_polynomial(
    const double *data_x,
    const double *data_y,
    const size_t data_size,
    const size_t degree_current,
    double **coefficients,
    size_t *coefficients_size,
    double *chisq
)
{
    /* by default, try a quadratic fit */
    gsl_matrix *X = gsl_matrix_alloc(data_size, degree_current + 1);
    gsl_vector *y = gsl_vector_alloc(data_size);
    gsl_vector *c = gsl_vector_alloc(degree_current + 1);
    gsl_matrix *cov = gsl_matrix_alloc(degree_current + 1, degree_current + 1);

    for (size_t i = 0; i < data_size; ++i){
        const double xi = data_x[i];
        const double yi = data_y[i];

        for (size_t j = 0; j < degree_current + 1; ++j)
            gsl_matrix_set(X, i, j, pow(xi, j));

        gsl_vector_set(y, i, yi);
    }

    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(data_size, degree_current + 1);

    gsl_multifit_linear(X, y, c, cov, chisq, work);
    gsl_multifit_linear_free(work);

    *coefficients_size = degree_current + 1;
    *coefficients = (double *)coffe_malloc(sizeof(double) * *coefficients_size);

    for (size_t i = 0; i < *coefficients_size; ++i)
        (*coefficients)[i] = gsl_vector_get(c, i);

    gsl_matrix_free(cov);
    gsl_matrix_free(X);
    gsl_vector_free(y);
    gsl_vector_free(c);

    return EXIT_SUCCESS;
}


int coffe_get_bias_coefficients(
    const coffe_interpolation *comoving_distance,
    const coffe_interpolation *z_as_chi,
    const double *z_mean,
    const size_t z_mean_size,
    const double *sep,
    const size_t sep_size,
    const coffe_interpolation *bias,
    const int degree,
    coffe_fit_coefficients_array_t *bias_coefficients
)
{
    bias_coefficients->size = z_mean_size;
    bias_coefficients->array = (coffe_fit_coefficients_t *)coffe_malloc(
        sizeof(coffe_fit_coefficients_t) * bias_coefficients->size
    );

    /* set everything to NULL so we don't free it accidentally */
    for (size_t i = 0; i < z_mean_size; ++i)
        coffe_new_fit_coefficients(&bias_coefficients->array[i]);

    for (size_t i = 0; i < z_mean_size; ++i){
        const double z_min = coffe_interp_spline(
            z_as_chi,
            coffe_interp_spline(comoving_distance, z_mean[i])
            -
            sep[sep_size - 1]
        );
        const double z_max = coffe_interp_spline(
            z_as_chi,
            coffe_interp_spline(comoving_distance, z_mean[i])
            +
            sep[sep_size - 1]
        );
        const double z_mean = (z_min + z_max) / 2;

        /* sampling the input spline */
        /* note that we cannot sample beyond z = 15 (hardcoded for now) */
        const size_t num = bias->spline->size;
        double *x = (double *)coffe_malloc(sizeof(double) * num);
        double *y = (double *)coffe_malloc(sizeof(double) * num);
        for (size_t j = 0; j < num; ++j){
            const double z = bias->spline->x[j];
            x[j] = coffe_interp_spline(comoving_distance, z)
                 - coffe_interp_spline(comoving_distance, z_mean);
            y[j] = coffe_interp_spline(bias, z);
        }

        bias_coefficients->array[i].degree = degree;
        bias_coefficients->array[i].z_min = z_min;
        bias_coefficients->array[i].z_max = z_max;
        coffe_fit_polynomial(
            x, y, num,
            degree,
            &bias_coefficients->array[i].coefficients,
            &bias_coefficients->array[i].size,
            &bias_coefficients->array[i].chisq
        );
        free(x);
        free(y);
    }

    return EXIT_SUCCESS;
}


struct integration_params
{
    double Omega0_cdm; /* omega parameter for cold dark matter */

    double Omega0_baryon; /* omega parameter for baryons */

    double Omega0_gamma; /* omega parameter for photons */

    double Omega0_de; /* present omega parameter for dark energy-like component */

    coffe_interpolation w; /* interpolator of w(z) */

    coffe_interpolation wint; /* result of exp(3*int((1 + w(z))/(1 + z))) */

    coffe_interpolation xint; /* result of Omega0_m/(1 - Omega0_m)*exp(-3*int(w(a)/a)) */
};


/* only needed here */

struct temp_background
{
    double *z; /* redshift */

    double *a; /* scale factor (normalized so that today a=1) */

    double *Hz; /* hubble parameter H(z) (in 1 / Mpc) */

    double *conformal_Hz; /* conformal hubble parameter (in 1 / Mpc) */

    double *conformal_Hz_prime; /* derivative of conformal hubble parameter wrt conformal time (in 1 / Mpc^2) */

    double *D1; /* growth rate D_1(a) */

    double *f; /* growth function f=d(log D)/d(log a) */

    double *G1, *G2;

    double *comoving_distance; /* comoving distance (in Mpc) */

};


/**
    differential equation for the growth rate D_1
**/

static int growth_rate_ode(
    double a,
    const double y[],
    double f[],
    void *params
)
{
    struct integration_params *par = (struct integration_params *) params;
    double z = 1./a - 1;
    double w = coffe_interp_spline(&par->w, z);
    double x = coffe_interp_spline(&par->xint, z);

    f[0] = y[1];
    f[1] = -3./2*(1 - w/(1 + x))*y[1]/a
           +
            3./2*x/(1 + x)*y[0]/pow(a, 2);
    return GSL_SUCCESS;
}


/**
    jacobian for the growth rate differential equation
**/

static int growth_rate_jac(
    double a,
    const double y[],
    double *dfdy,
    double dfdt[],
    void *params
)
{
    struct integration_params *par = (struct integration_params *) params;
    gsl_matrix_view dfdy_mat =
        gsl_matrix_view_array(dfdy, 2, 2);
    double z = 1./a - 1;
    double w = coffe_interp_spline(&par->w, z);
    double x = coffe_interp_spline(&par->xint, z);
    double x_der =
        gsl_spline_eval_deriv(par->xint.spline, z, par->xint.accel);
    gsl_matrix *m = &dfdy_mat.matrix;
    gsl_matrix_set(m, 0, 0, 0.0);
    gsl_matrix_set(m, 0, 1, 1.0);
    gsl_matrix_set(m, 1, 0, 3./2./pow(a, 2)*(x/(1 + x)));
    gsl_matrix_set(m, 1, 1, -3./a/2./(1. - w/(1 + x)));
    dfdt[0] = 0.0;
    dfdt[1] =
        3*y[1]/pow(a, 2)/2./(1 - w/(1 + x))
       -3*y[0]/pow(a, 3)*(x/(1 + x))
       -3*y[0]/2./pow(a, 2)*(x*x_der/pow(1 + x, 2))
       +3*y[0]/2./pow(a, 2)*(x_der/(1 + x));
    return GSL_SUCCESS;
}


/**
    integrand of (1 + w)/(1 + z)
**/

static double integrand_w(
    double z,
    void *p
)
{
    struct integration_params *par = (struct integration_params *) p;
    return (1 + coffe_interp_spline(&par->w, z))/(1 + z);
}

/**
    integrand w(a)/a
**/

static double integrand_x(
    double a,
    void *p
)
{
    struct integration_params *par = (struct integration_params *) p;
    double result = coffe_interp_spline(&par->w, 1/a - 1)/a;
    return result;
}


/**
    integrand of the comoving distance
**/

static double integrand_comoving(
    double z,
    void *p
)
{
    struct integration_params *par = (struct integration_params *) p;
    double integrand = 1./pow(
        (par->Omega0_cdm + par->Omega0_baryon)*pow(1 + z, 3)
       +par->Omega0_gamma*pow(1 + z, 4)
       +par->Omega0_de*coffe_interp_spline(&par->wint, z),
        1./2);
    return integrand;
}

/**
    computes and stores all the background functions
**/

int coffe_background_init(
    coffe_parameters_t *par,
    coffe_background_t *bg
)
{
    clock_t start, end;
    start = clock();

    if (par->verbose)
        printf("Initializing the background...\n");

    /* make sure we don't have memory leaks */
    coffe_background_free(bg);

    gsl_error_handler_t *default_handler =
        gsl_set_error_handler_off();

    struct temp_background *temp_bg =
        (struct temp_background *)coffe_malloc(sizeof(struct temp_background));
    temp_bg->z = (double *)coffe_malloc(sizeof(double)*par->background_bins);
    temp_bg->a = (double *)coffe_malloc(sizeof(double)*par->background_bins);
    temp_bg->Hz = (double *)coffe_malloc(sizeof(double)*par->background_bins);
    temp_bg->conformal_Hz = (double *)coffe_malloc(sizeof(double)*par->background_bins);
    temp_bg->conformal_Hz_prime = (double *)coffe_malloc(sizeof(double)*par->background_bins);
    temp_bg->D1 = (double *)coffe_malloc(sizeof(double)*par->background_bins);
    temp_bg->f = (double *)coffe_malloc(sizeof(double)*par->background_bins);
    temp_bg->G1 = (double *)coffe_malloc(sizeof(double)*par->background_bins);
    temp_bg->G2 = (double *)coffe_malloc(sizeof(double)*par->background_bins);
    temp_bg->comoving_distance = (double *)coffe_malloc(sizeof(double)*par->background_bins);

    struct integration_params ipar;
    ipar.Omega0_cdm = par->Omega0_cdm;
    ipar.Omega0_baryon = par->Omega0_baryon;
    ipar.Omega0_gamma = par->Omega0_gamma;
    ipar.Omega0_de = par->Omega0_de;
    {
        double z, z_max = 100.;
        const size_t bins = 16384;
        double *z_array = (double *)coffe_malloc(sizeof(double)*(bins + 1));

        double *w_array = (double *)coffe_malloc(sizeof(double)*(bins + 1));
        for (size_t i = 0; i <= bins; ++i){
            z = z_max*i/(double)bins;
            z_array[i] = z;
            w_array[i] = coffe_dark_energy_eos(z, par);
        }
        coffe_init_spline(&ipar.w, z_array, w_array, bins + 1, COFFE_INTERP_LINEAR);
        free(w_array);

        double *wint_array = (double *)coffe_malloc(sizeof(double)*(bins + 1));
        for (size_t i = 0; i <= bins; ++i){
            z = z_max*i/(double)bins;
            z_array[i] = z;
            wint_array[i] = exp(
                3 * coffe_integrate_1d(
                    &integrand_w,
                    &ipar,
                    0,
                    z
                )
            );
        }
        coffe_init_spline(&ipar.wint, z_array, wint_array, bins + 1, COFFE_INTERP_LINEAR);
        free(wint_array);

        double *xint_array = (double *)coffe_malloc(sizeof(double)*(bins + 1));
        for (size_t i = 0; i <= bins; ++i){
            z = z_max*i/(double)bins;
            if (i != 0){
                z_array[i] = z;
                xint_array[i] =
                    (ipar.Omega0_cdm + ipar.Omega0_baryon)
                   /(1 - (ipar.Omega0_cdm + ipar.Omega0_baryon))
                   *exp(
                        -3 * coffe_integrate_1d(
                            &integrand_x,
                            &ipar,
                            1. / (1 + z),
                            1.
                        )
                    );
            }
            else{
                z_array[i] = 0;
                xint_array[i] =
                    (ipar.Omega0_cdm + ipar.Omega0_baryon)
                   /(1 - (ipar.Omega0_cdm + ipar.Omega0_baryon));
            }
        }
        coffe_init_spline(&ipar.xint, z_array, xint_array, bins + 1, COFFE_INTERP_LINEAR);
        free(xint_array);

        free(z_array);
    }

    double h = 1E-6;

    double a_initial, z, w, wint;
    gsl_odeiv2_system sys =
        {growth_rate_ode, growth_rate_jac, 2, &ipar};

    const gsl_odeiv2_step_type *step_type =
        gsl_odeiv2_step_rk2;

    gsl_odeiv2_step *step =
        gsl_odeiv2_step_alloc(step_type, 2);
    gsl_odeiv2_control *control =
        gsl_odeiv2_control_y_new(0.0, h);
    gsl_odeiv2_evolve *evolve =
        gsl_odeiv2_evolve_alloc(2);

    /* initial values for the differential equation (D_1 and D_1') */
    double initial_values[2] = {0.05, 1.0};

    for (int i = 0; i < par->background_bins; ++i){
        a_initial = 0.05;
        initial_values[0] = 0.05;
        initial_values[1] = 1.0;

        z = 15.*i/(double)(par->background_bins - 1);

        w = coffe_interp_spline(&ipar.w, z);
        wint = coffe_interp_spline(&ipar.wint, z);

        (temp_bg->z)[i] = z;
        (temp_bg->a)[i] = 1./(1. + z);
        (temp_bg->Hz)[i] = sqrt(
             (par->Omega0_cdm + par->Omega0_baryon)*pow(1 + z, 3)
            +par->Omega0_gamma*pow(1 + z, 4)
            +par->Omega0_de*wint) * par->h * COFFE_H0; // in units 1 / Mpc
        (temp_bg->conformal_Hz)[i] = (temp_bg->a)[i]*(temp_bg->Hz)[i]; // in units 1 / Mpc
        (temp_bg->conformal_Hz_prime)[i] = -(
            pow(1 + z, 3)
           *(
                2 * (1 + z) * par->Omega0_gamma
                +
                (par->Omega0_cdm + par->Omega0_baryon)
            )
            +
            (1 + 3 * w) * par->Omega0_de * wint
        )/pow(1 + z, 2)/2. * pow(par->h * COFFE_H0, 2); // in units 1 / Mpc^2

        while (a_initial < (temp_bg->a)[i]){
            gsl_odeiv2_evolve_apply(
                evolve, control, step,
                &sys, &a_initial, (temp_bg->a)[i],
                &h, initial_values
            );
        }

        (temp_bg->D1)[i] = initial_values[0];
        (temp_bg->f)[i] = initial_values[1] * (temp_bg->a)[i] / (temp_bg->D1)[i];

        (temp_bg->comoving_distance)[i] = coffe_integrate_1d(
            &integrand_comoving,
            &ipar,
            0.,
            z
        ) / (par->h * COFFE_H0); // in units Mpc

        if (z > 1E-10){
            (temp_bg->G1)[i] =
                (temp_bg->conformal_Hz_prime)[i]
                    /pow((temp_bg->conformal_Hz)[i], 2)
                    +
                    (2 - 5*coffe_interp_spline(&par->magnification_bias1, z))
                   /(
                        (temp_bg->comoving_distance[i])
                       *(temp_bg->conformal_Hz)[i]
                    )
                    +
                    5*coffe_interp_spline(&par->magnification_bias1, z)
                    -
                    coffe_interp_spline(&par->evolution_bias1, z);

            (temp_bg->G2)[i] =
                (temp_bg->conformal_Hz_prime)[i]
                    /pow((temp_bg->conformal_Hz)[i], 2)
                    +
                    (2 - 5*coffe_interp_spline(&par->magnification_bias2, z))
                   /(
                        (temp_bg->comoving_distance[i])
                       *(temp_bg->conformal_Hz)[i]
                    )
                    +
                    5*coffe_interp_spline(&par->magnification_bias2, z)
                    -
                    coffe_interp_spline(&par->evolution_bias2, z);
        }
        else{
            (temp_bg->G1)[i] = 0;
            (temp_bg->G2)[i] = 0;
        }
    }
    /* normalizing D1 so it's 1 at z = 0 */
    {
        const double D10 = temp_bg->D1[0];
        for (int i = 0; i < par->background_bins; ++i)
            temp_bg->D1[i] /= D10;
    }

    /* initializing the splines; all splines are a function of z */
    coffe_init_spline(
        &bg->a,
        temp_bg->z,
        temp_bg->a,
        par->background_bins,
        par->interp_method
    );
    coffe_init_spline(
        &bg->Hz,
        temp_bg->z,
        temp_bg->Hz,
        par->background_bins,
        par->interp_method
    );
    coffe_init_spline(
        &bg->conformal_Hz,
        temp_bg->z,
        temp_bg->conformal_Hz,
        par->background_bins,
        par->interp_method
    );
    coffe_init_spline(
        &bg->conformal_Hz_prime,
        temp_bg->z,
        temp_bg->conformal_Hz_prime,
        par->background_bins,
        par->interp_method
    );
    coffe_init_spline(
        &bg->D1,
        temp_bg->z,
        temp_bg->D1,
        par->background_bins,
        par->interp_method
    );
    coffe_init_spline(
        &bg->f,
        temp_bg->z,
        temp_bg->f,
        par->background_bins,
        par->interp_method
    );
    coffe_init_spline(
        &bg->comoving_distance,
        temp_bg->z,
        temp_bg->comoving_distance,
        par->background_bins,
        par->interp_method
    );
    coffe_init_spline(
        &bg->G1,
        temp_bg->z,
        temp_bg->G1,
        par->background_bins,
        par->interp_method
    );
    coffe_init_spline(
        &bg->G2,
        temp_bg->z,
        temp_bg->G2,
        par->background_bins,
        par->interp_method
    );

    /* inverse of the z, chi(z) spline (only one we need to invert) */
    coffe_init_spline(
        &bg->z_as_chi,
        temp_bg->comoving_distance,
        temp_bg->z,
        par->background_bins,
        par->interp_method
    );

    /* memory cleanup */
    free(temp_bg->z);
    free(temp_bg->a);
    free(temp_bg->Hz);
    free(temp_bg->conformal_Hz);
    free(temp_bg->conformal_Hz_prime);
    free(temp_bg->D1);
    free(temp_bg->f);
    free(temp_bg->G1);
    free(temp_bg->G2);
    free(temp_bg->comoving_distance);
    free(temp_bg);
    gsl_odeiv2_step_free(step);
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_evolve_free(evolve);
    gsl_spline_free(ipar.w.spline), gsl_interp_accel_free(ipar.w.accel);
    gsl_spline_free(ipar.wint.spline), gsl_interp_accel_free(ipar.wint.accel);
    gsl_spline_free(ipar.xint.spline), gsl_interp_accel_free(ipar.xint.accel);

    /* this should maybe be placed in some other module as to not bloat this one */
    if (par->flatsky_nonlocal && par->z_mean_len){

        coffe_get_bias_coefficients(
            &bg->comoving_distance, &bg->z_as_chi,
            par->z_mean, par->z_mean_len,
            par->sep, par->sep_len,
            &par->galaxy_bias1, par->degree_galaxy_bias1, &par->galaxy_bias1_coefficients
        );

        coffe_get_bias_coefficients(
            &bg->comoving_distance, &bg->z_as_chi,
            par->z_mean, par->z_mean_len,
            par->sep, par->sep_len,
            &par->galaxy_bias2, par->degree_galaxy_bias2, &par->galaxy_bias2_coefficients
        );

        coffe_get_bias_coefficients(
            &bg->comoving_distance, &bg->z_as_chi,
            par->z_mean, par->z_mean_len,
            par->sep, par->sep_len,
            &par->magnification_bias1, par->degree_magnification_bias1, &par->magnification_bias1_coefficients
        );

        coffe_get_bias_coefficients(
            &bg->comoving_distance, &bg->z_as_chi,
            par->z_mean, par->z_mean_len,
            par->sep, par->sep_len,
            &par->magnification_bias2, par->degree_magnification_bias2, &par->magnification_bias2_coefficients
        );

        coffe_get_bias_coefficients(
            &bg->comoving_distance, &bg->z_as_chi,
            par->z_mean, par->z_mean_len,
            par->sep, par->sep_len,
            &par->evolution_bias1, par->degree_evolution_bias1, &par->evolution_bias1_coefficients
        );

        coffe_get_bias_coefficients(
            &bg->comoving_distance, &bg->z_as_chi,
            par->z_mean, par->z_mean_len,
            par->sep, par->sep_len,
            &par->evolution_bias2, par->degree_evolution_bias2, &par->evolution_bias2_coefficients
        );

    }

    gsl_set_error_handler(default_handler);

    end = clock();

    if (par->verbose)
        printf("Background initialized in %.2f s\n",
            (double)(end - start) / CLOCKS_PER_SEC);

    bg->flag = 1;

    return EXIT_SUCCESS;
}

int coffe_background_free(
    coffe_background_t *bg
)
{
    if (bg->flag){
        coffe_free_spline(&bg->z_as_chi);
        coffe_free_spline(&bg->a);
        coffe_free_spline(&bg->Hz);
        coffe_free_spline(&bg->conformal_Hz);
        coffe_free_spline(&bg->conformal_Hz_prime);
        coffe_free_spline(&bg->D1);
        coffe_free_spline(&bg->f);
        coffe_free_spline(&bg->G1);
        coffe_free_spline(&bg->G2);
        coffe_free_spline(&bg->comoving_distance);

    }
    bg->flag = 0;

    return EXIT_SUCCESS;
}
