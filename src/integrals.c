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
#include <time.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>

#ifdef HAVE_DOUBLE_EXPONENTIAL
#include "tanhsinh.h"
#endif

#include "common.h"
#include "background.h"
#include "integrals.h"
#include "twofast.h"

#ifdef HAVE_CLASS
#include "class.h"
#endif


#ifndef NORM
#define NORM(X) (X*COFFE_H0)
#endif

/**
    parameters for the integrand of the form P(k) k^2 j_l(k r)/(k r)^n
**/

struct integrals_params
{
    struct coffe_interpolation *result;
    double r;
    int n, l;
    enum coffe_integer_state state_n, state_l;
};


/**
    structure needed for integration of renormalization term
**/

struct integrals_divergent_params
{
    struct coffe_interpolation *result;
    double chi1, chi2;
    /* placeholder: in flat-sky, the non-renormalizable terms are nonzero */
    enum coffe_integer_state state;
};


/**
    list of separations (in Mpc/h!) for the integral I^4_0 so it doesn't need to rely
    on user input
**/

static const double coffe_sep[] = {
    0.000001, 0.000005, 0.0001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1., 1.5,
    2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20.,
    21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 32., 34., 36., 38., 40., 42., 44., 46.,
    48., 50., 52., 54., 56., 58., 60., 62., 64., 66., 68., 70., 72., 74., 76., 78., 80., 81., 82., 83.,
    84., 85., 86., 87., 88., 89., 90., 91., 92., 93., 94., 95., 95.5, 96., 96.5, 97., 97.5, 98.,
    98.5, 99., 99.5, 100., 100.5, 101., 101.5, 102., 102.5, 103., 103.5, 104., 104.5, 105.,
    106., 107., 108., 109., 110., 112., 114., 116., 118., 120., 124., 128., 132., 136., 140.,
    144., 148., 152., 156., 160., 164., 168., 172., 176., 180., 185., 190., 195., 200., 205.,
    210., 215., 220., 225., 230., 235., 240., 250., 260., 270., 280., 290., 300., 320., 340.,
    360., 380., 400., 420., 440., 460., 480., 500., 520., 540., 560., 580., 600., 620., 640.,
    660., 680., 700., 750., 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350,
    1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2050, 2100,
    2150, 2200, 2250, 2300, 2350, 2400, 2450, 2500, 2550, 2600, 2650, 2700, 2750, 2800, 2850,
    2900, 2950, 3000, 3050, 3100, 3150, 3200, 3250, 3300, 3350, 3400, 3450, 3500, 3550, 3600,
    3650, 3700, 3750, 3800, 3850, 3900, 3950, 4000, 4050, 4100, 4150, 4200, 4250, 4300, 4350,
    4400, 4450, 4500, 4550, 4600, 4650, 4700, 4750, 4800, 4850, 4900, 4950, 5000, 5050, 5100,
    5150, 5200, 5250, 5300, 5350, 5400, 5450, 5500, 5550, 5600, 5650, 5700, 5750, 5800, 5850,
    5900, 5950, 6000, 6050, 6100, 6150, 6200, 6250, 6300, 6350, 6400, 6450, 6500, 6550, 6600,
    7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000
};


/**
    length of the above
**/

static const size_t coffe_sep_len = sizeof(coffe_sep)/sizeof(coffe_sep[0]);


/**
    list of separations below 1 Mpc/h for which we need to evaluate
    the convergent integrals
**/

static const double min_sep[] = {
    NORM(1E-6), NORM(1E-5), NORM(1E-4),
    NORM(1E-3), NORM(2E-3), NORM(5E-3),
    NORM(7E-3), NORM(8E-3), NORM(9E-3),
    NORM(1E-2), NORM(1.1E-2), NORM(1.2E-2),
    NORM(1.25E-2), NORM(1.3E-2), NORM(1.35E-2),
    NORM(1.5E-2), NORM(2E-2), NORM(2.5E-2),
    NORM(3E-2), NORM(3.5E-2), NORM(5E-2),
    NORM(7E-2), NORM(8E-2), NORM(9E-2),
    NORM(1E-1), NORM(1.1E-1), NORM(1.2E-1),
    NORM(1.25E-1), NORM(1.3E-1), NORM(1.35E-1),
    NORM(1.5E-1), NORM(1.75E-1), NORM(2E-1), NORM(2.5E-1),
    NORM(3E-1), NORM(3.5E-1), NORM(4E-1),
    NORM(5E-1), NORM(6E-1),
    NORM(7E-1), NORM(8E-1), NORM(9E-1)
};


/**
    length of the above
**/

static const size_t min_sep_len = sizeof(min_sep) / sizeof(min_sep[0]);


/**
    find a given integral
**/

struct coffe_integral_t *coffe_find_integral(
    const struct coffe_integral_array_t *integral,
    const int n,
    const int l,
    const enum coffe_integer_state state_n,
    const enum coffe_integer_state state_l
)
{
    if (integral != NULL){
        if (integral->size != 0){
            for (size_t i = 0; i < integral->size; ++i){
                if (
                    integral->array[i].n == n &&
                    integral->array[i].l == l &&
                    integral->array[i].state_n == state_n &&
                    integral->array[i].state_l == state_l
                )
                    return &integral->array[i];
            }
            /* hasn't found anything, return NULL */
            return NULL;
        }
        else{
            return NULL;
        }
    }
    else{
        return NULL;
    }
}


/**
    checks if n and l are ints, and returns the local ones (half-ints if necessary)
**/

static int integrals_check_parameters(
    const int n,
    const int l,
    const enum coffe_integer_state state_n,
    const enum coffe_integer_state state_l,
    double *local_n,
    double *local_l
)
{
    if (state_n == COFFE_HALF_INTEGER)
        *local_n = (double)n / 2.;
    else
        *local_n = n;
    if (state_l == COFFE_HALF_INTEGER)
        *local_l = (double)l / 2.;
    else
        *local_l = l;

    return EXIT_SUCCESS;
}


/**
    coefficients of J_l and j_l for the r -> 0 limit
**/

static double integrals_coefficients(
    const int l,
    const enum coffe_integer_state state
)
{
    double result;
    double local_l;
    if (state == COFFE_HALF_INTEGER){
        local_l = (double)l / 2.;
        result = 1. / pow(2, local_l + 0.5) / gsl_sf_gamma(local_l + 1.5);
    }
    else{
        local_l = (double)l;
        result = sqrt(M_PI / 2.) / pow(2, local_l + 0.5) / gsl_sf_gamma(local_l + 1.5);
    }
    return result;
}


/**
    integrand k^(2 + l - n)*P(k)
**/

static double integrals_prefactor(
    double k,
#ifdef HAVE_DOUBLE_EXPONENTIAL
    const void *p
#else
    void *p
#endif
)
{
    struct integrals_params *test = (struct integrals_params *) p;

    double local_n, local_l;

    integrals_check_parameters(
        test->n,
        test->l,
        test->state_n,
        test->state_l,
        &local_n,
        &local_l
    );

    double result;

    if (
        test->state_n == COFFE_HALF_INTEGER &&
        test->state_l == COFFE_HALF_INTEGER
    )
        result = sqrt(M_PI / 2.)
          *pow(k, 2 + local_l - local_n)
          *coffe_interp_spline(test->result, k);
    else
        result = pow(k, 2 + local_l - local_n)
          *coffe_interp_spline(test->result, k);
    return result;
}


/**
    integrand of:
        k^2 P(k) j_l(k r)/(k r)^n
    or:
        k^2 P(k) J_l(k r)/(k r)^n
    if `state` == `coffe_integer_state`
**/

static double integrals_bessel_integrand(
    double k,
#ifdef HAVE_DOUBLE_EXPONENTIAL
    const void *p
#else
    void *p
#endif
)
{
    const struct integrals_params *integrand = (const struct integrals_params *) p;

    double local_n, local_l;

    integrals_check_parameters(
        integrand->n,
        integrand->l,
        integrand->state_n,
        integrand->state_l,
        &local_n,
        &local_l
    );

    double bessel_function;
    if (integrand->state_l == COFFE_HALF_INTEGER)
        bessel_function = sqrt(M_PI / 2.) * gsl_sf_bessel_Jn(
            (int)((integrand->l + 1) / 2),
            k * integrand->r
        )
       /pow(
            k * integrand->r,
            local_n + 0.5
        );
    else
        bessel_function = gsl_sf_bessel_jl(
            integrand->l,
            k * integrand->r
        )
       /pow(
            k * integrand->r,
            local_n
        );

    return k * k
            /* usually P(k) */
           *coffe_interp_spline(
                integrand->result,
                k
            )
           *bessel_function;
}


/**
    integrates any function (such as the one above)
    and returns the integral divided by 2 pi^2
**/

static double integrals_integrate_function(
    double (*func)(
        double,
#ifdef HAVE_DOUBLE_EXPONENTIAL
        const void*
#else
        void*
#endif
    ),
    const struct coffe_interpolation *input,
    const int n,
    const int l,
    const enum coffe_integer_state state_n,
    const enum coffe_integer_state state_l,
    const double sep,
    const double kmin,
    const double kmax
)
{
    struct integrals_params test;
    test.result = (struct coffe_interpolation *) input;
    test.n = n;
    test.l = l;
    test.state_n = state_n;
    test.state_l = state_l;
    test.r = sep;

    double local_n, local_l;

    integrals_check_parameters(
        test.n,
        test.l,
        test.state_n,
        test.state_l,
        &local_n,
        &local_l
    );

    const double result = coffe_integrate_1d(
        func,
        &test,
        kmin,
        kmax
    );

    double output;
    /* the case sep = 0 would just give us zero, so we skip it */
    if (local_n >= local_l && sep > 0)
        output = result * pow(sep, local_n - local_l);
    else
        output = result;

    return output / 2. / M_PI / M_PI;
}


/**
    renormalized integrand of I^4_0 at r = 0
**/

static double integrals_renormalization0_integrand(
    double k,
#ifdef HAVE_DOUBLE_EXPONENTIAL
    const void *p
#else
    void *p
#endif
)
{
    struct integrals_params *integrand = (struct integrals_params *) p;
    return coffe_interp_spline(integrand->result, k)
        *(1. - pow(gsl_sf_bessel_j0(k * integrand->r), 2))
        / k / k;
}


/**
    integrand of the renormalization term (the divergent one)
**/

static double integrals_renormalization_integrand(
    double k,
#ifdef HAVE_DOUBLE_EXPONENTIAL
    const void *p
#else
    void *p
#endif
)
{
    struct integrals_divergent_params *integrand =
        (struct integrals_divergent_params *) p;
    return coffe_interp_spline(integrand->result, k)
           *gsl_sf_bessel_j0(k * integrand->chi1)
           *gsl_sf_bessel_j0(k * integrand->chi2)
            / k / k;
}


/**
    integrates the above (the general integrator doesn't work here)
**/

static double integrals_renormalization(
    const struct coffe_interpolation *input,
    const double chi1,
    const double chi2,
    const double kmin,
    const double kmax
)
{
    struct integrals_divergent_params test;
    test.result = (struct coffe_interpolation*) input;
    test.chi1 = chi1;
    test.chi2 = chi2;

    const double output = coffe_integrate_1d(
        &integrals_renormalization_integrand,
        &test,
        kmin,
        kmax
    );

    return output / 2. / M_PI / M_PI;
}


/**
    computes the integrals I^n_l(r) for r in [0, infinity]
    If n >= l, returns the combination I^n_l(r) * r^(n - l).
    output_{x,y} are pointers to pointers so we can do memory allocation
    INSIDE the function, thus the user doesn't need to care about it themselves.
    Note that there's two parameters, output_len, and real_output_len;
    the first is just the FFT size, the second is the FFT size + the min separations,
    and should be used after calling this function instead of the first
**/

int coffe_integrals_renormalizable(
    double **output_x,
    double **output_y,
    const size_t output_len,
    size_t *real_output_len,
    const struct coffe_interpolation *spectrum,
    const int n,
    const int l,
    const enum coffe_integer_state state_n,
    const enum coffe_integer_state state_l,
    const double x_min,
    const double x_max
)
{
    /* how big the real_output_len is */
    *real_output_len = output_len + min_sep_len + 1;

    /* memory allocation of output_{x,y} */
    *output_x = (double *)coffe_malloc(sizeof(double) * *real_output_len);
    *output_y = (double *)coffe_malloc(sizeof(double) * *real_output_len);

    /* memory alloc of FFT */
    double *fft_x = (double *)coffe_malloc(sizeof(double) * output_len);
    double *fft_y = (double *)coffe_malloc(sizeof(double) * output_len);

    double local_n, local_l;

    integrals_check_parameters(
        n,
        l,
        state_n,
        state_l,
        &local_n,
        &local_l
    );

    /* do the FFTlog transform first */
    twofast_1bessel(
        fft_x,
        fft_y,
        output_len,
        spectrum->spline->x,
        spectrum->spline->y,
        spectrum->spline->size,
        local_l,
        local_n,
        COFFE_H0,
        x_min,
        x_min,
        x_max,
        0
    );

    double x0_result;

    if (local_n >= local_l){
        /* r^(n - l) * I^n_l(r) */
        coffe_multiply_power_array(
            fft_y,
            fft_y,
            fft_x,
            output_len,
            local_n - local_l
        );

        /* the result I^n_l(0) */
        x0_result = integrals_integrate_function(
                &integrals_prefactor,
                spectrum,
                n,
                l,
                state_n,
                state_l,
                0,
                x_min,
                x_max
            )
           *integrals_coefficients(
                l,
                state_l
            );
    }
    else{
        x0_result = 0.0;
    }

    (*output_x)[0] = 0.0;
    (*output_y)[0] = x0_result;

    /* setting the output for the smallest separations (under 1 Mpc/h) */
    for (size_t i = 1; i <= min_sep_len; ++i){
        (*output_x)[i] = min_sep[i - 1]; // dimensionless!
        (*output_y)[i] = integrals_integrate_function(
                &integrals_bessel_integrand,
                spectrum,
                n,
                l,
                state_n,
                state_l,
                (*output_x)[i],
                x_min,
                x_max
            );
    }

    for (size_t i = *real_output_len - output_len; i < *real_output_len; ++i){
        (*output_x)[i] = fft_x[i - (*real_output_len - output_len)];
        (*output_y)[i] = fft_y[i - (*real_output_len - output_len)];
    }

    free(fft_x);
    free(fft_y);

    return EXIT_SUCCESS;
}


/**
    computes all the nonzero I^n_l integrals
**/

int coffe_integrals_init(
    const struct coffe_parameters_t *par,
    const struct coffe_background_t *bg,
    struct coffe_integral_array_t *integral
)
{
    clock_t start, end;
    start = clock();

    if (par->verbose)
        printf("Calculating integrals of Bessel functions...\n");

    gsl_error_handler_t *default_handler =
        gsl_set_error_handler_off();

    /* default size is zero, and points to nothing */
    integral->array = NULL;
    integral->size = 0;

    if (
        par->output_type == CORRFUNC ||
        par->output_type == MULTIPOLES ||
        par->output_type == AVERAGE_MULTIPOLES
    ){
        /* the default renormalizable linear theory integrals in full-sky */
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

        integral->array = (struct coffe_integral_t *)coffe_malloc(sizeof(struct coffe_integral_t) * sizeof(terms) / sizeof(*terms));

        /* those default renormalized integrals */
        for (size_t i = 0; i < sizeof(terms) / sizeof(*terms); ++i){

            const int n = terms[i].n;
            const int l = terms[i].l;

            const size_t current_index = integral->size;
            /* alloc the space */
            if (current_index == 0)
                integral->array = (struct coffe_integral_t *)coffe_malloc(sizeof(struct coffe_integral_t));
            else
                integral->array = (struct coffe_integral_t *)realloc(
                    integral->array,
                    sizeof(struct coffe_integral_t) * (current_index + 1)
                );
            integral->array[current_index].n = n;
            integral->array[current_index].l = l;
            integral->array[current_index].state_n = COFFE_INTEGER;
            integral->array[current_index].state_l = COFFE_INTEGER;
            integral->array[current_index].renormalization.spline = NULL;
            integral->array[current_index].renormalization.xaccel = NULL;
            integral->array[current_index].renormalization.yaccel = NULL;

            /* if integral is not divergent, use implementation of 2FAST */
            const size_t npoints = par->bessel_bins;
            double *final_sep = NULL;
            double *final_result = NULL;
            size_t output_real_len = 0;

            coffe_integrals_renormalizable(
                &final_sep,
                &final_result,
                npoints,
                &output_real_len,
                &par->power_spectrum_norm,
                integral->array[current_index].n,
                integral->array[current_index].l,
                integral->array[current_index].state_n,
                integral->array[current_index].state_l,
                par->k_min_norm,
                par->k_max_norm
            );

            coffe_init_spline(
                &integral->array[current_index].result,
                final_sep,
                final_result,
                output_real_len,
                par->interp_method
            );

            free(final_sep);
            free(final_result);

#ifdef HAVE_CLASS
            if (
                par->pk_type &&
                par->midpoint_approximation
            ){

                /* do we want the linear, or the nonlinear one? */
                enum pk_outputs pk_type = pk_linear;
                if (par->pk_type == 2 || par->pk_type == 3) pk_type = pk_nonlinear;

                /* TODO make the upper bound modular */
                const double z_min = 0.0, z_max = 3.0;
                /* TODO make this modular */
                const size_t z_size = 100;

                /* list of all the redshifts */
                double *z = (double *)coffe_malloc(sizeof(double) * z_size);

                for (size_t i = 0; i < z_size; ++i)
                    z[i] = z_min + z_max * (double)i / z_size;

                /* alloc memory for 2D interpolation */
                double *pk_at_z2d = (double *)coffe_malloc(sizeof(double) * z_size * output_real_len);

                for (size_t i = 0; i < z_size; ++i){
                    const size_t k_size = ((struct nonlinear *)par->class_nonlinear)->k_size;
                    /* alloc memory for k and pk */
                    double *k = (double *)coffe_malloc(sizeof(double) * k_size);
                    double *pk = (double *)coffe_malloc(sizeof(double) * k_size);

                    /* get (non)linear power spectrum at redshift z (and store it in pk) */
                    nonlinear_pk_at_z(
                        (struct background *)par->class_background,
                        (struct nonlinear *)par->class_nonlinear,
                        logarithmic,
                        pk_type,
                        z[i],
                        ((struct nonlinear *)par->class_nonlinear)->index_pk_total,
                        pk,
                        NULL
                    );

                    /* need to rescale since CLASS internally works in units of 1/Mpc */
                    /* NOTE k and pk are the DIMENSIONLESS spectra (i.e. in units COFFE_H0) */
                    for (size_t j = 0; j < k_size; ++j){
                        k[j] = ((struct nonlinear *)par->class_nonlinear)->k[j] / par->h / COFFE_H0;
                        pk[j] = exp(pk[j]) * pow(par->h, 3) * pow(COFFE_H0, 3);
                    }

                    /* setup the interpolation of the power spectrum */
                    struct coffe_interpolation pk_at_z;
                    coffe_init_spline(
                        &pk_at_z,
                        k,
                        pk,
                        k_size,
                        par->interp_method
                    );

                    /* memory cleanup */
                    free(k);
                    free(pk);

                    /* get the FFTlog of the power spectra at redshift z (I^n_ell(r, z)) */
                    coffe_integrals_renormalizable(
                        &final_sep,
                        &final_result,
                        npoints,
                        &output_real_len,
                        &pk_at_z,
                        integral->array[current_index].n,
                        integral->array[current_index].l,
                        integral->array[current_index].state_n,
                        integral->array[current_index].state_l,
                        par->k_min_norm,
                        par->k_max_norm
                    );

                    /* save everything into a big array */
                    for (size_t j = 0; j < output_real_len; ++j)
                        /* first index redshift, second separation */
                        pk_at_z2d[j * z_size + i] = final_result[j];

                    /* memory cleanup */
                    if (i != z_size - 1) free(final_sep);
                    free(final_result);

                }

                /* now setup the 2D interpolation */
                coffe_init_spline2d(
                    &integral->array[current_index].result2d,
                    z,
                    final_sep,
                    pk_at_z2d,
                    z_size,
                    output_real_len,
                    2
                );

                /* we didn't free it above, so we do it here */
                free(final_sep);

            }
#endif

            integral->size += 1;
        }

        /* first see if we need to handle the divergent integral */
        if (par->divergent){
            const size_t current_index = integral->size;
            /* alloc the space */
            if (current_index == 0)
                integral->array = (struct coffe_integral_t *)coffe_malloc(sizeof(struct coffe_integral_t));
            else
                integral->array = (struct coffe_integral_t *)realloc(
                    integral->array,
                    sizeof(struct coffe_integral_t) * (current_index + 1)
                );
            integral->array[current_index].n = 4;
            integral->array[current_index].l = 0;
            integral->array[current_index].state_n = COFFE_INTEGER;
            integral->array[current_index].state_l = COFFE_INTEGER;

            double *result =
                (double *)coffe_malloc(sizeof(double) * coffe_sep_len);
            double *separations =
                (double *)coffe_malloc(sizeof(double) * coffe_sep_len);
            double *result0 =
                (double *)coffe_malloc(sizeof(double) * coffe_sep_len);

            #pragma omp parallel for num_threads(par->nthreads)
            for (size_t i = 0; i < coffe_sep_len; ++i){
                /* dimensionless */
                separations[i] = NORM(coffe_sep[i]);

                result[i] = integrals_integrate_function(
                    &integrals_bessel_integrand,
                    &par->power_spectrum_norm,
                    4,
                    0,
                    COFFE_INTEGER,
                    COFFE_INTEGER,
                    separations[i],
                    par->k_min_norm,
                    par->k_max_norm
                );

                result0[i] = integrals_integrate_function(
                    &integrals_renormalization0_integrand,
                    &par->power_spectrum_norm,
                    /* n and l don't matter in this case */
                    4,
                    0,
                    COFFE_INTEGER,
                    COFFE_INTEGER,
                    separations[i],
                    par->k_min_norm,
                    par->k_max_norm
                );

            }

            coffe_init_spline(
                &integral->array[current_index].result,
                separations,
                result,
                coffe_sep_len,
                par->interp_method
            );

            separations[0] = 0.0;
            result0[0] = 0.0;

            coffe_init_spline(
                &integral->array[current_index].renormalization_zero_separation,
                separations,
                result0,
                coffe_sep_len,
                par->interp_method
            );
            free(separations);
            free(result);
            free(result0);

            const size_t nbins = 200;
            double *result2d = (double *)coffe_malloc(
                sizeof(double) * (nbins + 1) * (nbins + 1)
            );

            double chi_min = 0.;
            double chi_max = 0.;
            if (
                par->output_type == CORRFUNC ||
                par->output_type == MULTIPOLES
            ){
                /* dimensionless */
                double *temp = (double *)coffe_malloc(sizeof(double) * par->z_mean_len);
                for (size_t i = 0; i < par->z_mean_len; ++i)
                    temp[i] = par->z_mean[i] + par->deltaz[i];
                chi_max = coffe_interp_spline(
                    &bg->comoving_distance,
                    coffe_max_array_double(temp, par->z_mean_len)
                );
                free(temp);
            }
            else if (par->output_type == AVERAGE_MULTIPOLES){
                /* dimensionless */
                chi_max = coffe_interp_spline(
                    &bg->comoving_distance,
                    par->z_max
                );
            }
            double *chi_array = (double *)coffe_malloc(
                sizeof(double) * (nbins + 1)
            );
            for (size_t j = 0; j <= nbins; ++j){
                chi_array[j] = chi_min + (double)j / nbins * (chi_max - chi_min);
            }

            #pragma omp parallel for num_threads(par->nthreads) collapse(2)
            for (size_t j = 0; j <= nbins; ++j){
                for (size_t k = 0; k <= nbins; ++k){
                    result2d[k * (nbins + 1) + j] = integrals_renormalization(
                        &par->power_spectrum_norm,
                        chi_array[j],
                        chi_array[k],
                        par->k_min_norm,
                        par->k_max_norm
                    );
                }
            }

            coffe_init_spline2d(
                &integral->array[current_index].renormalization,
                chi_array,
                chi_array,
                result2d,
                nbins + 1,
                nbins + 1,
                2
            );

            free(chi_array);
            free(result2d);

            /* don't forget to increase the size of the container! */
            integral->size += 1;
        }

        if (
            par->flatsky_local_nonlocal ||
            par->flatsky_nonlocal
        ){
            const size_t current_index = integral->size;
            /* alloc the space */
            if (current_index == 0)
                integral->array = (struct coffe_integral_t *)coffe_malloc(sizeof(struct coffe_integral_t));
            else
                integral->array = (struct coffe_integral_t *)realloc(
                    integral->array,
                    sizeof(struct coffe_integral_t) * (current_index + 1)
                );
            integral->array[current_index].n = 1;
            integral->array[current_index].l = -1;
            integral->array[current_index].state_n = COFFE_HALF_INTEGER;
            integral->array[current_index].state_l = COFFE_HALF_INTEGER;
            integral->array[current_index].renormalization.spline = NULL;
            integral->array[current_index].renormalization.xaccel = NULL;
            integral->array[current_index].renormalization.yaccel = NULL;

            /* if integral is not divergent, use implementation of 2FAST */
            const size_t npoints = par->bessel_bins;
            double *final_sep = NULL;
            double *final_result = NULL;
            size_t output_real_len = 0;

            coffe_integrals_renormalizable(
                &final_sep,
                &final_result,
                npoints,
                &output_real_len,
                &par->power_spectrum_norm,
                integral->array[current_index].n,
                integral->array[current_index].l,
                integral->array[current_index].state_n,
                integral->array[current_index].state_l,
                par->k_min_norm,
                par->k_max_norm
            );

            /* we're off by a factor sqrt(2 / pi) * 2 * pi^2, so we rescale first */
            /* note that we're actually computing r * I(r) since that's well defined at r = 0 */
            coffe_rescale_array(
                final_result,
                output_real_len,
                sqrt(2. / M_PI) * 2 * M_PI * M_PI
            );

            coffe_init_spline(
                &integral->array[current_index].result,
                final_sep,
                final_result,
                output_real_len,
                par->interp_method
            );

            free(final_sep);
            free(final_result);

#ifdef HAVE_CLASS
            if (
                par->pk_type &&
                par->midpoint_approximation
            ){

                /* do we want the linear, or the nonlinear one? */
                enum pk_outputs pk_type = pk_linear;
                if (par->pk_type == 2 || par->pk_type == 3) pk_type = pk_nonlinear;

                /* TODO make the upper bound modular */
                const double z_min = 0.0, z_max = 3.0;
                /* TODO make this modular */
                const size_t z_size = 100;

                /* list of all the redshifts */
                double *z = (double *)coffe_malloc(sizeof(double) * z_size);

                for (size_t i = 0; i < z_size; ++i)
                    z[i] = z_min + z_max * (double)i / z_size;

                /* alloc memory for 2D interpolation */
                double *pk_at_z2d = (double *)coffe_malloc(sizeof(double) * z_size * output_real_len);

                for (size_t i = 0; i < z_size; ++i){
                    const size_t k_size = ((struct nonlinear *)par->class_nonlinear)->k_size;
                    /* alloc memory for k and pk */
                    double *k = (double *)coffe_malloc(sizeof(double) * k_size);
                    double *pk = (double *)coffe_malloc(sizeof(double) * k_size);

                    /* get (non)linear power spectrum at redshift z (and store it in pk) */
                    nonlinear_pk_at_z(
                        (struct background *)par->class_background,
                        (struct nonlinear *)par->class_nonlinear,
                        logarithmic,
                        pk_type,
                        z[i],
                        ((struct nonlinear *)par->class_nonlinear)->index_pk_total,
                        pk,
                        NULL
                    );

                    /* need to rescale since CLASS internally works in units of 1/Mpc */
                    /* NOTE k and pk are the DIMENSIONLESS spectra (i.e. in units COFFE_H0) */
                    for (size_t j = 0; j < k_size; ++j){
                        k[j] = ((struct nonlinear *)par->class_nonlinear)->k[j] / par->h / COFFE_H0;
                        pk[j] = exp(pk[j]) * pow(par->h, 3) * pow(COFFE_H0, 3);
                    }

                    /* setup the interpolation of the power spectrum */
                    struct coffe_interpolation pk_at_z;
                    coffe_init_spline(
                        &pk_at_z,
                        k,
                        pk,
                        k_size,
                        par->interp_method
                    );

                    /* memory cleanup */
                    free(k);
                    free(pk);

                    /* get the FFTlog of the power spectra at redshift z (I^n_ell(r, z)) */
                    coffe_integrals_renormalizable(
                        &final_sep,
                        &final_result,
                        npoints,
                        &output_real_len,
                        &pk_at_z,
                        integral->array[current_index].n,
                        integral->array[current_index].l,
                        integral->array[current_index].state_n,
                        integral->array[current_index].state_l,
                        par->k_min_norm,
                        par->k_max_norm
                    );

                    /* save everything into a big array */
                    for (size_t j = 0; j < output_real_len; ++j)
                        /* first index redshift, second separation */
                        pk_at_z2d[j * z_size + i] = final_result[j];

                    /* memory cleanup */
                    if (i != z_size - 1) free(final_sep);
                    free(final_result);

                }

                /* now setup the 2D interpolation */
                coffe_init_spline2d(
                    &integral->array[current_index].result2d,
                    z,
                    final_sep,
                    pk_at_z2d,
                    z_size,
                    output_real_len,
                    2
                );

                /* we didn't free it above, so we do it here */
                free(final_sep);

            }
#endif


            integral->size += 1;
        }
        /* lensing-lensing multipoles are special */
        if (
            par->flatsky_nonlocal &&
            par->output_type == MULTIPOLES
        ){
            for (size_t i = 0; i < par->multipole_values_len; ++i){
                const int l = par->multipole_values[i];

                const size_t current_index = integral->size;
                /* alloc the space */
                if (current_index == 0)
                    integral->array = (struct coffe_integral_t *)coffe_malloc(sizeof(struct coffe_integral_t));
                else
                    integral->array = (struct coffe_integral_t *)realloc(
                        integral->array,
                        sizeof(struct coffe_integral_t) * (current_index + 1)
                    );
                integral->array[current_index].n = 1;
                integral->array[current_index].l = l;
                integral->array[current_index].state_n = COFFE_INTEGER;
                integral->array[current_index].state_l = COFFE_INTEGER;
                integral->array[current_index].renormalization.spline = NULL;
                integral->array[current_index].renormalization.xaccel = NULL;
                integral->array[current_index].renormalization.yaccel = NULL;

                /* if integral is not divergent, use implementation of 2FAST */
                const size_t npoints = par->bessel_bins;
                double *final_sep = NULL;
                double *final_result = NULL;
                size_t output_real_len = 0;

                coffe_integrals_renormalizable(
                    &final_sep,
                    &final_result,
                    npoints,
                    &output_real_len,
                    &par->power_spectrum_norm,
                    integral->array[current_index].n,
                    integral->array[current_index].l,
                    integral->array[current_index].state_n,
                    integral->array[current_index].state_l,
                    par->k_min_norm,
                    par->k_max_norm
                );

                coffe_init_spline(
                    &integral->array[current_index].result,
                    final_sep,
                    final_result,
                    output_real_len,
                    par->interp_method
                );

                free(final_sep);
                free(final_result);

#ifdef HAVE_CLASS
            if (
                par->pk_type &&
                par->midpoint_approximation
            ){

                /* do we want the linear, or the nonlinear one? */
                enum pk_outputs pk_type = pk_linear;
                if (par->pk_type == 2 || par->pk_type == 3) pk_type = pk_nonlinear;

                /* TODO make the upper bound modular */
                const double z_min = 0.0, z_max = 3.0;
                /* TODO make this modular */
                const size_t z_size = 100;

                /* list of all the redshifts */
                double *z = (double *)coffe_malloc(sizeof(double) * z_size);

                for (size_t i = 0; i < z_size; ++i)
                    z[i] = z_min + z_max * (double)i / z_size;

                /* alloc memory for 2D interpolation */
                double *pk_at_z2d = (double *)coffe_malloc(sizeof(double) * z_size * output_real_len);

                for (size_t i = 0; i < z_size; ++i){
                    const size_t k_size = ((struct nonlinear *)par->class_nonlinear)->k_size;
                    /* alloc memory for k and pk */
                    double *k = (double *)coffe_malloc(sizeof(double) * k_size);
                    double *pk = (double *)coffe_malloc(sizeof(double) * k_size);

                    /* get (non)linear power spectrum at redshift z (and store it in pk) */
                    nonlinear_pk_at_z(
                        (struct background *)par->class_background,
                        (struct nonlinear *)par->class_nonlinear,
                        logarithmic,
                        pk_type,
                        z[i],
                        ((struct nonlinear *)par->class_nonlinear)->index_pk_total,
                        pk,
                        NULL
                    );

                    /* need to rescale since CLASS internally works in units of 1/Mpc */
                    /* NOTE k and pk are the DIMENSIONLESS spectra (i.e. in units COFFE_H0) */
                    for (size_t j = 0; j < k_size; ++j){
                        k[j] = ((struct nonlinear *)par->class_nonlinear)->k[j] / par->h / COFFE_H0;
                        pk[j] = exp(pk[j]) * pow(par->h, 3) * pow(COFFE_H0, 3);
                    }

                    /* setup the interpolation of the power spectrum */
                    struct coffe_interpolation pk_at_z;
                    coffe_init_spline(
                        &pk_at_z,
                        k,
                        pk,
                        k_size,
                        par->interp_method
                    );

                    /* memory cleanup */
                    free(k);
                    free(pk);

                    /* get the FFTlog of the power spectra at redshift z (I^n_ell(r, z)) */
                    coffe_integrals_renormalizable(
                        &final_sep,
                        &final_result,
                        npoints,
                        &output_real_len,
                        &pk_at_z,
                        integral->array[current_index].n,
                        integral->array[current_index].l,
                        integral->array[current_index].state_n,
                        integral->array[current_index].state_l,
                        par->k_min_norm,
                        par->k_max_norm
                    );

                    /* save everything into a big array */
                    for (size_t j = 0; j < output_real_len; ++j)
                        /* first index redshift, second separation */
                        pk_at_z2d[j * z_size + i] = final_result[j];

                    /* memory cleanup */
                    if (i != z_size - 1) free(final_sep);
                    free(final_result);

                }

                /* now setup the 2D interpolation */
                coffe_init_spline2d(
                    &integral->array[current_index].result2d,
                    z,
                    final_sep,
                    pk_at_z2d,
                    z_size,
                    output_real_len,
                    2
                );

                /* we didn't free it above, so we do it here */
                free(final_sep);

            }
#endif

                integral->size += 1;
            }
        }
        /* density-lensing multipoles are also special */
        if (
            par->flatsky_local_nonlocal &&
            par->output_type == MULTIPOLES
        ){
            for (size_t i = 0; i < par->multipole_values_len; ++i){
                const int l = par->multipole_values[i];
                for (size_t k = 0; k <= (size_t)l / 2; ++k){

                    const size_t current_index = integral->size;
                    /* alloc the space */
                    if (current_index == 0)
                        integral->array = (struct coffe_integral_t *)coffe_malloc(sizeof(struct coffe_integral_t));
                    else
                        integral->array = (struct coffe_integral_t *)realloc(
                            integral->array,
                            sizeof(struct coffe_integral_t) * (current_index + 1)
                        );
                    integral->array[current_index].n = l - 2 * k + 3;
                    integral->array[current_index].l = l - 2 * k + 1;
                    integral->array[current_index].state_n = COFFE_HALF_INTEGER;
                    integral->array[current_index].state_l = COFFE_HALF_INTEGER;
                    integral->array[current_index].renormalization.spline = NULL;
                    integral->array[current_index].renormalization.xaccel = NULL;
                    integral->array[current_index].renormalization.yaccel = NULL;

                    /* if integral is not divergent, use implementation of 2FAST */
                    const size_t npoints = par->bessel_bins;
                    double *final_sep = NULL;
                    double *final_result = NULL;
                    size_t output_real_len = 0;

                    coffe_integrals_renormalizable(
                        &final_sep,
                        &final_result,
                        npoints,
                        &output_real_len,
                        &par->power_spectrum_norm,
                        integral->array[current_index].n,
                        integral->array[current_index].l,
                        integral->array[current_index].state_n,
                        integral->array[current_index].state_l,
                        par->k_min_norm,
                        par->k_max_norm
                    );

                    coffe_init_spline(
                        &integral->array[current_index].result,
                        final_sep,
                        final_result,
                        output_real_len,
                        par->interp_method
                    );

                    free(final_sep);
                    free(final_result);

#ifdef HAVE_CLASS
            if (
                par->pk_type &&
                par->midpoint_approximation
            ){

                /* do we want the linear, or the nonlinear one? */
                enum pk_outputs pk_type = pk_linear;
                if (par->pk_type == 2 || par->pk_type == 3) pk_type = pk_nonlinear;

                /* TODO make the upper bound modular */
                const double z_min = 0.0, z_max = 3.0;
                /* TODO make this modular */
                const size_t z_size = 100;

                /* list of all the redshifts */
                double *z = (double *)coffe_malloc(sizeof(double) * z_size);

                for (size_t i = 0; i < z_size; ++i)
                    z[i] = z_min + z_max * (double)i / z_size;

                /* alloc memory for 2D interpolation */
                double *pk_at_z2d = (double *)coffe_malloc(sizeof(double) * z_size * output_real_len);

                for (size_t i = 0; i < z_size; ++i){
                    const size_t k_size = ((struct nonlinear *)par->class_nonlinear)->k_size;
                    /* alloc memory for k and pk */
                    double *k = (double *)coffe_malloc(sizeof(double) * k_size);
                    double *pk = (double *)coffe_malloc(sizeof(double) * k_size);

                    /* get (non)linear power spectrum at redshift z (and store it in pk) */
                    nonlinear_pk_at_z(
                        (struct background *)par->class_background,
                        (struct nonlinear *)par->class_nonlinear,
                        logarithmic,
                        pk_type,
                        z[i],
                        ((struct nonlinear *)par->class_nonlinear)->index_pk_total,
                        pk,
                        NULL
                    );

                    /* need to rescale since CLASS internally works in units of 1/Mpc */
                    /* NOTE k and pk are the DIMENSIONLESS spectra (i.e. in units COFFE_H0) */
                    for (size_t j = 0; j < k_size; ++j){
                        k[j] = ((struct nonlinear *)par->class_nonlinear)->k[j] / par->h / COFFE_H0;
                        pk[j] = exp(pk[j]) * pow(par->h, 3) * pow(COFFE_H0, 3);
                    }

                    /* setup the interpolation of the power spectrum */
                    struct coffe_interpolation pk_at_z;
                    coffe_init_spline(
                        &pk_at_z,
                        k,
                        pk,
                        k_size,
                        par->interp_method
                    );

                    /* memory cleanup */
                    free(k);
                    free(pk);

                    /* get the FFTlog of the power spectra at redshift z (I^n_ell(r, z)) */
                    coffe_integrals_renormalizable(
                        &final_sep,
                        &final_result,
                        npoints,
                        &output_real_len,
                        &pk_at_z,
                        integral->array[current_index].n,
                        integral->array[current_index].l,
                        integral->array[current_index].state_n,
                        integral->array[current_index].state_l,
                        par->k_min_norm,
                        par->k_max_norm
                    );

                    /* save everything into a big array */
                    for (size_t j = 0; j < output_real_len; ++j)
                        /* first index redshift, second separation */
                        pk_at_z2d[j * z_size + i] = final_result[j];

                    /* memory cleanup */
                    if (i != z_size - 1) free(final_sep);
                    free(final_result);

                }

                /* now setup the 2D interpolation */
                coffe_init_spline2d(
                    &integral->array[current_index].result2d,
                    z,
                    final_sep,
                    pk_at_z2d,
                    z_size,
                    output_real_len,
                    2
                );

                /* we didn't free it above, so we do it here */
                free(final_sep);

            }
#endif

                    integral->size += 1;
                }
            }
        }
    }

    gsl_set_error_handler(default_handler);
    end = clock();

    if (par->verbose)
        printf("Integrals of Bessel functions calculated in %.2f s\n",
            (double)(end - start) / CLOCKS_PER_SEC);

    return EXIT_SUCCESS;
}


int coffe_integrals_free(
    struct coffe_integral_array_t *integral
)
{
    return EXIT_SUCCESS;
    if (integral->size > 0){
        for (size_t i = 0; i < integral->size; ++i){
            coffe_free_spline(&integral->array[i].result);
            coffe_free_spline(&integral->array[i].renormalization_zero_separation);
            coffe_free_spline2d(&integral->array[i].renormalization);
        }
        /* so we don't free it again */
        integral->size = 0;
        free(integral->array);
    }
    return EXIT_SUCCESS;
}

#undef NORM
