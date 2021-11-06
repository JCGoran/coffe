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

#include <math.h>
#include <time.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include "common.h"
#include "background.h"
#include "covariance.h"
#include "twobessel.h"

#ifdef HAVE_CLASS
#include "class.h"
#endif


/**
    find the element in the array which corresponds to some values of the parameters (z_mean, l1, l2, sep1, sep2)
    TODO what if we can't find anything?
**/

coffe_covariance_t coffe_covariance_find(
    const coffe_covariance_array_t *cov,
    const double z_mean,
    const int l1,
    const int l2,
    const double sep1,
    const double sep2
)
{
    for (size_t i = 0; i < cov->size; ++i){
        if (
            cov->array[i].coords.l1 == l1 &&
            cov->array[i].coords.l2 == l2 &&
            coffe_approx_equal(
                cov->array[i].coords.separation1, sep1, 1e-3, 0
            ) &&
            coffe_approx_equal(
                cov->array[i].coords.separation2, sep2, 1e-3, 0
            ) &&
            coffe_approx_equal(
                cov->array[i].coords.z_mean, z_mean, 1e-3, 0
            )
        ){
            return cov->array[i];
        }
    }
    fprintf(stderr, "ERROR: cannot find the appropriate element in the covariance array.\n");
    exit(EXIT_FAILURE);
}



/**
    the struct for integration of the windowed spherical Bessel function
**/
struct covariance_params_windowed
{
    int ell;
    double k;
    double s_mean, s_delta; /* in Mpc/h */
};


/**
    contains the necessary integration parameters
**/
struct covariance_params
{
    struct coffe_interpolation *power_spectrum;
    double chi1, chi2;
    int l1, l2;
    struct covariance_params_windowed windowed1, windowed2;
};


/**
    contains the parameter necessary to calculate the volume for average multipoles
**/
struct covariance_volume_params
{
    struct coffe_interpolation *conformal_Hz;
    struct coffe_interpolation *comoving_distance;
};


/**
    integrand of the volume
**/
static double covariance_volume_integrand(
    double z,
    void *p
)
{
    struct covariance_volume_params *s =
        (struct covariance_volume_params *) p;
    double result = 1.
       /coffe_interp_spline(s->conformal_Hz, z)
       /pow(coffe_interp_spline(s->comoving_distance, z), 2)
       /(1 + z);
    return result;
}


/**
    helper function for computing the volume element (without the 4 pi)
**/
static inline double covariance_volume_no_4pi(
    const double r_min,
    const double r_max
)
{
    return (pow(r_max, 3) - pow(r_min, 3)) / 3.;
}


/**
    calculates i^(l1 - l2) when l1 - l2 is even
**/
static int covariance_complex(int l1, int l2)
{
    if ((l1 - l2) % 4 == 0) return 1;
    else if ((l1 - l2) % 2 == 0) return -1;
    else return 0;
}


/**
    integrand of the windowed spherical Bessel function
**/
static double covariance_windowed_bessel_integrand(
    const double s,
    void *p
)
{
    struct covariance_params_windowed *test =
        (struct covariance_params_windowed *) p;
    return s * s * gsl_sf_bessel_jl(test->ell, s * test->k);
}


/**
    analytic value of the integral of the above for k = 1.
    Obtained with Wolfram Mathematica 11.1.1.0 on Linux x86_64
    Generalizations described in arXiv:1703.06428
**/
static double integral_power_spherical_bessel1_unit(
    const int n,
    const int l,
    const double x
){
    if (n != 2){
        fprintf(stderr, "WARNING: n = %d not implemented\n", n);
        return 0;
    }
    switch (l){
        case 0:
            return -x * cos(x) + sin(x);
            break;
        case 1:
            return - 2. * cos(x) - x * sin(x);
            break;
        case 2:
            return x * cos(x) - 4. * sin(x) + 3. * gsl_sf_Si(x);
            break;
        case 3:
            return 7. * cos(x) + ((-15. + x * x) * sin(x)) / x;
            break;
        case 4:
            return ((105. / x - 2. * x) * cos(x) + (22. - 105. / pow(x, 2)) * sin(x) + 15. * gsl_sf_Si(x)) / 2.;
            break;
        default:
            fprintf(stderr, "WARNING: l = %d not implemented.\n", l);
            break;
    }
    return 0;
}


/**
    -||-, but with k != 1
**/
static double integral_power_spherical_bessel1(
    const int n,
    const int l,
    const double k,
    const double x
)
{
    return
        integral_power_spherical_bessel1_unit(n, l, k * x) / pow(k, n + 1);
}


/**
    windowed spherical Bessel function
**/
static double covariance_windowed_bessel(
    const double s,
    const struct covariance_params_windowed *test
)
{
    /* if it's zero, there's nothing to window, and we just return the
    spherical Bessel function itself */
    if (test->s_delta == 0){
        return gsl_sf_bessel_jl(test->ell, s * test->k);
    }
    return 1. / covariance_volume_no_4pi(
        test->s_mean - test->s_delta / 2,
        test->s_mean + test->s_delta / 2
        )
       *(
            integral_power_spherical_bessel1(
                2, test->ell, test->k,
                test->s_mean + test->s_delta / 2.
            )
            -
            integral_power_spherical_bessel1(
                2, test->ell, test->k,
                test->s_mean - test->s_delta / 2.
            )
        );
}


/**
    integrand P(k) (or P(k)^2) * k^2 * j_l1(k\chi_1) * j_l2(k\chi_2)
**/
static double covariance_integrand(
    double k,
    void *p
)
{
    struct covariance_params *test = (struct covariance_params *) p;

    test->windowed1.k = k;
    test->windowed1.ell = test->l1;
    test->windowed1.s_mean = test->chi1;

    test->windowed2.k = k;
    test->windowed2.ell = test->l2;
    test->windowed2.s_mean = test->chi2;

    return
        k * k
       *coffe_interp_spline(test->power_spectrum, k)
       *covariance_windowed_bessel(test->chi1, &test->windowed1)
       *covariance_windowed_bessel(test->chi2, &test->windowed2);
}


/**
    integrates the above
**/
static double covariance_integral(
    struct coffe_interpolation *power_spectrum,
    double chi1, double chi2,
    double chi1_delta, double chi2_delta,
    int l1, int l2,
    double kmin, double kmax
)
{
    struct covariance_params test;
    test.power_spectrum = power_spectrum;
    test.chi1 = chi1;
    test.chi2 = chi2;
    test.l1 = l1;
    test.l2 = l2;
    test.windowed1.s_delta = chi1_delta;
    test.windowed2.s_delta = chi2_delta;

    return 2 * coffe_integrate_1d(
        &covariance_integrand,
        &test,
        kmin,
        kmax
    ) / M_PI;
}


/**
    2d FFTlog integration
**/
static int covariance_integrate_fftlog(
    const struct coffe_interpolation *spline,
    const double power,
    const int l1,
    const int l2,
    const size_t sampling_points,
    const double *separations,
    const size_t npixels_max,
    const double k_min_norm,
    const double k_max_norm,
    const int interpolation_method,
    double *result
)
{

    /* setting the config */
    config covariance_config;
    covariance_config.nu1 = 1.01;
    covariance_config.nu2 = 1.01;
    covariance_config.c_window_width = 0.25;
    covariance_config.l1 = l1;
    covariance_config.l2 = l2;

    /* step size */
    const double step = log(
        k_max_norm / k_min_norm
    ) / sampling_points;

    /* first point of order: logarithmically sample k and P(k) */
    double *k_sampled, **pk_sampled;
    /* memory alloc */
    k_sampled = (double *)coffe_malloc(sizeof(double) * sampling_points);
    pk_sampled = (double **)coffe_malloc(sizeof(double *) * sampling_points);

    for (size_t i = 0; i < sampling_points; ++i){
        /* k sampled in log space */
        k_sampled[i] = k_min_norm * pow(k_max_norm / k_min_norm, (double)i / sampling_points);
        pk_sampled[i] = (double *)coffe_malloc(sizeof(double) * sampling_points);
        /* off-diagonal elements are 0 */
        for (size_t j = 0; j < sampling_points; ++j)
            pk_sampled[i][j] = 0;
        /* diagonal elements equal to k^3 P(k) / dlnk */
        pk_sampled[i][i] =
            pow(k_sampled[i], 3)
           *pow(coffe_interp_spline(spline, k_sampled[i]), power)
           /step;
    }

    /* separations */
    double *r1, *r2;
    r1 = (double *)coffe_malloc(sampling_points * sizeof(double));
    r2 = (double *)coffe_malloc(sampling_points * sizeof(double));

    /* the result as a 2D array */
    double **result_pk = (double **)coffe_malloc(sampling_points * sizeof(double *));
    for(size_t i = 0; i < sampling_points; ++i)
        result_pk[i] = (double *)coffe_malloc(sampling_points * sizeof(double));

    /* integral of P(k) */
    two_sph_bessel(
        k_sampled, k_sampled, pk_sampled,
        sampling_points, sampling_points, &covariance_config,
        r1, r2, result_pk
    );

    /* we don't need x_sampled and y_sampled anymore */
    free(k_sampled);
    for (size_t i = 0; i < sampling_points; ++i)
        free(pk_sampled[i]);
    free(pk_sampled);

    /* realigning 2D array into 1D for GSL interpolation */
    double *result2d_pk = (double *)coffe_malloc(sizeof(double) * sampling_points * sampling_points);

    for (size_t m = 0; m < sampling_points; ++m){
    for (size_t n = 0; n < sampling_points; ++n){
        result2d_pk[n * sampling_points + m] = result_pk[m][n];
    }}

    /* we don't need result_pk anymore */
    for (size_t i = 0; i < sampling_points; ++i)
        free(result_pk[i]);
    free(result_pk);

    struct coffe_interpolation2d pk;

    coffe_init_spline2d(
        &pk,
        r1, r2, result2d_pk,
        sampling_points, sampling_points,
        interpolation_method
    );

    /* we don't need r1, r2, result2d_pk anymore */
    free(r1);
    free(r2);
    free(result2d_pk);

    for (size_t m = 0; m < npixels_max; ++m){
    for (size_t n = 0; n < npixels_max; ++n){
        /* this one is dimensionless */
        result[npixels_max * n + m] =
            coffe_interp_spline2d(
                &pk,
                COFFE_H0 * separations[m],
                COFFE_H0 * separations[n]
            );
    }}

    /* final memory cleanup */
    coffe_free_spline2d(&pk);

    return EXIT_SUCCESS;
}


/**
    computes the covariance of either multipoles or redshift averaged
    multipoles
**/
int coffe_covariance_init(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_covariance_array_t *cov_mp,
    struct coffe_covariance_array_t *cov_ramp
)
{
    if (par->output_type == 4 || par->output_type == 5){
        coffe_covariance_free(cov_mp);
        coffe_covariance_free(cov_ramp);
        time_t start, end;
        start = clock();

        if (par->verbose){
            if (par->output_type == 4)
                printf("Calculating covariance of multipoles...\n");
            else
                printf("Calculating covariance of redshift averaged multipoles...\n");
        }

        gsl_error_handler_t *default_handler =
            gsl_set_error_handler_off();

        if (par->output_type == 4){
            cov_mp->size =
              par->multipole_values_len * par->multipole_values_len
            * par->sep_len * par->sep_len
            * par->covariance_z_mean_len;
            cov_mp->array = (coffe_covariance_t *)coffe_malloc(
                sizeof(coffe_covariance_t) * cov_mp->size
            );
        }
        else{
            cov_ramp->size =
              par->multipole_values_len * par->multipole_values_len
            * par->sep_len * par->sep_len
            * par->covariance_zmin_len;
            cov_ramp->array = (coffe_covariance_t *)coffe_malloc(
                sizeof(coffe_covariance_t) * cov_ramp->size
            );
        }

        double *separations = par->sep;
        const size_t redshifts_to_allocate =
            /* trigraph */
            par->pk_type ? par->covariance_z_mean_len : 1;
        /* allocating memory for the integrals of P(k) and P^2(k) (D_l1l2 and G_l1l2) */
        /* why triple pointers, you ask? First index is redshift, second is multipole, third are separations */
        double ***integral_pk =
            (double ***)coffe_malloc(sizeof(double **) * redshifts_to_allocate);
        double ***integral_pk2 =
            (double ***)coffe_malloc(sizeof(double **) * redshifts_to_allocate);

        for (size_t index_redshift = 0; index_redshift < redshifts_to_allocate; ++index_redshift){
            integral_pk[index_redshift] =
                (double **)coffe_malloc(sizeof(double *) * par->multipole_values_len * par->multipole_values_len);
            integral_pk2[index_redshift] =
                (double **)coffe_malloc(sizeof(double *) * par->multipole_values_len * par->multipole_values_len);
            for (size_t i = 0; i < par->multipole_values_len; ++i){
            for (size_t j = 0; j<par->multipole_values_len; ++j){
                integral_pk[index_redshift][i*par->multipole_values_len + j] =
                    (double *)coffe_malloc(sizeof(double) * par->sep_len * par->sep_len);
                integral_pk2[index_redshift][i*par->multipole_values_len + j] =
                    (double *)coffe_malloc(sizeof(double) * par->sep_len * par->sep_len);
            }}
        }
        for (size_t index_redshift = 0; index_redshift < redshifts_to_allocate; ++index_redshift){
            size_t k_size;
            double k_min, k_max;
            double *k, *pk;
#ifdef HAVE_CLASS
            if (par->have_class && par->pk_type){
                k_size = ((struct nonlinear *)par->class_nonlinear)->k_size;
                /* alloc memory for k and pk */
                k = (double *)coffe_malloc(sizeof(double) * k_size);
                pk = (double *)coffe_malloc(sizeof(double) * k_size);
                enum pk_outputs pk_type = pk_linear;
                if (par->pk_type == 2 || par->pk_type == 3) pk_type = pk_nonlinear;

                /* get (non)linear power spectrum at redshift z (and store it in pk) */
                nonlinear_pk_at_z(
                    (struct background *)par->class_background,
                    (struct nonlinear *)par->class_nonlinear,
                    logarithmic,
                    pk_type,

                    /* trigraph */
                    par->pk_type ? par->covariance_z_mean[index_redshift] : 0,
                    ((struct nonlinear *)par->class_nonlinear)->index_pk_total,
                    pk,
                    NULL
                );

                /* need to rescale since CLASS internally works in units of 1/Mpc */
                for (size_t j = 0; j < k_size; ++j){
                    k[j] = ((struct nonlinear *)par->class_nonlinear)->k[j] / par->h;
                    pk[j] = exp(pk[j]) * pow(par->h, 3);
                }

                k_min = k[0];
                k_max = k[k_size - 1];
            }
            else{
#endif
                k_size = par->power_spectrum.spline->size;
                k_min = par->k_min;
                k_max = par->k_max;
                k = (double *)coffe_malloc(sizeof(double) * par->power_spectrum.spline->size);
                pk = (double *)coffe_malloc(sizeof(double) * par->power_spectrum.spline->size);
                for (size_t j = 0; j < k_size; ++j){
                    k[j] = par->power_spectrum.spline->x[j];
                    pk[j] = par->power_spectrum.spline->y[j];
                }
#ifdef HAVE_CLASS
            }
#endif

            /* setup the interpolation of the power spectrum */
            struct coffe_interpolation pk_at_z;
            coffe_init_spline(
                &pk_at_z,
                k,
                pk,
                k_size,
                par->interp_method
            );

            double *temp_spectrum_pk =
                (double *)coffe_malloc(sizeof(double) * k_size);
            double *temp_spectrum_pk2 =
                (double *)coffe_malloc(sizeof(double) * k_size);
            struct coffe_interpolation integrand_pk, integrand_pk2;

            /* setting the power spectra P(k) and P^2(k) */
            for (size_t i = 0; i < k_size; ++i){
                temp_spectrum_pk[i] =
                    pk_at_z.spline->y[i];
                temp_spectrum_pk2[i] =
                    pow(pk_at_z.spline->y[i], 2);
            }

            /* interpolations of P(k) and P^2(k) */
            coffe_init_spline(
                &integrand_pk,
                k,
                temp_spectrum_pk,
                k_size,
                5
            );

            coffe_init_spline(
                &integrand_pk2,
                k,
                temp_spectrum_pk2,
                k_size,
                5
            );

            /* memory cleanup */
            free(temp_spectrum_pk);
            free(temp_spectrum_pk2);
            free(k);
            free(pk);

            if (par->covariance_integration_method == 1){
                /* calculating the integrals G_l1l2 and D_l1l2 (without the scale factor D1) */
                for (size_t i = 0; i < par->multipole_values_len; ++i){
                for (size_t j = i; j < par->multipole_values_len; ++j){
                    #pragma omp parallel for num_threads(par->nthreads) collapse(2)
                    for (size_t m = 0; m < par->sep_len; ++m){
                    for (size_t n = 0; n < par->sep_len; ++n){
                        integral_pk[index_redshift][i * par->multipole_values_len + j][par->sep_len * n + m] =
                            (2 * par->multipole_values[i] + 1)
                           *(2 * par->multipole_values[j] + 1)
                           *covariance_integral(
                                &integrand_pk,
                                separations[m],
                                separations[n],
                                /* trigraph */
                                (par->covariance_window ? par->covariance_pixelsize[index_redshift] : 0),
                                /* trigraph */
                                (par->covariance_window ? par->covariance_pixelsize[index_redshift] : 0),
                                par->multipole_values[i],
                                par->multipole_values[j],
                                k_min,
                                k_max
                            )
                           /M_PI;
                    }}
                    #pragma omp parallel for num_threads(par->nthreads) collapse(2)
                    for (size_t m = 0; m < par->sep_len; ++m){
                    for (size_t n = 0; n < par->sep_len; ++n){
                        integral_pk2[index_redshift][i * par->multipole_values_len + j][par->sep_len * n + m] =
                            (2 * par->multipole_values[i] + 1)
                           *(2 * par->multipole_values[j] + 1)
                           *covariance_integral(
                                &integrand_pk2,
                                separations[m],
                                separations[n],
                                /* trigraph */
                                (par->covariance_window ? par->covariance_pixelsize[index_redshift] : 0),
                                /* trigraph */
                                (par->covariance_window ? par->covariance_pixelsize[index_redshift] : 0),
                                par->multipole_values[i],
                                par->multipole_values[j],
                                k_min,
                                k_max
                            )
                           /2.
                           /M_PI;
                    }}
                }}
            }
            else if (par->covariance_integration_method == 2){
                fprintf(stderr, "Using 2D FFTlog\n");
                /* here we do the 2DFFTlog */

                /* first point of order: logarithmically sample k and P(k) */
                const size_t sampling_points = par->covariance_integration_bins;

                /* compute the whole thing (just the upper half to save time and space) */
                for (size_t i = 0; i < par->multipole_values_len; ++i){
                for (size_t j = i; j < par->multipole_values_len; ++j){
                    covariance_integrate_fftlog(
                        &par->power_spectrum_norm,
                        1,
                        par->multipole_values[i],
                        par->multipole_values[j],
                        sampling_points,
                        separations,
                        par->sep_len,
                        par->k_min_norm,
                        par->k_max_norm,
                        par->covariance_interpolation_method,
                        integral_pk[index_redshift][i * par->multipole_values_len + j]
                    );

                    covariance_integrate_fftlog(
                        &par->power_spectrum_norm,
                        2,
                        par->multipole_values[i],
                        par->multipole_values[j],
                        sampling_points,
                        separations,
                        par->sep_len,
                        par->k_min_norm,
                        par->k_max_norm,
                        par->covariance_interpolation_method,
                        integral_pk2[index_redshift][i * par->multipole_values_len + j]
                    );

                    for (size_t m = 0; m < par->sep_len; ++m){
                    for (size_t n = 0; n < par->sep_len; ++n){
                        /* this one is dimensionless */
                        integral_pk[index_redshift][i * par->multipole_values_len + j][par->sep_len * n + m] *=
                            (2 * par->multipole_values[i] + 1)
                           *(2 * par->multipole_values[j] + 1)
                           * 2. / M_PI / M_PI;
                        /* this one is not so we need to put a conversion factor */
                        integral_pk2[index_redshift][i * par->multipole_values_len + j][par->sep_len * n + m] *=
                            (2 * par->multipole_values[i] + 1)
                           *(2 * par->multipole_values[j] + 1)
                           / M_PI / M_PI / pow(COFFE_H0, 3);
                    }}
                }}
            }

            /* memory cleanup */
            coffe_free_spline(&integrand_pk);
            coffe_free_spline(&integrand_pk2);
            coffe_free_spline(&pk_at_z);
        }

        /* setting the transpose */
        for (size_t index_redshift = 0; index_redshift < redshifts_to_allocate; ++index_redshift){
        for (size_t i = 0; i < par->multipole_values_len; ++i){
        for (size_t j = 0; j < i; ++j){
        for (size_t m = 0; m < par->sep_len; ++m){
        for (size_t n = 0; n < par->sep_len; ++n){
            integral_pk[index_redshift][i*par->multipole_values_len + j][par->sep_len*n + m] =
                integral_pk[index_redshift][j*par->multipole_values_len + i][par->sep_len*m + n];
            integral_pk2[index_redshift][i*par->multipole_values_len + j][par->sep_len*n + m] =
                integral_pk2[index_redshift][j*par->multipole_values_len + i][par->sep_len*m + n];
        }}}}}

        double *volume =
            (double *)coffe_malloc(sizeof(double)*par->covariance_density_len);
        size_t index = 0;
        double z_mean = 0;
        for (size_t k = 0; k < par->covariance_density_len; ++k){

            if (par->output_type == 4){
                z_mean = par->z_mean[k];
                volume[k] = 4 * M_PI
                   *par->covariance_fsky[k]
                   *covariance_volume_no_4pi(
                        coffe_interp_spline(
                            &bg->comoving_distance,
                            par->z_mean[k] - par->covariance_deltaz[k]
                        ),
                        coffe_interp_spline(
                            &bg->comoving_distance,
                            par->z_mean[k] + par->covariance_deltaz[k]
                        )
                    )
                   /pow(COFFE_H0, 3);
            }
            else{
/*
                z_mean = (par->zmin[k] + par->zmax[k]) / 2.;
                struct covariance_volume_params test;
                test.conformal_Hz = &bg->conformal_Hz;
                test.comoving_distance = &bg->comoving_distance;

                volume[k] = 4 * M_PI
                   *cov_ramp->fsky[k]
                   /coffe_integrate_1d(
                        &covariance_volume_integrand,
                        &test,
                        cov_ramp->zmin[k],
                        cov_ramp->zmax[k]
                    )
                   /pow(COFFE_H0, 3);
*/
            }

            /* TODO implement covariance for two populations */
            double galaxy_bias1 = 0, growth_rate = 0;

            /* set bias to nonzero only if there is a density contribution */
            if (par->correlation_contrib.den)
                galaxy_bias1 = coffe_interp_spline(&par->galaxy_bias1, z_mean);
            /* set growth rate to nonzero only if there's an rsd contribution */
            if (par->correlation_contrib.rsd)
                growth_rate = coffe_interp_spline(&bg->f, z_mean);

            /* b^2 + 2/3 b growth_rate + growth_rate^2/5 */
            const double c0 = pow(galaxy_bias1, 2) + 2 * galaxy_bias1 * growth_rate / 3. + pow(growth_rate, 2) / 5.;
            /* 4/3 b growth_rate + 4/7 growth_rate^2 */
            const double c2 = 4 * galaxy_bias1 * growth_rate / 3. + 4 * pow(growth_rate, 2) / 7.;
            /* 8/35 growth_rate^2 */
            const double c4 = 8 * pow(growth_rate, 2) / 35.;

            const double c0bar =
                c0*c0 + c2*c2/5. + c4*c4/9.;
            const double c2bar =
                2*c2*(7*c0 + c2)/7. + 4*c2*c4/7. + 100*c4*c4/693.;
            const double c4bar =
                18*c2*c2/35. + 2*c0*c4 + 40*c2*c4/77. + 162*c4*c4/1001.;
            const double c6bar =
                10*c4*(9*c2 + 2*c4)/99.;
            const double c8bar =
                490*c4*c4/1287.;

            const double D1z = coffe_interp_spline(&bg->D1, z_mean);

            const double coeff_array[] = {c0, c2, c4};
            const double coeffbar_array[] = {c0bar, c2bar, c4bar, c6bar, c8bar};

            for (size_t i = 0; i<par->multipole_values_len; ++i){
            for (size_t j = 0; j<par->multipole_values_len; ++j){

                /* the sums c_i wigner3j^2(l1, l2, i) */
                double coeff_sum = 0;
                for (size_t cnt = 0; cnt < COFFE_ARRAY_SIZE(coeff_array); ++cnt){
                    coeff_sum +=
                        coeff_array[cnt]
                       *pow(gsl_sf_coupling_3j(2*par->multipole_values[i], 2*par->multipole_values[j], 4*cnt, 0, 0, 0), 2);
                }

                double coeffbar_sum = 0;
                for (size_t cnt = 0; cnt < COFFE_ARRAY_SIZE(coeffbar_array); ++cnt){
                    coeffbar_sum +=
                        coeffbar_array[cnt]
                       *pow(gsl_sf_coupling_3j(2*par->multipole_values[i], 2*par->multipole_values[j], 4*cnt, 0, 0, 0), 2);
                }

                /* Kronecker delta (l1, l2) */
                double deltal1l2 =
                    (par->multipole_values[i] == par->multipole_values[j]) ? 1 : 0;

                for (size_t m = 0; m < par->sep_len; ++m){
                for (size_t n = 0; n < par->sep_len; ++n){
                    double deltaij = (m == n) ? 1 : 0;
                    /* flat-sky covariance */
                    const double result_mp_or_ramp =
                        covariance_complex(
                            par->multipole_values[i],
                            par->multipole_values[j]
                        )
                       *(
                           (2 * par->multipole_values[i] + 1)
                           *deltal1l2
                           *deltaij
                           / 2. / M_PI
                            /* trigraph */
                           *(
                                par->covariance_window
                                ?
                                1.
                                / covariance_volume_no_4pi(
                                    separations[m] - par->covariance_pixelsize[k] / 2.,
                                    separations[m] + par->covariance_pixelsize[k] / 2.
                                )
                                :
                                1.
                                / separations[n]
                                / separations[m]
                                / par->covariance_pixelsize[k]
                            )
                            /par->covariance_density[k]
                            /par->covariance_density[k]
                            +
                            /* trigraph */
                            (par->pk_type ? 1 : D1z * D1z)
                            /* trigraph */
                           *(
                                par->pk_type
                                ?
                                integral_pk[k][i*par->multipole_values_len + j][par->sep_len*n + m]
                                :
                                integral_pk[0][i*par->multipole_values_len + j][par->sep_len*n + m]
                            )
                           *coeff_sum/par->covariance_density[k]
                            +
                            /* trigraph */
                            (par->pk_type ? 1 : D1z * D1z * D1z * D1z)
                            /* TODO fix this abomination */
                            /* trigraph */
                           *(
                                par->pk_type
                                ?
                                integral_pk2[k][i*par->multipole_values_len + j][par->sep_len*n + m]
                                :
                                integral_pk2[0][i*par->multipole_values_len + j][par->sep_len*n + m]
                            )
                           *coeffbar_sum
                        )
                       /volume[k];
                    if (par->output_type == 4){
                        cov_mp->array[index].value = result_mp_or_ramp;
                    }
                    else{
                        cov_ramp->array[index].value = result_mp_or_ramp;
                    }
                    ++index;
                }
                }
            }
            }
        }

        for (size_t i = 0; i < par->covariance_coords.size; ++i){
            cov_mp->array[i].coords.z_mean = par->covariance_coords.array[i].z_mean;
            cov_mp->array[i].coords.separation1 = par->covariance_coords.array[i].separation1;
            cov_mp->array[i].coords.separation2 = par->covariance_coords.array[i].separation2;
            cov_mp->array[i].coords.l1 = par->covariance_coords.array[i].l1;
            cov_mp->array[i].coords.l2 = par->covariance_coords.array[i].l2;
        }

        /* memory cleanup */
        for (size_t index_redshift = 0; index_redshift < redshifts_to_allocate; ++index_redshift){
            for (size_t i = 0; i<par->multipole_values_len*par->multipole_values_len; ++i){
                free(integral_pk[index_redshift][i]);
                free(integral_pk2[index_redshift][i]);
            }
            free(integral_pk[index_redshift]);
            free(integral_pk2[index_redshift]);
        }
        free(separations);
        free(integral_pk);
        free(integral_pk2);
        free(volume);

        end = clock();

        if (par->verbose)
            printf("Covariance calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);

        gsl_set_error_handler(default_handler);
    }

    return EXIT_SUCCESS;
}

int coffe_covariance_free(
    coffe_covariance_array_t *cov
)
{
    if (cov->size)
        free(cov->array);
    cov->array = NULL;
    cov->size = 0;
    return EXIT_SUCCESS;
}
