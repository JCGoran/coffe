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

#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_errno.h>
#include "common.h"
#include "background.h"
#include "covariance.h"

/**
    contains the necessary integration parameters
**/
struct covariance_params
{
    struct coffe_interpolation *power_spectrum;
    double chi1, chi2, resolution;
    int l1, l2;
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
       /interp_spline(s->conformal_Hz, z)
       /pow(interp_spline(s->comoving_distance, z), 2)
       /(1 + z);
    return result;
}


/**
    finds the maximum of a double array
**/
static size_t covariance_find_maximum(size_t array[], size_t len)
{
    size_t max = array[0];

    for (size_t i = 1; i<len; ++i){
        if (array[i] > max){
           max = array[i];
        }
    }
    return max;
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
    the window function for the resolution cutoff
**/
static double covariance_window(
    double x
)
{
    return 3.0 * gsl_sf_bessel_j1(x) / x;
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
    return
        interp_spline(test->power_spectrum, k)*k*k
       *gsl_sf_bessel_jl(test->l1, k*test->chi1)
       *gsl_sf_bessel_jl(test->l2, k*test->chi2)
       *pow(covariance_window(k * test->resolution), 4);
}


/**
    integrates the above
**/
static double covariance_integral(
    struct coffe_interpolation *power_spectrum,
    double chi1, double chi2,
    int l1, int l2,
    double kmin, double kmax,
    double resolution
)
{
    struct covariance_params test;
    test.power_spectrum = power_spectrum;
    test.chi1 = chi1;
    test.chi2 = chi2;
    test.l1 = l1;
    test.l2 = l2;
    test.resolution = resolution;

    double result, error, prec = 1E-9;
    gsl_function integrand;
    integrand.function = &covariance_integrand;
    integrand.params = &test;

    gsl_integration_workspace *wspace =
        gsl_integration_workspace_alloc(COFFE_MAX_INTSPACE);
    gsl_integration_qag(
        &integrand, kmin, kmax, 0,
        prec, COFFE_MAX_INTSPACE,
        GSL_INTEG_GAUSS61, wspace,
        &result, &error
    );

    gsl_integration_workspace_free(wspace);
    return 2*result/M_PI;
}


/**
    computes the covariance of either multipoles or redshift averaged
    multipoles
**/
int coffe_covariance_init(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_covariance_t *cov_mp,
    struct coffe_covariance_t *cov_ramp
)
{
    cov_mp->flag = 0;
    cov_ramp->flag = 0;
    if (par->output_type == 4 || par->output_type == 5){
        if (par->output_type == 4){
            cov_mp->flag = 1;
        }
        else{
            cov_ramp->flag = 1;
        }
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
            cov_mp->z_mean = (double *)par->covariance_z_mean;
            cov_mp->deltaz = (double *)par->covariance_deltaz;

            cov_mp->density = (double *)par->covariance_density;
            cov_mp->fsky = (double *)par->covariance_fsky;
            cov_mp->pixelsize = (double)par->covariance_pixelsize;
            cov_mp->list_len = (size_t)par->covariance_density_len;

            cov_mp->sep = (double **)coffe_malloc(sizeof(double *)*cov_mp->list_len);
            cov_mp->sep_len = (size_t *)coffe_malloc(sizeof(size_t)*cov_mp->list_len);

            cov_mp->l_len = (size_t)par->multipole_values_len;
            cov_mp->l = (int *)coffe_malloc(sizeof(int)*cov_mp->l_len);
            for (size_t i = 0; i<cov_mp->l_len; ++i){
                cov_mp->l[i] = (int)par->multipole_values[i];
            }
        }
        else{
            cov_ramp->zmin = (double *)par->covariance_zmin;
            cov_ramp->zmax = (double *)par->covariance_zmax;

            cov_ramp->density = (double *)par->covariance_density;
            cov_ramp->fsky = (double *)par->covariance_fsky;
            cov_ramp->pixelsize = (double)par->covariance_pixelsize;
            cov_ramp->list_len = (size_t)par->covariance_density_len;

            cov_ramp->sep = (double **)coffe_malloc(sizeof(double *)*cov_ramp->list_len);
            cov_ramp->sep_len = (size_t *)coffe_malloc(sizeof(size_t)*cov_ramp->list_len);

            cov_ramp->l_len = (size_t)par->multipole_values_len;
            cov_ramp->l = (int *)coffe_malloc(sizeof(int)*cov_ramp->l_len);
            for (size_t i = 0; i<cov_ramp->l_len; ++i){
                cov_ramp->l[i] = (int)par->multipole_values[i];
            }
        }

        double *temp_spectrum_pk =
            (double *)coffe_malloc(sizeof(double)*par->power_spectrum.spline->size);
        double *temp_spectrum_pk2 =
            (double *)coffe_malloc(sizeof(double)*par->power_spectrum.spline->size);
        struct coffe_interpolation integrand_pk, integrand_pk2;

        /* setting the power spectra P(k) and P^2(k) */
        for (size_t i = 0; i<par->power_spectrum.spline->size; ++i){
            temp_spectrum_pk[i] =
                par->power_spectrum.spline->y[i];
            temp_spectrum_pk2[i] =
                par->power_spectrum.spline->y[i]*par->power_spectrum.spline->y[i];
        }

        /* interpolations of P(k) and P^2(k) */
        init_spline(
            &integrand_pk,
            par->power_spectrum.spline->x,
            temp_spectrum_pk,
            par->power_spectrum.spline->size,
            5
        );

        init_spline(
            &integrand_pk2,
            par->power_spectrum.spline->x,
            temp_spectrum_pk2,
            par->power_spectrum.spline->size,
            5
        );

        /* memory cleanup */
        free(temp_spectrum_pk);
        free(temp_spectrum_pk2);

        /* finding the largest separation */
        double *upper_limit =
            (double *)coffe_malloc(sizeof(double)*par->covariance_density_len);
        for (size_t i = 0; i<par->covariance_density_len; ++i){
            if (par->output_type == 4){
                upper_limit[i] = 2*(
                    interp_spline(
                        &bg->comoving_distance,
                        cov_mp->z_mean[i] + cov_mp->deltaz[i]
                    )
                   -interp_spline(
                        &bg->comoving_distance,
                        cov_mp->z_mean[i]
                    )
                )/COFFE_H0;
            }
            else{
                upper_limit[i] = (
                    interp_spline(
                        &bg->comoving_distance,
                        cov_ramp->zmax[i]
                    )
                   -interp_spline(
                        &bg->comoving_distance,
                        cov_ramp->zmin[i]
                    )
                )/COFFE_H0;
            }
        }

        /* max number of pixels possible for each redshift z_mean and range deltaz */
        size_t *npixels =
            (size_t *)coffe_malloc(sizeof(size_t)*par->covariance_density_len);

        for (size_t i = 0; i<par->covariance_density_len; ++i){
            npixels[i] = (size_t)(
                (upper_limit[i] - par->covariance_minimum_separation) / par->covariance_pixelsize
            );

            if (par->output_type == 4){
                cov_mp->sep[i] =
                    (double *)coffe_malloc(sizeof(double)*npixels[i]);
                cov_mp->sep_len[i] = npixels[i];
            }
            else{
                cov_ramp->sep[i] =
                    (double *)coffe_malloc(sizeof(double)*npixels[i]);
                cov_ramp->sep_len[i] = npixels[i];
            }
        }
        size_t npixels_max = covariance_find_maximum(npixels, par->covariance_density_len);

        /* all of the separations */
        double *separations = (double *)coffe_malloc(sizeof(double) * npixels_max);
        for (size_t i = 0; i<npixels_max; ++i)
            separations[i] = par->covariance_minimum_separation + i * par->covariance_pixelsize;

        /* allocating memory for the integrals of P(k) and P^2(k) (D_l1l2 and G_l1l2) */
        double **integral_pk =
            (double **)coffe_malloc(sizeof(double *)*par->multipole_values_len*par->multipole_values_len);
        double **integral_pk2 =
            (double **)coffe_malloc(sizeof(double *)*par->multipole_values_len*par->multipole_values_len);
        for (size_t i = 0; i<par->multipole_values_len; ++i){
            for (size_t j = 0; j<par->multipole_values_len; ++j){
            integral_pk[i*par->multipole_values_len + j] =
                (double *)coffe_malloc(sizeof(double)*npixels_max*npixels_max);
            integral_pk2[i*par->multipole_values_len + j] =
                (double *)coffe_malloc(sizeof(double)*npixels_max*npixels_max);
            }
        }

        /* calculating the integrals G_l1l2 and D_l1l2 (without the scale factor D1) */
        for (size_t i = 0; i<par->multipole_values_len; ++i){
            for (size_t j = i; j<par->multipole_values_len; ++j){
                #pragma omp parallel for num_threads(par->nthreads) collapse(2)
                for (size_t m = 0; m<npixels_max; ++m){
                    for (size_t n = 0; n<npixels_max; ++n){
                        integral_pk[i*par->multipole_values_len + j][npixels_max*n + m] =
                            (2*par->multipole_values[i] + 1)*(2*par->multipole_values[j] + 1)
                           *covariance_integral(
                                &integrand_pk,
                                separations[m], separations[n],
                                par->multipole_values[i], par->multipole_values[j],
                                par->k_min, par->k_max, par->covariance_pixelsize
                            )/M_PI;
                    }
                }
                #pragma omp parallel for num_threads(par->nthreads) collapse(2)
                for (size_t m = 0; m<npixels_max; ++m){
                    for (size_t n = 0; n<npixels_max; ++n){
                        integral_pk2[i*par->multipole_values_len + j][npixels_max*n + m] =
                            (2*par->multipole_values[i] + 1)*(2*par->multipole_values[j] + 1)
                           *covariance_integral(
                                &integrand_pk2,
                                separations[m], separations[n],
                                par->multipole_values[i], par->multipole_values[j],
                                par->k_min, par->k_max, par->covariance_pixelsize
                            )/2./M_PI;
                    }
                }
            }
        }

        /* allocating memory for the final result */
        if (par->output_type == 4){
            cov_mp->result = (double ***)coffe_malloc(sizeof(double **)*cov_mp->list_len);
            for (size_t k = 0; k<cov_mp->list_len; ++k){
                cov_mp->result[k] =
                    (double **)coffe_malloc(sizeof(double *)*cov_mp->l_len*cov_mp->l_len);
                for (size_t i = 0; i<cov_mp->l_len; ++i){
                    for (size_t j = 0; j<cov_mp->l_len; ++j){
                        cov_mp->result[k][i*cov_mp->l_len + j] =
                            (double *)coffe_malloc(sizeof(double)*npixels[k]*npixels[k]);
                    }
                }
            }
        }
        else{
            cov_ramp->result = (double ***)coffe_malloc(sizeof(double **)*cov_ramp->list_len);
            for (size_t k = 0; k<cov_ramp->list_len; ++k){
                cov_ramp->result[k] =
                    (double **)coffe_malloc(sizeof(double *)*cov_ramp->l_len*cov_ramp->l_len);
                for (size_t i = 0; i<cov_ramp->l_len; ++i){
                    for (size_t j = 0; j<cov_ramp->l_len; ++j){
                        cov_ramp->result[k][i*cov_ramp->l_len + j] =
                            (double *)coffe_malloc(sizeof(double)*npixels[k]*npixels[k]);
                    }
                }
            }
        }

        /* setting the transpose */
        for (size_t i = 0; i<par->multipole_values_len; ++i){
            for (size_t j = 0; j<i; ++j){
                for (size_t m = 0; m<npixels_max; ++m){
                    for (size_t n = 0; n<npixels_max; ++n){
                        integral_pk[i*par->multipole_values_len + j][npixels_max*n + m] =
                            integral_pk[j*par->multipole_values_len + i][npixels_max*m + n];
                        integral_pk2[i*par->multipole_values_len + j][npixels_max*n + m] =
                            integral_pk2[j*par->multipole_values_len + i][npixels_max*m + n];
                    }
                }
            }
        }

        double *volume =
            (double *)coffe_malloc(sizeof(double)*par->covariance_density_len);
        for (size_t k = 0; k<par->covariance_density_len; ++k){
            double z_mean;

            if (par->output_type == 4){
                z_mean = cov_mp->z_mean[k];
                volume[k] = 4*M_PI*cov_mp->fsky[k]
                   *(
                        pow(interp_spline(
                            &bg->comoving_distance, cov_mp->z_mean[k] + cov_mp->deltaz[k]), 3)
                       -pow(interp_spline(
                            &bg->comoving_distance, cov_mp->z_mean[k] - cov_mp->deltaz[k]), 3)
                    )/3./pow(COFFE_H0, 3);
            }
            else{
                z_mean = (cov_ramp->zmin[k] + cov_ramp->zmax[k])/2.;
                struct covariance_volume_params test;
                test.conformal_Hz = &bg->conformal_Hz;
                test.comoving_distance = &bg->comoving_distance;

                gsl_integration_workspace *space =
                    gsl_integration_workspace_alloc(COFFE_MAX_INTSPACE);
                gsl_function integrand;
                integrand.function = &covariance_volume_integrand;
                integrand.params = &test;
                double prec = 1E-6;
                double integral_result, integral_error;
                gsl_integration_qag(
                    &integrand, cov_ramp->zmin[k], cov_ramp->zmax[k], 0, prec,
                    COFFE_MAX_INTSPACE,
                    GSL_INTEG_GAUSS61, space,
                    &integral_result, &integral_error
                );
                gsl_integration_workspace_free(space);
                volume[k] = 4*M_PI*cov_ramp->fsky[k]/integral_result/pow(COFFE_H0, 3);
            }

            /* TODO implement covariance for two populations */
            /* TODO maybe pick a different variable name for the growth rate? */
            double matter_bias1 = 0, f = 0;

            /* set bias to nonzero only if there is a density contribution */
            if (par->correlation_contrib.den)
                matter_bias1 = interp_spline(&par->matter_bias1, z_mean);
            /* set growth rate to nonzero only if there's an rsd contribution */
            if (par->correlation_contrib.rsd)
                f = interp_spline(&bg->f, z_mean);

            /* b^2 + 2/3 b f + f^2/5 */
            double c0 = pow(matter_bias1, 2) + 2 * matter_bias1 * f / 3. + pow(f, 2) / 5.;
            /* 4/3 b f + 4/7 f^2 */
            double c2 = 4 * matter_bias1 * f / 3. + 4 * pow(f, 2) / 7.;
            /* 8/35 f^2 */
            double c4 = 8 * pow(f, 2) / 35.;

            double c0bar =
                c0*c0 + c2*c2/5. + c4*c4/9.;
            double c2bar =
                2*c2*(7*c0 + c2)/7. + 4*c2*c4/7. + 100*c4*c4/693.;
            double c4bar =
                18*c2*c2/35. + 2*c0*c4 + 40*c2*c4/77. + 162*c4*c4/1001.;
            double c6bar =
                10*c4*(9*c2 + 2*c4)/99.;
            double c8bar =
                490*c4*c4/1287.;

            double D10 = interp_spline(&bg->D1, 0);
            double D1z = interp_spline(&bg->D1, z_mean);

            double coeff_array[] = {c0, c2, c4};
            double coeffbar_array[] = {c0bar, c2bar, c4bar, c6bar, c8bar};
            double coeff_sum = 0;
            double coeffbar_sum = 0;

            for (size_t i = 0; i<par->multipole_values_len; ++i){
                for (size_t j = 0; j<par->multipole_values_len; ++j){
                    /* the sums c_i wigner3j^2(l1, l2, i) */
                    for (size_t cnt = 0; cnt<sizeof(coeff_array)/sizeof(coeff_array[0]); ++cnt){
                        coeff_sum +=
                            coeff_array[cnt]
                           *pow(gsl_sf_coupling_3j(2*par->multipole_values[i], 2*par->multipole_values[j], 4*cnt, 0, 0, 0), 2);
                    }
                    for (size_t cnt = 0; cnt<sizeof(coeffbar_array)/sizeof(coeffbar_array[0]); ++cnt){
                        coeffbar_sum +=
                            coeffbar_array[cnt]
                           *pow(gsl_sf_coupling_3j(2*par->multipole_values[i], 2*par->multipole_values[j], 4*cnt, 0, 0, 0), 2);
                    }

                    double deltal1l2;
                    if (par->multipole_values[i] == par->multipole_values[j]) deltal1l2 = 1;
                    else deltal1l2 = 0;

                    for (size_t m = 0; m<npixels[k]; ++m){
                        for (size_t n = 0; n<npixels[k]; ++n){
                            double deltaij;
                            if (m == n) deltaij = 1;
                            else deltaij = 0;
                            /* flat-sky covariance */
                            const double result_mp_or_ramp =
                                covariance_complex(
                                    par->multipole_values[i],
                                    par->multipole_values[j]
                            )
                               *(
                                    (2*par->multipole_values[i] + 1)*deltaij*deltal1l2
                                   / 2. / M_PI
                                   /par->covariance_density[k] / par->covariance_density[k]
                                   /par->covariance_pixelsize
                                   /separations[n]
                                   /separations[m]
                                   +D1z * D1z / D10 / D10
                                   *integral_pk[i*par->multipole_values_len + j][npixels_max*n + m]
                                   *coeff_sum/par->covariance_density[k]
                                   +D1z * D1z * D1z * D1z / D10 / D10 / D10 / D10
                                   *integral_pk2[i*par->multipole_values_len + j][npixels_max*n + m]
                                   *coeffbar_sum
                                )
                               /volume[k];
                            if (par->output_type == 4){
                                cov_mp->sep[k][m] = separations[m];
                                cov_mp->sep[k][n] = separations[n];
                                cov_mp->result[k][par->multipole_values_len*i + j][npixels[k]*n + m] = result_mp_or_ramp;
                            }
                            else{
                                cov_ramp->sep[k][m] = separations[m];
                                cov_ramp->sep[k][n] = separations[n];
                                cov_ramp->result[k][par->multipole_values_len*i + j][npixels[k]*n + m] = result_mp_or_ramp;
                            }
                        }
                    }
                    coeff_sum = 0, coeffbar_sum = 0;
                }
            }
        }

        /* memory cleanup */
        for (size_t i = 0; i<par->multipole_values_len*par->multipole_values_len; ++i){
            free(integral_pk[i]);
            free(integral_pk2[i]);
        }
        free(separations);
        free(integral_pk);
        free(integral_pk2);
        gsl_spline_free(integrand_pk.spline);
        gsl_interp_accel_free(integrand_pk.accel);
        gsl_spline_free(integrand_pk2.spline);
        gsl_interp_accel_free(integrand_pk2.accel);

        end = clock();

        if (par->verbose)
            printf("Covariance calculated in %.2f s\n",
                (double)(end - start) / CLOCKS_PER_SEC);

        gsl_set_error_handler(default_handler);
    }

    return EXIT_SUCCESS;
}

int coffe_covariance_free(
    struct coffe_covariance_t *cov
)
{
    if (cov->flag){
        for (size_t k = 0; k<cov->list_len; ++k){
            for (size_t i = 0; i<cov->l_len; ++i){
                for (size_t j = 0; j<cov->l_len; ++j){
                    free(cov->result[k][i*cov->l_len + j]);
                }
            }
            free(cov->result[k]);
        }
        free(cov->result);

        for (size_t i = 0; i<cov->list_len; ++i){
            free(cov->sep[i]);
        }
        free(cov->sep);
        free(cov->sep_len);
        free(cov->l);
        cov->flag = 0;
    }
    return EXIT_SUCCESS;
}
