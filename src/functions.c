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
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include "common.h"
#include "background.h"
#include "integrals.h"
#include "functions.h"

#ifdef HAVE_CLASS
#include "class.h"
#endif


/*
    explanation for the future reader:
    There are several cases to consider for each of the `functions_*`, depending on whether or
    not CLASS is included (hence the barrage of `#ifdef`s everywhere), and whether
    or not the user requested the flatsky approximation, and even whether or not the output_type
    is multipoles or not:
    1. flatsky
        1.a. with CLASS
            1.a.I. power spectrum generated at z = 0
            1.a.II. let CLASS generate the power spectrum at the _MEAN_
                    redshift (so NOT P(k, z₁, z₂), but P(k, (z₁ + z₂) / 2).
                    This case includes 3 ways to generate it: linear (pk_type = 1),
                    nonlinear halofit (pk_type = 2), and nonlinear HMcode (pk_type = 2).
                    Note that P(k) is not computed here (to save
                    computational time), but rather in `integrals.c`, where we make
                    a bicubic spline in (k, z), and here we just interpolate it.
                    In practice, we don't interpolate P(k, z), but rather
                    its Fourier-Bessel transform, I^n_l(r, z).
                    You may ask why not simply interpolate P(k, z₁, z₂) instead of P(k, z).
                    There are a couple of problems with this:
                    - CLASS can only generate the nonlinear POWER spectrum P(k, z),
                      but not the CROSS spectrum P(k, z₁, z₂). In principle, we may assume
                      that P(k, z₁, z₂) ≈ √(P(k, z₁)×P(k, z₂), (or even better, use the Zeldovich
                      approximation, see arXiv:1905.02078, eq. (13)), but this
                      brings us to the below problem:
                    - as of the time of writing, GSL, the C library used for 1D
                      and 2D interpolation, does not support 3D interpolation.
                      Even 2D interpolation is very computationally expensive
                      and suffers from worse numerical errors compared to its 1D counterpart,
                      and I imagine the 3D version would be even worse
                    - finally, the (midpoint) approximation P(k, z₁, z₂) ≈ P(k, (z₁ + z₂) / 2),
                      turns out to be almost unreasonably effective (even for lensing),
                      see arXiv:2011.06185, notably their figures 9 and 10
        1.b. without CLASS. In this case, since COFFE can't generate
             P(k) on its own, it assumes that the user gave it a file with
             P(k, z = 0), and setting any nonlinear parameters doesn't work.
    2. full-sky
        2.a. with CLASS
            2.a.I. same as 1.a.I.
            2.a.II. same as 1.a.II.
        2.b. without CLASS. Same as 1.b., i.e. nonlinear parameters do NOT have any effect.

    Note that in flat-sky, the multipoles of lensing-lensing and density-lensing can be
    computed analytically (see arXiv:2011.01878), so they are in their own separate functions.
    In principle, all of the nonintegrated terms have straightforward flat-sky
    expressions for multipoles, but the 1D integrator is so fast that it currently
    isn't practical to make separate functions (though it would be nice from a
    software design perspective, but alas, there are only so many hours in a day).
*/


/**
    all the nonintegrated terms in one place
**/

double functions_nonintegrated(
    const coffe_parameters_t *par,
    const coffe_background_t *bg,
    const coffe_integral_array_t *integral,
    const double z_mean,
    const double mu,
    const double sep
)
{
    const double chi_mean = coffe_interp_spline(&bg->comoving_distance, z_mean);
    const double chi1 = chi_mean - sep*mu/2.;
    const double chi2 = chi_mean + sep*mu/2.;
    const double costheta =
        (2.*chi_mean*chi_mean - sep*sep + mu*mu*sep*sep/2.)
       /(2.*chi_mean*chi_mean - mu*mu*sep*sep/2.);

    double result = 0;
    const double z1 = coffe_interp_spline(&bg->z_as_chi, chi1);
    const double z2 = coffe_interp_spline(&bg->z_as_chi, chi2);
    const double f1 = coffe_interp_spline(&bg->f, z1);
    const double f2 = coffe_interp_spline(&bg->f, z2);
    const double fmean = coffe_interp_spline(&bg->f, z_mean);
    const double curlyH1 = coffe_interp_spline(&bg->conformal_Hz, z1); // dimensionless
    const double curlyH2 = coffe_interp_spline(&bg->conformal_Hz, z2); // dimensionless
    const double curlyH_mean = coffe_interp_spline(&bg->conformal_Hz, z_mean); // dimensionless
    const double b1 = coffe_interp_spline(&par->galaxy_bias1, z1);
    const double b2 = coffe_interp_spline(&par->galaxy_bias2, z2);
    const double bz_mean1 = coffe_interp_spline(&par->galaxy_bias1, z_mean);
    const double bz_mean2 = coffe_interp_spline(&par->galaxy_bias2, z_mean);
    const double G1 = coffe_interp_spline(&bg->G1, z1);
    const double G2 = coffe_interp_spline(&bg->G2, z2);
    const double Gz_mean1 = coffe_interp_spline(&bg->G1, z_mean);
    const double Gz_mean2 = coffe_interp_spline(&bg->G2, z_mean);
    const double s1 = coffe_interp_spline(&par->magnification_bias1, z1);
    const double s2 = coffe_interp_spline(&par->magnification_bias2, z2);
    const double fevo1 = coffe_interp_spline(&par->evolution_bias1, z1);
    const double fevo2 = coffe_interp_spline(&par->evolution_bias2, z2);
    const double a1 = coffe_interp_spline(&bg->a, z1);
    const double a2 = coffe_interp_spline(&bg->a, z2);

    /* den-den term */
    if (
        par->correlation_contrib.den &&
        !par->only_cross_correlations
    ){
        /* den-den modified by flatsky */
        if (par->flatsky_local){
#ifdef HAVE_CLASS
            if (
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){
                result +=
                    bz_mean1
                   *bz_mean2
                   *coffe_interp_spline2d(
                        &coffe_find_integral(
                            integral,
                            0,
                            0,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->result2d,
                        z_mean,
                        sep
                    );
            }
            else{
#endif
            if (par->b_derivative || par->b_tilde_derivative){
            result += 2 * bz_mean1
               *coffe_interp_spline(
                    &coffe_find_integral(
                        integral,
                        0,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->result,
                    sep
                );
            }
            else if (par->f_derivative || par->f_tilde_derivative)
                result += 0;
            else{
            result +=
                bz_mean1
               *bz_mean2
               *coffe_interp_spline(
                    &coffe_find_integral(
                        integral,
                        0,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->result,
                    sep
                );
            }
#ifdef HAVE_CLASS
            }
#endif
        }
        else{
#ifdef HAVE_CLASS
            if (
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){
                result +=
                    b1
                   *b2
                   *coffe_interp_spline2d(
                        &coffe_find_integral(
                            integral,
                            0,
                            0,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->result2d,
                        z_mean,
                        sep
                    );
            }
            else
#endif
                result +=
                    b1
                   *b2
                   *coffe_interp_spline(
                        &coffe_find_integral(
                        integral,
                            0,
                            0,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->result,
                        sep
                    );
        }
    }
    /* rsd-rsd term */
    if (
        par->correlation_contrib.rsd &&
        !par->only_cross_correlations
    ){
        /* rsd-rsd modified by flatsky */
        if (par->flatsky_local){
#ifdef HAVE_CLASS
            if (
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){
                result +=
                    fmean * fmean
                   *coffe_interp_spline2d(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result2d, z_mean, sep) / 5.
                    -
                    4 * fmean * fmean
                   *coffe_interp_spline2d(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result2d, z_mean, sep) / 7.
                   *gsl_sf_legendre_P2(mu)
                    +
                    8 * fmean * fmean
                   *coffe_interp_spline2d(&coffe_find_integral(integral, 0, 4, COFFE_INTEGER, COFFE_INTEGER)->result2d, z_mean, sep) / 35.
                   *gsl_sf_legendre_Pl(4, mu);
            }
            else{
#endif
            if (par->b_derivative || par->b_tilde_derivative)
                result += 0;
            else if (par->f_derivative || par->f_tilde_derivative){
            result +=
                2 * fmean
               *coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep) / 5.
                -
                8 * fmean
               *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sep) / 7.
               *gsl_sf_legendre_P2(mu)
                +
                16 * fmean
               *coffe_interp_spline(&coffe_find_integral(integral, 0, 4, COFFE_INTEGER, COFFE_INTEGER)->result, sep) / 35.
               *gsl_sf_legendre_Pl(4, mu);
            }
            else{
            result +=
                fmean * fmean
               *coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep) / 5.
                -
                4 * fmean * fmean
               *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sep) / 7.
               *gsl_sf_legendre_P2(mu)
                +
                8 * fmean * fmean
               *coffe_interp_spline(&coffe_find_integral(integral, 0, 4, COFFE_INTEGER, COFFE_INTEGER)->result, sep) / 35.
               *gsl_sf_legendre_Pl(4, mu);
            }
#ifdef HAVE_CLASS
            }
#endif
        }
        else{
#ifdef HAVE_CLASS
            if (
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){
                result +=
                    f1*f2*(1 + 2*pow(costheta, 2))/15
                   *coffe_interp_spline2d(
                        &coffe_find_integral(
                            integral,
                            0,
                            0,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->result2d,
                        z_mean,
                        sep
                    )
                    -
                    f1*f2/21.*(
                        (1 + 11.*pow(costheta, 2)) + 18*costheta*(pow(costheta, 2) - 1)*chi1*chi2/sep/sep
                    )
                   *coffe_interp_spline2d(
                        &coffe_find_integral(
                            integral,
                            0,
                            2,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->result2d,
                        z_mean,
                        sep
                    )
                    +
                    f1*f2*(
                        4*(3*pow(costheta, 2) - 1)*(pow(chi1, 4) + pow(chi2, 4))/35./pow(sep, 4)
                        +
                        chi1*chi2*(3 + pow(costheta, 2))*(
                            3*(3 + pow(costheta, 2))*chi1*chi2 - 8*(pow(chi1, 2) + pow(chi2, 2))*costheta
                        )/35./pow(sep, 4)
                    )
                   *coffe_interp_spline2d(
                        &coffe_find_integral(
                            integral,
                            0,
                            4,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->result2d,
                        z_mean,
                        sep
                    );
            }
            else
#endif
                result +=
                    f1*f2*(1 + 2*pow(costheta, 2))/15
                   *coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                    -
                    f1*f2/21.*(
                        (1 + 11.*pow(costheta, 2)) + 18*costheta*(pow(costheta, 2) - 1)*chi1*chi2/sep/sep
                    )
                   *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                    +
                    f1*f2*(
                        4*(3*pow(costheta, 2) - 1)*(pow(chi1, 4) + pow(chi2, 4))/35./pow(sep, 4)
                        +
                        chi1*chi2*(3 + pow(costheta, 2))*(
                            3*(3 + pow(costheta, 2))*chi1*chi2 - 8*(pow(chi1, 2) + pow(chi2, 2))*costheta
                        )/35./pow(sep, 4)
                    )
                   *coffe_interp_spline(&coffe_find_integral(integral, 0, 4, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
        }
    }
    /* d1-d1 term */
    if (
        par->correlation_contrib.d1 &&
        !par->only_cross_correlations
    ){
        /* TODO should this be yet another flag, or go under standard? */
        if (par->flatsky_local){
#ifdef HAVE_CLASS
            if (
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){
            result +=
             (
                curlyH_mean
               *curlyH_mean
               *fmean
               *fmean
               *Gz_mean1
               *Gz_mean2
               *coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        2,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->result2d,
                    z_mean,
                    sep
                )
               /3.
                +
                2
               *curlyH_mean
               *curlyH_mean
               *fmean
               *fmean
               *Gz_mean1
               *Gz_mean2
               *pow(sep, 2)
               *coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        2,
                        2,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->result2d,
                    z_mean,
                    sep
                )
               *gsl_sf_legendre_P2(mu)
               /3.
            );
            }
            else{
#endif
            result +=
             (
                curlyH_mean
               *curlyH_mean
               *fmean
               *fmean
               *Gz_mean1
               *Gz_mean2
               *coffe_interp_spline(
                    &coffe_find_integral(
                        integral,
                        2,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->result,
                    sep
                )
               /3.
                +
                2
               *curlyH_mean
               *curlyH_mean
               *fmean
               *fmean
               *Gz_mean1
               *Gz_mean2
               *pow(sep, 2)
               *coffe_interp_spline(
                    &coffe_find_integral(
                        integral,
                        2,
                        2,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->result,
                    sep
                )
               *gsl_sf_legendre_P2(mu)
               /3.
            );
#ifdef HAVE_CLASS
            }
#endif
        }
        else{
#ifdef HAVE_CLASS
            if(
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){
                result +=
                    (
                        curlyH1
                       *curlyH2
                       *f1
                       *f2
                       *G1
                       *G2
                       *costheta
                       *coffe_interp_spline2d(
                            &coffe_find_integral(
                                integral,
                                2,
                                0,
                                COFFE_INTEGER,
                                COFFE_INTEGER
                            )->result2d,
                            z_mean,
                            sep
                        )
                       /3.
                       -
                        curlyH1
                       *curlyH2
                       *f1
                       *f2
                       *G1
                       *G2
                       *(
                            (chi2 - chi1 * costheta)
                           *(chi1 - chi2 * costheta)
                            +
                            pow(sep, 2)
                           *costheta
                           /3.
                        )
                       *coffe_interp_spline2d(
                            &coffe_find_integral(
                                integral,
                                2,
                                2,
                                COFFE_INTEGER,
                                COFFE_INTEGER
                            )->result2d,
                            z_mean,
                            sep
                        )
                    );
            }
            else{
#endif
            result +=
                (
                    curlyH1
                   *curlyH2
                   *f1
                   *f2
                   *G1
                   *G2
                   *costheta
                   *coffe_interp_spline(
                        &coffe_find_integral(
                            integral,
                            2,
                            0,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->result,
                        sep
                    )
                   /3.
                   -
                    curlyH1
                   *curlyH2
                   *f1
                   *f2
                   *G1
                   *G2
                   *(
                        (chi2 - chi1 * costheta)
                       *(chi1 - chi2 * costheta)
                        +
                        pow(sep, 2)
                       *costheta
                       /3.
                    )
                   *coffe_interp_spline(
                        &coffe_find_integral(
                            integral,
                            2,
                            2,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->result,
                        sep
                    )
                );
#ifdef HAVE_CLASS
            }
#endif
        }
    }
    /* d2-d2 term */
    if (
        par->correlation_contrib.d2 &&
        !par->only_cross_correlations
    ){
        result +=
            (3 - fevo1)*(3 - fevo2)*pow(curlyH1, 2)*pow(curlyH2, 2)*f1*f2
           *(
                coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                /* renormalization term */
                -
                coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        4,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->renormalization,
                    chi1,
                    chi2
                )
            );
    }
    /* g1-g1 term */
    if (
        par->correlation_contrib.g1 &&
        !par->only_cross_correlations
    ){
        result += 9*pow((par->Omega0_cdm + par->Omega0_baryon), 2)
           *(1 + G1)*(1 + G2)/4/a1/a2
           *(
                coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                /* renormalization term */
                -
                coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        4,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->renormalization,
                    chi1,
                    chi2
                )
            );
    }
    /* g2-g2 term */
    if (
        par->correlation_contrib.g2 &&
        !par->only_cross_correlations
    ){
        result += 9*pow((par->Omega0_cdm + par->Omega0_baryon), 2)
           *(5*s1 - 2)*(5*s2 - 2)/4/a1/a2
           *(
                coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                /* renormalization term */
                -
                coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        4,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->renormalization,
                    chi1,
                    chi2
                )
            );
    }
    /* g3-g3 term */
    if (
        par->correlation_contrib.g3 &&
        !par->only_cross_correlations
    ){
        result += 9*pow((par->Omega0_cdm + par->Omega0_baryon), 2)
           *(f1 - 1)*(f2 - 1)/4/a1/a2
           *(
                coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                /* renormalization term */
                -
                coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        4,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->renormalization,
                    chi1,
                    chi2
                )
            );
    }
    /* den-rsd + rsd-den term */
    if (
        par->correlation_contrib.den &&
        par->correlation_contrib.rsd
    ){
        /* den-rsd modified by flatsky */
        if (par->flatsky_local){
#ifdef HAVE_CLASS
            if (
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){
            result +=
                (bz_mean1 * fmean / 3. + bz_mean2 * fmean / 3.)
               *coffe_interp_spline2d(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result2d, z_mean, sep)
               -
                (2 * bz_mean1 * fmean / 3. + 2 * bz_mean2 * fmean / 3.)
               *coffe_interp_spline2d(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result2d, z_mean, sep)
               *gsl_sf_legendre_P2(mu);
            }
            else{
#endif
            if (par->b_derivative || par->b_tilde_derivative){
            result +=
                (fmean/3. + fmean/3.)
               *coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
               -
                (2*fmean/3. + 2*fmean/3.)
               *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
               *gsl_sf_legendre_P2(mu);
            }
            else if (par->f_derivative || par->f_tilde_derivative){
            result +=
                (bz_mean1/3. + bz_mean2/3.)
               *coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
               -
                (2*bz_mean1/3. + 2*bz_mean2/3.)
               *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
               *gsl_sf_legendre_P2(mu);
            }
            else{
            result +=
                (bz_mean1*fmean/3. + bz_mean2*fmean/3.)
               *coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
               -
                (2*bz_mean1*fmean/3. + 2*bz_mean2*fmean/3.)
               *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
               *gsl_sf_legendre_P2(mu);
            }
#ifdef HAVE_CLASS
            }
#endif
        }
        else{
#ifdef HAVE_CLASS
            if (
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){
                result += (b1 * f2 / 3. + b2 * f1 / 3.)
                   *coffe_interp_spline2d(
                        &coffe_find_integral(
                            integral,
                            0,
                            0,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->result2d,
                        z_mean,
                        sep
                    )
                    -
                    (
                        b1 * f2 * (2. / 3. - (1. - pow(costheta, 2)) * pow(chi1 / sep, 2))
                        +
                        b2*f1*(2. / 3. - (1. - pow(costheta, 2))*pow(chi2 / sep, 2))
                    )
                   *coffe_interp_spline2d(
                        &coffe_find_integral(
                            integral,
                            0,
                            2,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->result2d,
                        z_mean,
                        sep
                    );
            }
            else
#endif
                result += (b1 * f2 / 3. + b2 * f1 / 3.)
                   *coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                   -
                    (
                        b1 * f2 * (2. / 3. - (1. - pow(costheta, 2)) * pow(chi1 / sep, 2))
                        +
                        b2 * f1 * (2. / 3. - (1. - pow(costheta, 2)) * pow(chi2 / sep, 2))
                    )
                   *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
        }
    }
    /* den-d1 + d1-den term */
    if (
        par->correlation_contrib.den &&
        par->correlation_contrib.d1
    ){
        if (par->flatsky_local){
#ifdef HAVE_CLASS
            if (
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){
            result += (
                    bz_mean1
                   *fmean
                   *curlyH_mean
                   *Gz_mean2
                   *gsl_sf_legendre_P1(mu)
                    +
                    bz_mean2
                   *fmean
                   *curlyH_mean
                   *Gz_mean1
                   *gsl_sf_legendre_P1(-mu)
                )
               *sep
               *coffe_interp_spline2d(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result2d, z_mean, sep);
            }
            else
#endif
            result += (
                    bz_mean1
                   *fmean
                   *curlyH_mean
                   *Gz_mean2
                   *gsl_sf_legendre_P1(mu)
                    +
                    bz_mean2
                   *fmean
                   *curlyH_mean
                   *Gz_mean1
                   *gsl_sf_legendre_P1(-mu)
                )
               *sep
               *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
        }
        else{
#ifdef HAVE_CLASS
            if (
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){
            result += -(
                    b1
                   *f2
                   *curlyH2
                   *G2
                   *(chi1 * costheta - chi2)
                    +
                    b2
                   *f1
                   *curlyH1
                   *G1
                   *(chi2 * costheta - chi1)
                )
               *coffe_interp_spline2d(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result2d, z_mean, sep);
            }
            else
#endif
            result += -(
                    b1
                   *f2
                   *curlyH2
                   *G2
                   *(chi1 * costheta - chi2)
                    +
                    b2
                   *f1
                   *curlyH1
                   *G1
                   *(chi2 * costheta - chi1)
                )
               *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
        }
    }
    /* den-d2 + d2-den term */
    if (
        par->correlation_contrib.den &&
        par->correlation_contrib.d2
    ){
        result += (
                (3 - fevo2)*b1*f2*pow(curlyH2, 2)
                +
                (3 - fevo1)*b2*f1*pow(curlyH1, 2)
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
    }
    /* den-g1 + g1-den term */
    if (
        par->correlation_contrib.den &&
        par->correlation_contrib.g1
    ){
        result += -(
                b1*3*(par->Omega0_cdm + par->Omega0_baryon)/2/a2*(1 + G2)
                +
                b2*3*(par->Omega0_cdm + par->Omega0_baryon)/2/a1*(1 + G1)
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
    }
    /* den-g2 + g2-den term */
    if (
        par->correlation_contrib.den &&
        par->correlation_contrib.g2
    ){
        result += -(
                b1*3*(par->Omega0_cdm + par->Omega0_baryon)/2/a2*(5*s2 - 2)
                +
                b2*3*(par->Omega0_cdm + par->Omega0_baryon)/2/a1*(5*s1 - 2)
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
    }
    /* den-g3 + g3-den term */
    if (
        par->correlation_contrib.den &&
        par->correlation_contrib.g3
    ){
        result += -(
                b1*3*(par->Omega0_cdm + par->Omega0_baryon)/2/a2*(f2 - 1)
                +
                b2*3*(par->Omega0_cdm + par->Omega0_baryon)/2/a1*(f1 - 1)
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
    }
    /* rsd-d1 + d1-rsd term */
    if (
        par->correlation_contrib.rsd &&
        par->correlation_contrib.d1
    ){
        if (par->flatsky_local){
#ifdef HAVE_CLASS
            if (
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){
                result += (
                    -3. / 5. * (
                        fmean
                       *fmean
                       *curlyH_mean
                       *Gz_mean2
                       *gsl_sf_legendre_P1(-mu)
                        +
                        fmean
                       *fmean
                       *curlyH_mean
                       *Gz_mean1
                       *gsl_sf_legendre_P1(mu)
                    )
                   *coffe_interp_spline2d(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result2d, z_mean, sep)
                   *sep
                    + (
                        fmean
                       *fmean
                       *curlyH_mean
                       *Gz_mean2
                       *2. / 5. * gsl_sf_legendre_Pl(3, -mu)
                        +
                        fmean
                       *fmean
                       *curlyH_mean
                       *Gz_mean1
                       *2. / 5. * gsl_sf_legendre_Pl(3, mu)
                    )
                   *coffe_interp_spline2d(&coffe_find_integral(integral, 1, 3, COFFE_INTEGER, COFFE_INTEGER)->result2d, z_mean, sep)
                   *sep
                );
            }
            else
#endif
            result += (
                -3. / 5. * (
                    fmean
                   *fmean
                   *curlyH_mean
                   *Gz_mean2
                   *gsl_sf_legendre_P1(-mu)
                    +
                    fmean
                   *fmean
                   *curlyH_mean
                   *Gz_mean1
                   *gsl_sf_legendre_P1(mu)
                )
               *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
               *sep
                + (
                    fmean
                   *fmean
                   *curlyH_mean
                   *Gz_mean2
                   *2. / 5. * gsl_sf_legendre_Pl(3, -mu)
                    +
                    fmean
                   *fmean
                   *curlyH_mean
                   *Gz_mean1
                   *2. / 5. * gsl_sf_legendre_Pl(3, mu)
                )
               *coffe_interp_spline(&coffe_find_integral(integral, 1, 3, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
               *sep
            );
        }
        else{
#ifdef HAVE_CLASS
            if (
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){
                result += (
                    (
                        f1
                       *f2
                       *curlyH2
                       *G2
                       *(
                            (1. + 2 * pow(costheta, 2)) * chi2
                           -3 * chi1 * costheta
                        )
                       /5.
                        +
                        f2
                       *f1
                       *curlyH1
                       *G1
                       *(
                            (1. + 2 * pow(costheta, 2)) * chi1
                           -3 * chi2 * costheta
                        )
                       /5.
                    )
                   *coffe_interp_spline2d(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result2d, z_mean, sep)
                    + (
                        f1
                       *f2
                       *curlyH2
                       *G2
                       *(
                            (1. - 3 * costheta * costheta)
                           *pow(chi2, 3)
                            +
                            costheta
                           *(5. + pow(costheta, 2))
                           *pow(chi2, 2)
                           *chi1
                            -
                            2
                           *(2. + pow(costheta, 2))
                           *chi2
                           *pow(chi1, 2)
                            +
                            2
                           *pow(chi1, 3)
                           *costheta
                        )
                       /5
                        +
                        f2
                       *f1
                       *curlyH1
                       *G1
                       *(
                            (1. - 3 * costheta * costheta)
                           *pow(chi1, 3)
                            +
                            costheta
                           *(5. + pow(costheta, 2))
                           *pow(chi1, 2)
                           *chi2
                            -
                            2
                           *(2. + pow(costheta, 2))
                           *chi1
                           *pow(chi2, 2)
                            +
                            2
                           *pow(chi2, 3)
                           *costheta
                        )
                       /5
                    )
                   *coffe_interp_spline2d(&coffe_find_integral(integral, 1, 3, COFFE_INTEGER, COFFE_INTEGER)->result2d, z_mean, sep)
                   /pow(sep, 2)
                );
            }
            else
#endif
            result += (
                (
                    f1
                   *f2
                   *curlyH2
                   *G2
                   *(
                        (1. + 2 * pow(costheta, 2)) * chi2
                       -3 * chi1 * costheta
                    )
                   /5.
                    +
                    f2
                   *f1
                   *curlyH1
                   *G1
                   *(
                        (1. + 2 * pow(costheta, 2)) * chi1
                       -3 * chi2 * costheta
                    )
                   /5.
                )
               *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                + (
                    f1
                   *f2
                   *curlyH2
                   *G2
                   *(
                        (1. - 3 * costheta * costheta)
                       *pow(chi2, 3)
                        +
                        costheta
                       *(5. + pow(costheta, 2))
                       *pow(chi2, 2)
                       *chi1
                        -
                        2
                       *(2. + pow(costheta, 2))
                       *chi2
                       *pow(chi1, 2)
                        +
                        2
                       *pow(chi1, 3)
                       *costheta
                    )
                   /5
                    +
                    f2
                   *f1
                   *curlyH1
                   *G1
                   *(
                        (1. - 3 * costheta * costheta)
                       *pow(chi1, 3)
                        +
                        costheta
                       *(5. + pow(costheta, 2))
                       *pow(chi1, 2)
                       *chi2
                        -
                        2
                       *(2. + pow(costheta, 2))
                       *chi1
                       *pow(chi2, 2)
                        +
                        2
                       *pow(chi2, 3)
                       *costheta
                    )
                   /5
                )
               *coffe_interp_spline(&coffe_find_integral(integral, 1, 3, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
               /pow(sep, 2)
            );
        }
    }
    /* rsd-d2 + d2-rsd term */
    if (
        par->correlation_contrib.rsd &&
        par->correlation_contrib.d2
    ){
        result += (
            (
                (3 - fevo2)/3*f1*f2*pow(curlyH2, 2)
                +
                (3 - fevo1)/3*f2*f1*pow(curlyH1, 2)
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
            - (
                (3 - fevo2)*f1*f2*pow(curlyH2, 2)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi2, 2))
                +
                (3 - fevo1)*f2*f1*pow(curlyH1, 2)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi1, 2))
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
        );
    }
    /* rsd-g1 + g1-rsd term */
    if (
        par->correlation_contrib.rsd &&
        par->correlation_contrib.g1
    ){
        result += -(
                (par->Omega0_cdm + par->Omega0_baryon)/2./a2*f1*(1 + G2)
                +
                (par->Omega0_cdm + par->Omega0_baryon)/2./a1*f2*(1 + G1)
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
            + (
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a2*f1*(1 + G2)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi2, 2))
                +
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a1*f2*(1 + G1)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi1, 2))
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
    }
    /* rsd-g2 + g2-rsd term */
    if (
        par->correlation_contrib.rsd &&
        par->correlation_contrib.g2
    ){
        result += -(
                (par->Omega0_cdm + par->Omega0_baryon)/2./a2*f1*(5*s2 - 2)
                +
                (par->Omega0_cdm + par->Omega0_baryon)/2./a1*f2*(5*s1 - 2)
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
            + (
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a2*f1*(5*s2 - 2)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi2, 2))
                +
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a1*f2*(5*s1 - 2)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi1, 2))
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
    }
    /* rsd-g3 + g3-rsd term */
    if (
        par->correlation_contrib.rsd &&
        par->correlation_contrib.g3
    ){
        result += -(
                (par->Omega0_cdm + par->Omega0_baryon)/2./a2*f1*(f2 - 1)
                +
                (par->Omega0_cdm + par->Omega0_baryon)/2./a1*f2*(f1 - 1)
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
            + (
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a2*f1*(f2 - 1)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi2, 2))
                +
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a1*f2*(f1 - 1)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi1, 2))
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sep);

    }
    /* d1-d2 + d2-d1 term */
    if (
        par->correlation_contrib.d1 &&
        par->correlation_contrib.d2
    ){
        result += -(
                (3 - fevo2)*curlyH1*pow(curlyH2, 2)*f1*f2*(chi2*costheta - chi1)
                +
                (3 - fevo1)*curlyH2*pow(curlyH1, 2)*f2*f1*(chi1*costheta - chi2)
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
    }
    /* d1-g1 + g1-d1 term */
    if (
        par->correlation_contrib.d1 &&
        par->correlation_contrib.g1
    ){
        result += (
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a2*curlyH1*f1*(1 + G2)*(chi2*costheta - chi1)
                +
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a1*curlyH2*f2*(1 + G1)*(chi1*costheta - chi2)
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
    }
    /* d1-g2 + g2-d1 term */
    if (
        par->correlation_contrib.d1 &&
        par->correlation_contrib.g2
    ){
        result += (
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a2*curlyH1*f1*(5*s2 - 2)*(chi2*costheta - chi1)
                +
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a1*curlyH2*f2*(5*s1 - 2)*(chi1*costheta - chi2)
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
    }
    /* d1-g3 + g3-d1 term */
    if (
        par->correlation_contrib.d1 &&
        par->correlation_contrib.g3
    ){
        result += (
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a2*curlyH1*f1*(f2 - 1.)*(chi2*costheta - chi1)
                +
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a1*curlyH2*f2*(f1 - 1.)*(chi1*costheta - chi2)
            )
           *coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sep);
    }
    /* d2-g1 + g1-d2 term */
    if (
        par->correlation_contrib.d2 &&
        par->correlation_contrib.g1
    ){
        result += -(
                3*(3 - fevo1)*(par->Omega0_cdm + par->Omega0_baryon)/2./a2*pow(curlyH1, 2)*f1*(1 + G2)
                +
                3*(3 - fevo2)*(par->Omega0_cdm + par->Omega0_baryon)/2./a1*pow(curlyH2, 2)*f2*(1 + G1)
            )
           *(
                coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                /* renormalization term */
                -
                coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        4,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->renormalization,
                    chi1,
                    chi2
                )
            );
    }
    /* d2-g2 + g2-d2 term */
    if (
        par->correlation_contrib.d2 &&
        par->correlation_contrib.g2
    ){
        result += -(
                3*(3 - fevo1)*(par->Omega0_cdm + par->Omega0_baryon)/2./a2*pow(curlyH1, 2)*f1*(5*s2 - 2)
                +
                3*(3 - fevo2)*(par->Omega0_cdm + par->Omega0_baryon)/2./a1*pow(curlyH2, 2)*f2*(5*s1 - 2)
            )
           *(
                coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                /* renormalization term */
                -
                coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        4,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->renormalization,
                    chi1,
                    chi2
                )
            );
    }
    /* d2-g3 + g3-d2 term */
    if (
        par->correlation_contrib.d2 &&
        par->correlation_contrib.g3
    ){
        result += -(
                3*(3 - fevo1)*(par->Omega0_cdm + par->Omega0_baryon)/2./a2*pow(curlyH1, 2)*f1*(f2 - 1)
                +
                3*(3 - fevo2)*(par->Omega0_cdm + par->Omega0_baryon)/2./a1*pow(curlyH2, 2)*f2*(f1 - 1)
            )
           *(
                coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                /* renormalization term */
                -
                coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        4,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->renormalization,
                    chi1,
                    chi2
                )
            );
    }
    /* g1-g2 + g2-g1 term */
    if (
        par->correlation_contrib.g1 &&
        par->correlation_contrib.g2
    ){
        result += (
                9*pow((par->Omega0_cdm + par->Omega0_baryon), 2)/4./a1/a2*(1 + G1)*(5*s2 - 2)
                +
                9*pow((par->Omega0_cdm + par->Omega0_baryon), 2)/4./a2/a1*(1 + G2)*(5*s1 - 2)
            )
           *(
                coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                /* renormalization term */
                -
                coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        4,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->renormalization,
                    chi1,
                    chi2
                )
            );
    }
    /* g1-g3 + g3-g1 term */
    if (
        par->correlation_contrib.g1 &&
        par->correlation_contrib.g3
    ){
        result += (
                9*pow((par->Omega0_cdm + par->Omega0_baryon), 2)/4./a1/a2*(1 + G1)*(f2 - 1)
                +
                9*pow((par->Omega0_cdm + par->Omega0_baryon), 2)/4./a2/a1*(1 + G2)*(f1 - 1)
            )
           *(
                coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
                /* renormalization term */
                -
                coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        4,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->renormalization,
                    chi1,
                    chi2
                )
            );
    }
    /* g2-g3 + g3-g2 term */
    if (
        par->correlation_contrib.g2 &&
        par->correlation_contrib.g3
    ){
        result += 9*pow((par->Omega0_cdm + par->Omega0_baryon), 2)/4.*(
                (5*s1 - 2)*(f2 - 1)/a1/a2
                +
                (5*s2 - 2)*(f1 - 1)/a2/a1
            )
           *(
                coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sep)
            /* renormalization term */
                -
                coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        4,
                        0,
                        COFFE_INTEGER,
                        COFFE_INTEGER
                    )->renormalization,
                    chi1,
                    chi2
                )
            );
    }
    if (gsl_finite(result)){
#ifdef HAVE_CLASS
        if (
            par->pk_type != COFFE_PK_LINEAR &&
            par->midpoint_approximation
        )
            return result;
        else
#endif
            if (
                par->b_derivative
                ||
                par->f_derivative
            ){
                return
                    result
                   *pow(coffe_interp_spline(&bg->D1, z_mean), 2);
            }
            else if (
                par->b_tilde_derivative
                ||
                par->f_tilde_derivative
            ){
                return
                    result
                   *pow(coffe_interp_spline(&bg->D1, z_mean), 1);
            }
            else{
                if (par->flatsky_local){
                    return
                        result
                       *pow(coffe_interp_spline(&bg->D1, z_mean), 2);
                }
                else{
                    return
                        result
                       *coffe_interp_spline(&bg->D1, z1)
                       *coffe_interp_spline(&bg->D1, z2);
                }
            }
    }
    else{
        fprintf(stderr,
            "ERROR: in function %s, values:\n"
            "mu = %e\n"
            "z_mean = %e\n"
            "chi_mean = %e\n"
            "sep = %e\n"
            "z1 = %e\n"
            "z2 = %e\n"
            "chi1 = %e\n"
            "chi2 = %e\n",
            __func__, mu, z_mean, chi_mean, sep, z1, z2, chi1, chi2);
        exit(EXIT_FAILURE);
    }
}

double functions_single_integrated(
    const coffe_parameters_t *par,
    const coffe_background_t *bg,
    const coffe_integral_array_t *integral,
    const double z_mean,
    const double mu,
    const double sep,
    const double x
)
{
    double result = 0;

    double chi_mean = coffe_interp_spline(&bg->comoving_distance, z_mean);
    const double chi1 = chi_mean - sep*mu/2.;
    const double chi2 = chi_mean + sep*mu/2.;
    const double costheta =
        (2*chi_mean*chi_mean - sep*sep + mu*mu*sep*sep/2.)
       /(2*chi_mean*chi_mean - mu*mu*sep*sep/2.);
    const double lambda1 = chi1*x, lambda2 = chi2*x;

    double r21 = lambda2*lambda2 + chi1*chi1 - 2*chi1*lambda2*costheta;
    double r22 = lambda1*lambda1 + chi2*chi2 - 2*chi2*lambda1*costheta;
    if (r21 < 0) r21 = 0;
    if (r22 < 0) r22 = 0;
    const double z1_const = coffe_interp_spline(&bg->z_as_chi, chi1);
    const double z2_const = coffe_interp_spline(&bg->z_as_chi, chi2);
    const double z1 = coffe_interp_spline(&bg->z_as_chi, lambda1);
    const double z2 = coffe_interp_spline(&bg->z_as_chi, lambda2);

    const double s1 = coffe_interp_spline(&par->magnification_bias1, z1_const);
    const double s2 = coffe_interp_spline(&par->magnification_bias2, z2_const);
    const double sz_mean1 = coffe_interp_spline(&par->magnification_bias1, z_mean);
    const double sz_mean2 = coffe_interp_spline(&par->magnification_bias2, z_mean);
    const double b1 = coffe_interp_spline(&par->galaxy_bias1, z1_const);
    const double b2 = coffe_interp_spline(&par->galaxy_bias2, z2_const);
    const double bz_mean1 = coffe_interp_spline(&par->galaxy_bias1, z_mean);
    const double bz_mean2 = coffe_interp_spline(&par->galaxy_bias2, z_mean);

    double ren1 = 0, ren2 = 0;
    if (par->divergent){
        if (r21 == 0.0) ren1 = coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->renormalization_zero_separation, lambda2);
        else ren1 = coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                    /* renormalization term */
                    -
                    coffe_interp_spline2d(
                        &coffe_find_integral(
                            integral,
                            4,
                            0,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->renormalization,
                        lambda2,
                        chi1
                    );
        if (r22 == 0.0) ren2 = coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->renormalization_zero_separation, lambda1);
        else ren2 = coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                    /* renormalization term */
                    -
                    coffe_interp_spline2d(
                        &coffe_find_integral(
                            integral,
                            4,
                            0,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->renormalization,
                        lambda1,
                        chi2
                    );
    }

    /* den-len + len-den term */
    if (
        par->correlation_contrib.den &&
        par->correlation_contrib.len
    ){
        /* den-len + len-den modified in flatsky */
        if (
            par->flatsky_local_nonlocal &&
            par->output_type != 2
        ){
#ifdef HAVE_CLASS
            if (
                par->pk_type != COFFE_PK_LINEAR &&
                par->midpoint_approximation
            ){

            result +=
               -3 * (par->Omega0_cdm + par->Omega0_baryon) / M_PI / 8
               *(1 + z_mean)
               *fabs(mu)
               *sep
               *(
                    (2 - 5 * sz_mean1) * bz_mean2
                    +
                    (2 - 5 * sz_mean2) * bz_mean1
                )
               *coffe_interp_spline2d(
                    &coffe_find_integral(
                        integral,
                        1,
                        -1,
                        COFFE_HALF_INTEGER,
                        COFFE_HALF_INTEGER
                    )->result2d,
                    z_mean,
                    sep * sqrt(1 - mu * mu)
                );
            }
            else{
#endif
            result +=
               -3 * (par->Omega0_cdm + par->Omega0_baryon) / M_PI / 8
               *coffe_interp_spline(&bg->D1, z_mean)
               *coffe_interp_spline(&bg->D1, z_mean)
               *(1 + z_mean)
               *fabs(mu)
               *sep
               *(
                    (2 - 5 * sz_mean1) * bz_mean2
                    +
                    (2 - 5 * sz_mean2) * bz_mean1
                )
               *coffe_interp_spline(
                    &coffe_find_integral(
                        integral,
                        1,
                        -1,
                        COFFE_HALF_INTEGER,
                        COFFE_HALF_INTEGER
                    )->result,
                    sep * sqrt(1 - mu * mu)
                );
#ifdef HAVE_CLASS
            }
#endif
        }
        else if (!par->flatsky_local_nonlocal){
            if (r21 != 0.0 && r22 != 0.0){
#ifdef HAVE_CLASS
                if (
                    par->pk_type != COFFE_PK_LINEAR &&
                    par->midpoint_approximation
                ){
                    result +=
                       -3
                       *(par->Omega0_cdm + par->Omega0_baryon)
                       /2.
                       *(
                            b1 * (2 - 5 * s2)
                           *chi2
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z2)
                           *(
                                2 * chi1 * costheta
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        1,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1_const + z2) / 2.,
                                    sqrt(r21)
                                )
                                -
                                chi1 * chi1 * lambda2 * (1 - costheta * costheta)
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        0,
                                        2,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1_const + z2) / 2.,
                                    sqrt(r21)
                                )
                               /r21
                            )
                            +
                            b2 * (2 - 5 * s1)
                           *chi1
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z1)
                           *(
                                2 * chi2 * costheta
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        1,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z2_const + z1) / 2.,
                                    sqrt(r22)
                                )
                                -
                                chi2 * chi2 * lambda1
                               *(1 - costheta * costheta)
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        0,
                                        2,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z2_const + z1) / 2.,
                                    sqrt(r22)
                                )
                               /r22
                            )
                        );
                }
                else
#endif
                    result +=
                       -3 * (par->Omega0_cdm + par->Omega0_baryon) / 2.
                       *(
                            b1 * (2 - 5 * s2)
                           *coffe_interp_spline(&bg->D1, z1_const)
                           *chi2
                            /* integrand */
                           *(1 - x)
                           *coffe_interp_spline(&bg->D1, z2)
                           /coffe_interp_spline(&bg->a, z2)
                           *(
                                2 * chi1 * costheta
                               *coffe_interp_spline(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        1,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result,
                                    sqrt(r21)
                                )
                                -
                                chi1 * chi1 * lambda2
                               *(1 - costheta * costheta)
                               *coffe_interp_spline(
                                    &coffe_find_integral(
                                        integral,
                                        0,
                                        2,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result,
                                    sqrt(r21)
                                )
                               /r21
                            )
                            +
                            b2 * (2 - 5 * s1)
                           *coffe_interp_spline(&bg->D1, z2_const)
                           *chi1
                            /* integrand */
                           *(1 - x)
                           *coffe_interp_spline(&bg->D1, z1)
                           /coffe_interp_spline(&bg->a, z1)
                           *(
                                2 * chi2 * costheta
                               *coffe_interp_spline(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        1,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result,
                                    sqrt(r22)
                                )
                                -
                                chi2 * chi2 * lambda1
                               *(1 - costheta * costheta)
                               *coffe_interp_spline(
                                    &coffe_find_integral(
                                        integral,
                                        0,
                                        2,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result,
                                    sqrt(r22)
                                )
                               /r22
                            )
                        );
            }
            else if (r21 == 0.0 && r22 != 0){
#ifdef HAVE_CLASS
                if (par->pk_type != COFFE_PK_LINEAR){
                    result +=
                       -3*(par->Omega0_cdm + par->Omega0_baryon)/2.
                       *(
                            b1 * (2 - 5 * s2)
                           *chi2
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z2)
                           *(
                                2 * chi1
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        1,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1_const + z2) / 2.,
                                    0.0
                                )
                            )
                            +
                            b2 * (2 - 5 * s1)
                           *chi1
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z1)
                           *(
                                2 * chi2 * costheta
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        1,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z2_const + z1) / 2.,
                                    sqrt(r22)
                                )
                                -
                                chi2 * chi2 * lambda1 * (1 - costheta * costheta)
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        0,
                                        2,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z2_const + z1) / 2.,
                                    sqrt(r22)
                                )
                               /r22
                            )
                        );
                }
                else
#endif
                    result +=
                       -3 * (par->Omega0_cdm + par->Omega0_baryon) / 2.
                       *(
                            b1 * (2 - 5 * s2)
                           *chi2
                           *coffe_interp_spline(&bg->D1, z1_const)
                            /* integrand */
                           *(1 - x)
                           *coffe_interp_spline(&bg->D1, z2)
                           /coffe_interp_spline(&bg->a, z2)
                           *(
                                2 * chi1
                               *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)
                            )
                            +
                            b2 * (2 - 5 * s1)
                           *chi1
                           *coffe_interp_spline(&bg->D1, z2_const)
                            /* integrand */
                           *(1 - x)
                           *coffe_interp_spline(&bg->D1, z1)
                           /coffe_interp_spline(&bg->a, z1)
                           *(
                                2 * chi2 * costheta
                               *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                                -
                                chi2 * chi2 * lambda1 * (1 - costheta * costheta)
                               *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                               /r22
                            )
                        );
            }
            else if (r21 != 0 && r22 == 0.0){
#ifdef HAVE_CLASS
                if (
                    par->pk_type != COFFE_PK_LINEAR &&
                    par->midpoint_approximation
                ){
                    result +=
                       -3 * (par->Omega0_cdm + par->Omega0_baryon) / 2.
                       *(
                            b1 * (2 - 5 * s2)
                           *chi2
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z2)
                           *(
                                2 * chi1 * costheta
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        1,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1_const + z2) / 2.,
                                    sqrt(r21)
                                )
                                -
                                chi1 * chi1 * lambda2 * (1 - costheta * costheta)
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        0,
                                        2,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1_const + z2) / 2.,
                                    sqrt(r21)
                                )
                               /r21
                            )
                            +
                            b2 * (2 - 5 * s1)
                           *chi1
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z1)
                           *(
                                2 * chi2
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        0,
                                        2,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z2_const + z1) / 2.,
                                    0.0
                                )
                            )
                        );
                }
                else
#endif
                    result +=
                       -3 * (par->Omega0_cdm + par->Omega0_baryon) / 2.
                       *(
                            b1 * (2 - 5 * s2)
                           *coffe_interp_spline(&bg->D1, z1_const)
                           *chi2
                            /* integrand */
                           *(1 - x)
                           *coffe_interp_spline(&bg->D1, z2)
                           /coffe_interp_spline(&bg->a, z2)
                           *(
                                2 * chi1 * costheta
                               *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                                -
                                chi1 * chi1 * lambda2 * (1 - costheta * costheta)
                               *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                               /r21
                            )
                            +
                            b2 * (2 - 5 * s1)
                           *coffe_interp_spline(&bg->D1, z2_const)
                           *chi1
                            /* integrand */
                           *(1 - x)
                           *coffe_interp_spline(&bg->D1, z1)
                           /coffe_interp_spline(&bg->a, z1)
                           *(
                                2 * chi2
                               *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)
                            )
                        );
            }
            else{
#ifdef HAVE_CLASS
                if (
                    par->pk_type != COFFE_PK_LINEAR &&
                    par->midpoint_approximation
                ){
                    result +=
                       -3 * (par->Omega0_cdm + par->Omega0_baryon) / 2.
                       *(
                            b1 * (2 - 5 * s2)
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z2)
                           *(
                                2 * chi1 * chi2
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        0,
                                        2,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1_const + z2) / 2.,
                                    0.0
                                )
                            )
                            +
                            b2 * (2 - 5 * s1)
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z1)
                           *(
                                2 * chi2 * chi1
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        0,
                                        2,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z2_const + z1) / 2.,
                                    0.0
                                )
                            )
                        );
                }
                else
#endif
                    result +=
                       -3 * (par->Omega0_cdm + par->Omega0_baryon) / 2.
                       *(
                            b1 * (2 - 5 * s2)
                           *coffe_interp_spline(&bg->D1, z1_const)
                            /* integrand */
                           *(1 - x)
                           *coffe_interp_spline(&bg->D1, z2)
                           /coffe_interp_spline(&bg->a, z2)
                           *(
                                2 * chi1 * chi2
                               *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)
                            )
                            +
                            b2 * (2 - 5 * s1)
                           *coffe_interp_spline(&bg->D1, z2_const)
                            /* integrand */
                           *(1 - x)
                           *coffe_interp_spline(&bg->D1, z1)
                           /coffe_interp_spline(&bg->a, z1)
                           *(
                                2 * chi2 * chi1
                               *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)
                            )
                        );
            }
        }
    }
    /* rsd-len + len-rsd term */
    if (
        par->correlation_contrib.rsd &&
        par->correlation_contrib.len
    ){
        /* zero in flatsky */
        if (!par->flatsky_local_nonlocal){
            if (r21 != 0 && r22 != 0){
#ifdef HAVE_CLASS
                if (
                    par->pk_type != COFFE_PK_LINEAR &&
                    par->midpoint_approximation
                ){
                result +=
                    /* constant in front */
                    3 * (par->Omega0_cdm + par->Omega0_baryon) / 2.
                   *(
                        chi2
                       *coffe_interp_spline(&bg->f, z1_const)
                       *(2 - 5 * s2)
                        /* integrand */
                       *(1 - x)
                       /coffe_interp_spline(&bg->a, z2)
                       *(
                            (lambda2 - 6 * chi1 * costheta + 3 * lambda2 * (2 * costheta * costheta - 1))
                           *coffe_interp_spline2d(
                                &coffe_find_integral(
                                    integral,
                                        0,
                                        0,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1_const + z2) / 2.,
                                    sqrt(r21)
                            )
                            /15.
                           -(
                                6 * chi1 * chi1 * chi1 * costheta
                                -
                                chi1 * chi1 * lambda2
                               *(9 * costheta * costheta + 11)
                                +
                                chi1 * lambda2 * lambda2
                               *costheta * (3 * (2 * costheta * costheta - 1) + 19)
                                -
                                2 * lambda2 * lambda2 * lambda2
                               *(3 * (2 * costheta * costheta - 1) + 1)
                            )
                           *coffe_interp_spline2d(
                                &coffe_find_integral(
                                    integral,
                                    0,
                                    2,
                                    COFFE_INTEGER,
                                    COFFE_INTEGER
                                )->result2d,
                                (z1_const + z2) / 2.,
                                sqrt(r21)
                            )
                            /r21
                            /21.
                        )
                        +
                        chi1
                       *coffe_interp_spline(&bg->f, z2_const)
                       *(2 - 5 * s1)
                        /* integrand */
                       *(1 - x)
                       /coffe_interp_spline(&bg->a, z1)
                       *(
                            (lambda1 - 6 * chi2 * costheta + 3 * lambda1 * (2 * costheta * costheta - 1))
                           *coffe_interp_spline2d(
                                &coffe_find_integral(
                                    integral,
                                    0,
                                    0,
                                    COFFE_INTEGER,
                                    COFFE_INTEGER
                                )->result2d,
                                (z2_const + z1) / 2.,
                                sqrt(r22)
                            )
                            /15.
                           -(
                                6 * chi2 * chi2 * chi2 * costheta
                                -
                                chi2 * chi2 * lambda1
                               *(9 * costheta * costheta + 11)
                                +
                                chi2 * lambda1 * lambda1
                               *costheta * (3 * (2 * costheta * costheta - 1) + 19)
                                -
                                2 * lambda1 * lambda1 * lambda1
                               *(3 * (2 * costheta * costheta - 1) + 1)
                            )
                           *coffe_interp_spline2d(
                                &coffe_find_integral(
                                    integral,
                                    0,
                                    2,
                                    COFFE_INTEGER,
                                    COFFE_INTEGER
                                )->result2d,
                                (z2_const + z1) / 2.,
                                sqrt(r22)
                            )
                           /r22
                           /21.
                        )
                    );
                }
                else
#endif
                result +=
                    /* constant in front */
                    3*(par->Omega0_cdm + par->Omega0_baryon)/2.
                   *(
                        chi2*coffe_interp_spline(&bg->f, z1_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                       *(
                            (lambda2 - 6*chi1*costheta + 3*lambda2*(2*costheta*costheta - 1))
                           *coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/15.
                           -(
                                6*chi1*chi1*chi1*costheta - chi1*chi1*lambda2*(9*costheta*costheta + 11)
                               +chi1*lambda2*lambda2*costheta*(3*(2*costheta*costheta - 1) + 19)
                               -2*lambda2*lambda2*lambda2*(3*(2*costheta*costheta - 1) + 1)
                            )
                           *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/r21/21.
                        )
                        +
                        chi1*coffe_interp_spline(&bg->f, z2_const)*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                       *(
                            (lambda1 - 6*chi2*costheta + 3*lambda1*(2*costheta*costheta - 1))
                           *coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/15.
                           -(
                                6*chi2*chi2*chi2*costheta - chi2*chi2*lambda1*(9*costheta*costheta + 11)
                               +chi2*lambda1*lambda1*costheta*(3*(2*costheta*costheta - 1) + 19)
                               -2*lambda1*lambda1*lambda1*(3*(2*costheta*costheta - 1) + 1)
                            )
                           *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/r22/21.
                        )
                    );
                if (fabs(mu) < 0.999){
#ifdef HAVE_CLASS
                    if (
                        par->pk_type != COFFE_PK_LINEAR &&
                        par->midpoint_approximation
                    ){
                    result +=
                    3 * (par->Omega0_cdm + par->Omega0_baryon) / 2.
                   *(
                        chi2
                       *coffe_interp_spline(&bg->f, z1_const)
                       *(2 - 5 * s2)
                        /* integrand */
                       *(1 - x)
                       /coffe_interp_spline(&bg->a, z2)
                       *(
                           -(
                               -4 * pow(chi1, 5) * costheta
                                -
                                pow(chi1, 3)*pow(lambda2, 2)
                               *costheta * ((2 * costheta * costheta - 1) + 7)
                                +
                                pow(chi1, 2) * pow(lambda2, 3)
                               *(pow(costheta, 4) + 12 * costheta * costheta - 21)
                                -
                                3 * chi1 * pow(lambda2, 4)
                               *costheta * ((2 * costheta * costheta - 1) - 5)
                                -
                                pow(lambda2, 5)
                               *(3 * (2 * costheta * costheta - 1) + 1)
                                +
                                12 * pow(chi1, 4) * lambda2
                            )
                           *coffe_interp_spline2d(
                                &coffe_find_integral(
                                    integral,
                                    0,
                                    4,
                                    COFFE_INTEGER,
                                    COFFE_INTEGER
                                )->result2d,
                                (z1_const + z2) / 2.,
                                sqrt(r21)
                            )
                           /r21
                           /r21
                           /35.
                        )
                        +
                        chi1
                       *coffe_interp_spline(&bg->f, z2_const)
                       *(2 - 5 * s1)
                        /* integrand */
                       *(1 - x)
                       /coffe_interp_spline(&bg->a, z1)
                       *(
                           -(
                               -4 * pow(chi2, 5) * costheta
                                -
                                pow(chi2, 3) * pow(lambda1, 2)
                               *costheta * ((2 * costheta * costheta - 1) + 7)
                                +
                                pow(chi2, 2) * pow(lambda1, 3)
                               *(pow(costheta, 4) + 12 * costheta * costheta - 21)
                                -
                                3 * chi2 * pow(lambda1, 4)
                               *costheta * ((2 * costheta * costheta - 1) - 5)
                                -
                                pow(lambda1, 5)
                               *(3 * (2 * costheta * costheta - 1) + 1)
                                +
                                12 * pow(chi2, 4) * lambda1
                            )
                           *coffe_interp_spline2d(
                                &coffe_find_integral(
                                    integral,
                                    0,
                                    4,
                                    COFFE_INTEGER,
                                    COFFE_INTEGER
                                )->result2d,
                                (z2_const + z1) / 2.,
                                sqrt(r22)
                            )
                           /r22
                           /r22
                           /35.
                        )
                    );
                    }
                    else
#endif
                    result +=
                    3*(par->Omega0_cdm + par->Omega0_baryon)/2.
                   *(
                        chi2*coffe_interp_spline(&bg->f, z1_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                       *(
                           -(
                               -4*pow(chi1, 5)*costheta
                               -pow(chi1, 3)*pow(lambda2, 2)*costheta*((2*costheta*costheta - 1) + 7)
                               +pow(chi1, 2)*pow(lambda2, 3)*(pow(costheta, 4) + 12*costheta*costheta - 21)
                               -3*chi1*pow(lambda2, 4)*costheta*((2*costheta*costheta - 1) - 5)
                               -pow(lambda2, 5)*(3*(2*costheta*costheta - 1) + 1)
                               +12*pow(chi1, 4)*lambda2
                            )
                           *coffe_interp_spline(&coffe_find_integral(integral, 0, 4, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/r21/r21/35.
                        )
                        +
                        chi1*coffe_interp_spline(&bg->f, z2_const)*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                       *(
                           -(
                               -4*pow(chi2, 5)*costheta
                               -pow(chi2, 3)*pow(lambda1, 2)*costheta*((2*costheta*costheta - 1) + 7)
                               +pow(chi2, 2)*pow(lambda1, 3)*(pow(costheta, 4) + 12*costheta*costheta - 21)
                               -3*chi2*pow(lambda1, 4)*costheta*((2*costheta*costheta - 1) - 5)
                               -pow(lambda1, 5)*(3*(2*costheta*costheta - 1) + 1)
                               +12*pow(chi2, 4)*lambda1
                            )
                           *coffe_interp_spline(&coffe_find_integral(integral, 0, 4, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/r22/r22/35.
                        )
                    );
                }
                else{
#ifdef HAVE_CLASS
                    if (
                        par->pk_type != COFFE_PK_LINEAR &&
                        par->midpoint_approximation
                    ){
                    result +=
                        3 * (par->Omega0_cdm + par->Omega0_baryon) / 2.
                       *(
                            chi2
                           *coffe_interp_spline(&bg->f, z1_const)
                           *(2 - 5 * s2)
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z2)
                           *4. * (lambda2 + chi1)
                           *coffe_interp_spline2d(
                                &coffe_find_integral(
                                    integral,
                                    0,
                                    4,
                                    COFFE_INTEGER,
                                    COFFE_INTEGER
                                )->result2d,
                                (z1_const + z2) / 2.,
                                sqrt(r21)
                            )
                           /35.
                            +
                            chi1
                           *coffe_interp_spline(&bg->f, z2_const)
                           *(2 - 5 * s1)
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z1)
                           *4. * (lambda1 + chi2)
                           *coffe_interp_spline2d(
                                &coffe_find_integral(
                                    integral,
                                    0,
                                    4,
                                    COFFE_INTEGER,
                                    COFFE_INTEGER
                                )->result2d,
                                (z2_const + z1) / 2.,
                                sqrt(r22)
                            )
                           /35.
                        );
                    }
                    else
#endif
                    result +=
                    3*(par->Omega0_cdm + par->Omega0_baryon)/2.
                   *(
                        chi2*coffe_interp_spline(&bg->f, z1_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                       *4.*(lambda2 + chi1)*coffe_interp_spline(&coffe_find_integral(integral, 0, 4, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/35.
                        +
                        chi1*coffe_interp_spline(&bg->f, z2_const)*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                       *4.*(lambda1 + chi2)*coffe_interp_spline(&coffe_find_integral(integral, 0, 4, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/35.
                    );
                }
            }
            else if (r21 == 0 && r22 != 0){
                result +=
                    3*(par->Omega0_cdm + par->Omega0_baryon)/2.
                   *(
                        chi2*coffe_interp_spline(&bg->f, z1_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
                        /* integrand */
                       *(1 - chi1/chi2)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                       *(
                           -2*chi1*coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)/15.
                        )
                        +
                        chi1*coffe_interp_spline(&bg->f, z2_const)*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                       *(
                            (lambda1 - 6*chi2*costheta + 3*lambda1*(2*costheta*costheta - 1))
                           *coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/15.
                           -(
                                6*chi2*chi2*chi2*costheta - chi2*chi2*lambda1*(9*costheta*costheta + 11)
                               +chi2*lambda1*lambda1*costheta*(3*(2*costheta*costheta - 1) + 19)
                               -2*lambda1*lambda1*lambda1*(3*(2*costheta*costheta - 1) + 1)
                            )
                           *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/r22/21.
                           -(
                               -4*pow(chi2, 5)*costheta
                               -pow(chi2, 3)*pow(lambda1, 2)*costheta*((2*costheta*costheta - 1) + 7)
                               +pow(chi2, 2)*pow(lambda1, 3)*(pow(costheta, 4) + 12*costheta*costheta - 21)
                               -3*chi2*pow(lambda1, 4)*costheta*((2*costheta*costheta - 1) - 5)
                               -pow(lambda1, 5)*(3*(2*costheta*costheta - 1) + 1)
                               +12*pow(chi2, 4)*lambda1
                            )
                           *coffe_interp_spline(&coffe_find_integral(integral, 0, 4, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/r22/r22/35.
                        )
                    );
            }
            else if (r21 != 0 && r22 == 0){
                result +=
                    /* constant in front */
                    3*(par->Omega0_cdm + par->Omega0_baryon)/2.
                   *(
                        chi2*coffe_interp_spline(&bg->f, z1_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                       *(
                            (lambda2 - 6*chi1*costheta + 3*lambda2*(2*costheta*costheta - 1))
                           *coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/15.
                           -(
                                6*chi1*chi1*chi1*costheta - chi1*chi1*lambda2*(9*costheta*costheta + 11)
                               +chi1*lambda2*lambda2*costheta*(3*(2*costheta*costheta - 1) + 19)
                               -2*lambda2*lambda2*lambda2*(3*(2*costheta*costheta - 1) + 1)
                            )
                           *coffe_interp_spline(&coffe_find_integral(integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/r21/21.
                           -(
                               -4*pow(chi1, 5)*costheta
                               -pow(chi1, 3)*pow(lambda2, 2)*costheta*((2*costheta*costheta - 1) + 7)
                               +pow(chi1, 2)*pow(lambda2, 3)*(pow(costheta, 4) + 12*costheta*costheta - 21)
                               -3*chi1*pow(lambda2, 4)*costheta*((2*costheta*costheta - 1) - 5)
                               -pow(lambda2, 5)*(3*(2*costheta*costheta - 1) + 1)
                               +12*pow(chi1, 4)*lambda2
                            )
                           *coffe_interp_spline(&coffe_find_integral(integral, 0, 4, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/r21/r21/35.
                        )
                        +
                        chi1*coffe_interp_spline(&bg->f, z2_const)*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
                        /* integrand */
                       *(1 - chi2/chi1)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                       *(
                          -2*chi2*coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)/15.
                        )
                    );
            }
            else{
                result +=
                    3*(par->Omega0_cdm + par->Omega0_baryon)/2.
                   *(
                        chi2*coffe_interp_spline(&bg->f, z1_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
                        /* integrand */
                       *(1 - chi1/chi2)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                       *(
                          -2*chi1*coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)/15.
                        )
                        +
                        chi1*coffe_interp_spline(&bg->f, z2_const)*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
                        /* integrand */
                       *(1 - chi2/chi1)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                       *(
                           -2*chi2*coffe_interp_spline(&coffe_find_integral(integral, 0, 0, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)/15.
                        )
                    );
            }
        }
    }
    /* d1-len + len-d1 term */
    if (
        par->correlation_contrib.d1 &&
        par->correlation_contrib.len
    ){
        /* zero in flatsky */
        if (!par->flatsky_local_nonlocal){
            if (r21 != 0 && r22 != 0){
#if HAVE_CLASS
                if (
                    par->pk_type != COFFE_PK_LINEAR &&
                    par->midpoint_approximation
                ){
                    result +=
                        /* constant in front */
                        3 * (par->Omega0_cdm + par->Omega0_baryon) / 2.
                       *(
                            chi2
                           *coffe_interp_spline(&bg->conformal_Hz, z1_const)
                           *coffe_interp_spline(&bg->f, z1_const)
                           *coffe_interp_spline(&bg->G1, z1_const)
                           *(2 - 5 * s2)
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z2)
                           *(
                                2 * (costheta * (lambda2 * lambda2 - 2 * chi1 * chi1) + chi1 * lambda2 * (2 * (2 * costheta * costheta - 1) - 1))
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        1,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1_const + z2) / 2.,
                                    sqrt(r21))/15.
                                +
                                2 * costheta
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        2,
                                        0,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1_const + z2) / 2.,
                                    sqrt(r21)) / 3.
                                -
                                (
                                    4 * pow(chi1, 4) * costheta
                                    -
                                    pow(chi1, 3) * lambda2 * (costheta * costheta + 9)
                                    +
                                    chi1 * chi1 * lambda2 * lambda2 * costheta * (costheta * costheta + 5)
                                    -
                                    2 * chi1 * pow(lambda2, 3) * ((2 * costheta * costheta - 1) - 2)
                                    -
                                    2 * pow(lambda2, 4) * costheta
                                )
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        3,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1_const + z2) / 2.,
                                    sqrt(r21))
                               /r21
                               /15.
                            )
                            +
                            chi1
                           *coffe_interp_spline(&bg->conformal_Hz, z2_const)
                           *coffe_interp_spline(&bg->f, z2_const)
                           *coffe_interp_spline(&bg->G1, z2_const)
                           *(2 - 5 * s2)
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z1)
                           *(
                                2 * (costheta * (lambda1 * lambda1 - 2 * chi2 * chi2) + chi2 * lambda1 * (2 * (2 * costheta * costheta - 1) - 1))
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        1,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1 + z2_const) / 2.,
                                    sqrt(r22)
                                )
                               /15.
                                +
                                2 * costheta
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        2,
                                        0,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1 + z2_const) / 2.,
                                    sqrt(r22)
                                )
                               /3.
                                -
                                (
                                    4 * pow(chi2, 4) * costheta
                                    -
                                    pow(chi2, 3) * lambda1 * (costheta * costheta + 9)
                                    +
                                    chi2 * chi2 * lambda1 * lambda1 * costheta * (costheta * costheta + 5)
                                    -
                                    2 * chi2 * pow(lambda1, 3) * ((2 * costheta * costheta - 1) - 2)
                                    -
                                    2 * pow(lambda1, 4) * costheta
                                )
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        3,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1 + z2_const) / 2.,
                                    sqrt(r22))
                               /r22
                               /15.
                            )
                        );
                }
                else
#endif
                result +=
                    /* constant in front */
                    3*(par->Omega0_cdm + par->Omega0_baryon)/2.
                   *(
                        chi2*coffe_interp_spline(&bg->conformal_Hz, z1_const)*coffe_interp_spline(&bg->f, z1_const)
                       *coffe_interp_spline(&bg->G1, z1_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                       *(
                            2*(costheta*(lambda2*lambda2 - 2*chi1*chi1) + chi1*lambda2*(2*(2*costheta*costheta - 1) - 1))
                           *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/15.
                           +2*costheta*coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/3.
                           -(
                                4*pow(chi1, 4)*costheta
                               -pow(chi1, 3)*lambda2*(costheta*costheta + 9)
                               +chi1*chi1*lambda2*lambda2*costheta*(costheta*costheta + 5)
                               -2*chi1*pow(lambda2, 3)*((2*costheta*costheta - 1) - 2)
                               -2*pow(lambda2, 4)*costheta
                            )*coffe_interp_spline(&coffe_find_integral(integral, 1, 3, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/r21/15.
                        )
                        +
                        chi1*coffe_interp_spline(&bg->conformal_Hz, z2_const)*coffe_interp_spline(&bg->f, z2_const)
                       *coffe_interp_spline(&bg->G1, z2_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z2_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                       *(
                            2*(costheta*(lambda1*lambda1 - 2*chi2*chi2) + chi2*lambda1*(2*(2*costheta*costheta - 1) - 1))
                           *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/15.
                           +2*costheta*coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/3.
                           -(
                                4*pow(chi2, 4)*costheta
                               -pow(chi2, 3)*lambda1*(costheta*costheta + 9)
                               +chi2*chi2*lambda1*lambda1*costheta*(costheta*costheta + 5)
                               -2*chi2*pow(lambda1, 3)*((2*costheta*costheta - 1) - 2)
                               -2*pow(lambda1, 4)*costheta
                            )*coffe_interp_spline(&coffe_find_integral(integral, 1, 3, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/r22/15.
                        )
                    );
            }
            else if (r21 == 0 && r22 != 0){
#if HAVE_CLASS
                if (
                    par->pk_type != COFFE_PK_LINEAR &&
                    par->midpoint_approximation
                ){
                    result +=
                        /* constant in front */
                        3*(par->Omega0_cdm + par->Omega0_baryon)/2.
                       *(
                            chi2
                           *coffe_interp_spline(&bg->conformal_Hz, z1_const)
                           *coffe_interp_spline(&bg->f, z1_const)
                           *coffe_interp_spline(&bg->G1, z1_const)
                           *(2 - 5 * s2)
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z2)
                           *(
                               2
                               *coffe_interp_spline2d(
                                        &coffe_find_integral(
                                        integral,
                                            2,
                                            0,
                                            COFFE_INTEGER,
                                            COFFE_INTEGER
                                        )->result2d,
                                        (z1_const + z2) / 2.,
                                        0.0
                                )
                               /3.
                            )
                            +
                            chi1
                           *coffe_interp_spline(&bg->conformal_Hz, z2_const)
                           *coffe_interp_spline(&bg->f, z2_const)
                           *coffe_interp_spline(&bg->G1, z2_const)
                           *(2 - 5 * s2)
                            /* integrand */
                           *(1 - x)
                           /coffe_interp_spline(&bg->a, z1)
                           *(
                                2 * (costheta * (lambda1 * lambda1 - 2 * chi2 * chi2) + chi2 * lambda1 * (2 * (2 * costheta * costheta - 1) - 1))
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        1,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1 + z2_const) / 2.,
                                    sqrt(r22)
                                )
                               /15.
                                +
                                2 * costheta
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        2,
                                        0,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1 + z2_const) / 2.,
                                    sqrt(r22)
                                )
                               /3.
                                -
                                (
                                    4 * pow(chi2, 4) * costheta
                                    -
                                    pow(chi2, 3) * lambda1 * (costheta * costheta + 9)
                                    +
                                    chi2 * chi2 * lambda1 * lambda1 * costheta * (costheta * costheta + 5)
                                    -
                                    2 * chi2 * pow(lambda1, 3) * ((2 * costheta * costheta - 1) - 2)
                                    -
                                    2 * pow(lambda1, 4) * costheta
                                )
                               *coffe_interp_spline2d(
                                    &coffe_find_integral(
                                        integral,
                                        1,
                                        3,
                                        COFFE_INTEGER,
                                        COFFE_INTEGER
                                    )->result2d,
                                    (z1 + z2_const) / 2.,
                                    sqrt(r22)
                                )
                               /r22
                               /15.
                            )
                        );
                }
                else
#endif
                result +=
                    /* constant in front */
                    3*(par->Omega0_cdm + par->Omega0_baryon)/2.
                   *(
                        chi2*coffe_interp_spline(&bg->conformal_Hz, z1_const)*coffe_interp_spline(&bg->f, z1_const)
                       *coffe_interp_spline(&bg->G1, z1_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                       *(
                           2*coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)/3.
                        )
                        +
                        chi1*coffe_interp_spline(&bg->conformal_Hz, z2_const)*coffe_interp_spline(&bg->f, z2_const)
                       *coffe_interp_spline(&bg->G1, z2_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z2_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                       *(
                            2*(costheta*(lambda1*lambda1 - 2*chi2*chi2) + chi2*lambda1*(2*(2*costheta*costheta - 1) - 1))
                           *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/15.
                           +2*costheta*coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/3.
                           -(
                                4*pow(chi2, 4)*costheta
                               -pow(chi2, 3)*lambda1*(costheta*costheta + 9)
                               +chi2*chi2*lambda1*lambda1*costheta*(costheta*costheta + 5)
                               -2*chi2*pow(lambda1, 3)*((2*costheta*costheta - 1) - 2)
                               -2*pow(lambda1, 4)*costheta
                            )*coffe_interp_spline(&coffe_find_integral(integral, 1, 3, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/r22/15.
                        )
                    );
            }
            else if (r21 != 0 && r22 == 0){
#if HAVE_CLASS
                if (
                    par->pk_type != COFFE_PK_LINEAR &&
                    par->midpoint_approximation
                ){
                result +=
                    /* constant in front */
                    3 * (par->Omega0_cdm + par->Omega0_baryon) / 2.
                   *(
                        chi2
                       *coffe_interp_spline(&bg->conformal_Hz, z1_const)
                       *coffe_interp_spline(&bg->f, z1_const)
                       *coffe_interp_spline(&bg->G1, z1_const)
                       *(2 - 5 * s2)
                        /* integrand */
                       *(1 - x)
                       /coffe_interp_spline(&bg->a, z2)
                       *(
                            2 * (costheta * (lambda2 * lambda2 - 2 * chi1 * chi1) + chi1 * lambda2 * (2 * (2 * costheta * costheta - 1) - 1))
                           *coffe_interp_spline2d(
                               &coffe_find_integral(
                                   integral,
                                   1,
                                   1,
                                   COFFE_INTEGER,
                                   COFFE_INTEGER
                                   )->result2d,
                               (z1_const + z2) / 2.,
                               sqrt(r21)
                               )
                           /15.
                            +
                            2 * costheta
                           *coffe_interp_spline2d(
                               &coffe_find_integral(
                                   integral,
                                   2,
                                   0,
                                   COFFE_INTEGER,
                                   COFFE_INTEGER
                                   )->result2d,
                               (z1_const + z2) / 2.,
                               sqrt(r21)
                               )
                           /3.
                           -(
                                4 * pow(chi1, 4) * costheta
                               -pow(chi1, 3) * lambda2 * (costheta * costheta + 9)
                               +chi1 * chi1 * lambda2 * lambda2 * costheta * (costheta * costheta + 5)
                               -2 * chi1 * pow(lambda2, 3) * ((2 * costheta * costheta - 1) - 2)
                               -2 * pow(lambda2, 4) * costheta
                            )
                           *coffe_interp_spline2d(
                                   &coffe_find_integral(
                                       integral,
                                       1,
                                       3,
                                       COFFE_INTEGER,
                                       COFFE_INTEGER)->result2d,
                                   (z1_const + z2) / 2.,
                                   sqrt(r21)
                                   )
                           /r21
                           /15.
                        )
                        +
                        chi1
                       *coffe_interp_spline(&bg->conformal_Hz, z2_const)
                       *coffe_interp_spline(&bg->f, z2_const)
                       *coffe_interp_spline(&bg->G1, z2_const)
                       *(2 - 5 * s2)
                        /* integrand */
                       *(1 - x)
                       /coffe_interp_spline(&bg->a, z1)
                       *(
                           2
                          *coffe_interp_spline2d(
                               &coffe_find_integral(
                                   integral,
                                   2,
                                   0,
                                   COFFE_INTEGER,
                                   COFFE_INTEGER
                                )->result2d,
                                (z1 + z2_const) / 2.,
                                0.0
                           )
                          /3.
                        )
                    );
                }
                else
#endif
                result +=
                    /* constant in front */
                    3*(par->Omega0_cdm + par->Omega0_baryon)/2.
                   *(
                        chi2*coffe_interp_spline(&bg->conformal_Hz, z1_const)*coffe_interp_spline(&bg->f, z1_const)
                       *coffe_interp_spline(&bg->G1, z1_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                       *(
                            2*(costheta*(lambda2*lambda2 - 2*chi1*chi1) + chi1*lambda2*(2*(2*costheta*costheta - 1) - 1))
                           *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/15.
                           +2*costheta*coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/3.
                           -(
                                4*pow(chi1, 4)*costheta
                               -pow(chi1, 3)*lambda2*(costheta*costheta + 9)
                               +chi1*chi1*lambda2*lambda2*costheta*(costheta*costheta + 5)
                               -2*chi1*pow(lambda2, 3)*((2*costheta*costheta - 1) - 2)
                               -2*pow(lambda2, 4)*costheta
                            )*coffe_interp_spline(&coffe_find_integral(integral, 1, 3, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/r21/15.
                        )
                        +
                        chi1*coffe_interp_spline(&bg->conformal_Hz, z2_const)*coffe_interp_spline(&bg->f, z2_const)
                       *coffe_interp_spline(&bg->G1, z2_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z2_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                       *(
                           2*coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)/3.
                        )
                    );
            }
            else{
#if HAVE_CLASS
                if (
                    par->pk_type != COFFE_PK_LINEAR &&
                    par->midpoint_approximation
                ){
                result +=
                    /* constant in front */
                    3 * (par->Omega0_cdm + par->Omega0_baryon) / 2.
                   *(
                        chi2
                       *coffe_interp_spline(&bg->conformal_Hz, z1_const)
                       *coffe_interp_spline(&bg->f, z1_const)
                       *coffe_interp_spline(&bg->G1, z1_const)
                       *(2 - 5 * s2)
                        /* integrand */
                       *(1 - x)
                       /coffe_interp_spline(&bg->a, z2)
                       *(
                           2
                          *coffe_interp_spline2d(
                              &coffe_find_integral(
                                  integral,
                                  2,
                                  0,
                                  COFFE_INTEGER,
                                  COFFE_INTEGER)->result2d,
                              (z1_const + z2) / 2.,
                              0.0
                              )
                          /3.
                        )
                        +
                        chi1
                       *coffe_interp_spline(&bg->conformal_Hz, z2_const)
                       *coffe_interp_spline(&bg->f, z2_const)
                       *coffe_interp_spline(&bg->G1, z2_const)
                       *(2 - 5 * s2)
                        /* integrand */
                       *(1 - x)
                       /coffe_interp_spline(&bg->a, z1)
                       *(
                           2
                          *coffe_interp_spline2d(
                              &coffe_find_integral(
                                    integral,
                                    2,
                                    0,
                                    COFFE_INTEGER,
                                    COFFE_INTEGER
                                )->result2d,
                                (z1 + z2_const) / 2.,
                                0.0
                              )
                          /3.
                        )
                    );
                }
                else
#endif
                result +=
                    /* constant in front */
                    3*(par->Omega0_cdm + par->Omega0_baryon)/2.
                   *(
                        chi2*coffe_interp_spline(&bg->conformal_Hz, z1_const)*coffe_interp_spline(&bg->f, z1_const)
                       *coffe_interp_spline(&bg->G1, z1_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                       *(
                           2*coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)/3.
                        )
                        +
                        chi1*coffe_interp_spline(&bg->conformal_Hz, z2_const)*coffe_interp_spline(&bg->f, z2_const)
                       *coffe_interp_spline(&bg->G1, z2_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z2_const)
                        /* integrand */
                       *(1 - x)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                       *(
                           2*coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)/3.
                        )
                    );
            }
        }
    }
    /* d2-len + len-d2 term */
    if (
        par->correlation_contrib.d2 &&
        par->correlation_contrib.len
    ){
        result +=
            /* constant in front */
           -3*(par->Omega0_cdm + par->Omega0_baryon)/2.
           *(
                chi2*(3 - coffe_interp_spline(&par->evolution_bias1, z1_const))*coffe_interp_spline(&bg->f, z1_const)
               *pow(coffe_interp_spline(&bg->conformal_Hz, z1_const), 2)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
               *(
                    /* integrand */
                    (1 - x)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                   *(
                        2*chi1*costheta*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                       -chi1*chi1*lambda2*(1 - costheta*costheta)*coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                    )
                )
                +
                chi1*(3 - coffe_interp_spline(&par->evolution_bias2, z2_const))*coffe_interp_spline(&bg->f, z2_const)
               *pow(coffe_interp_spline(&bg->conformal_Hz, z2_const), 2)*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
               *(
                    /* integrand */
                    (1 - x)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                   *(
                        2*chi2*costheta*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                       -chi2*chi2*lambda1*(1 - costheta*costheta)*coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                    )
                )
            );
    }
    /* g1-len + len-g1 term */
    if (
        par->correlation_contrib.g1 &&
        par->correlation_contrib.len
    ){
        result +=
            /* constant in front */
            9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/4.
           *(
                chi2*(1 + coffe_interp_spline(&bg->G1, z1_const))*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
               *(
                    /* integrand */
                    (1 - x)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                   *(
                        2*chi1*costheta*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                       -chi1*chi1*lambda2*(1 - costheta*costheta)*coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                    )
                )
                +
                chi1*(1 + coffe_interp_spline(&bg->G2, z2_const))*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
               *(
                    /* integrand */
                    (1 - x)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                   *(
                        2*chi2*costheta*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                       -chi2*chi2*lambda1*(1 - costheta*costheta)*coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                    )
                )
            );
    }
    /* g2-len + len-g2 term */
    if (
        par->correlation_contrib.g2 &&
        par->correlation_contrib.len
    ){
        result +=
            /* constant in front */
            9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/4.
           *(
                chi2*(5*s1 - 2)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
               *(
                    /* integrand */
                    (1 - x)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                   *(
                        2*chi1*costheta*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                       -chi1*chi1*lambda2*(1 - costheta*costheta)*coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                    )
                )
                +
                chi1*(5*s2 - 2)*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
               *(
                    /* integrand */
                    (1 - x)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                   *(
                        2*chi2*costheta*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                       -chi2*chi2*lambda1*(1 - costheta*costheta)*coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                    )
                )
            );
    }
    /* g3-len + len-g3 term */
    if (
        par->correlation_contrib.g3 &&
        par->correlation_contrib.len
    ){
        result +=
            /* constant in front */
            9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/4.
           *(
                chi2*(coffe_interp_spline(&bg->f, z1_const) - 1)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
               *(
                    /* integrand */
                    (1 - x)*coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
                   *(
                        2*chi1*costheta*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                       -chi1*chi1*lambda2*(1 - costheta*costheta)*coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                    )
                )
                +
                chi1*(coffe_interp_spline(&bg->f, z2_const) - 1)*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
               *(
                    /* integrand */
                    (1 - x)*coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
                   *(
                        2*chi2*costheta*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                       -chi2*chi2*lambda1*(1 - costheta*costheta)*coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                    )
                )
            );
    }
    /* den-g4 + g4-den term */
    if (
        par->correlation_contrib.den &&
        par->correlation_contrib.g4
    ){
        result +=
            /* constant in front */
           -3*(par->Omega0_cdm + par->Omega0_baryon)
           *(
                b1*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
                /* integrand */
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
               *coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
               +
                b2*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
                /* integrand */
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
               *coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
            );
    }
    /* den-g5 + g5-den term */
    if (
        par->correlation_contrib.den &&
        par->correlation_contrib.g5
    ){
        result +=
            /* constant in front */
           -3*(par->Omega0_cdm + par->Omega0_baryon)
           *(
                chi2*b1*coffe_interp_spline(&bg->G2, z2_const)*coffe_interp_spline(&bg->D1, z1_const)
                /* integrand */
               *coffe_interp_spline(&bg->conformal_Hz, z2)*(coffe_interp_spline(&bg->f, z2) - 1)
               *coffe_interp_spline(&bg->D1, z2)*coffe_interp_spline(&bg->a, z2)
               *coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
               +
                chi1*b2*coffe_interp_spline(&bg->G1, z1_const)*coffe_interp_spline(&bg->D1, z2_const)
                /* integrand */
               *coffe_interp_spline(&bg->conformal_Hz, z1)*(coffe_interp_spline(&bg->f, z1) - 1)
               *coffe_interp_spline(&bg->D1, z1)*coffe_interp_spline(&bg->a, z1)
               *coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
            );
    }
    /* rsd-g4 + g4-rsd term */
    if (
        par->correlation_contrib.rsd &&
        par->correlation_contrib.g4
    ){
        result +=
            3*(par->Omega0_cdm + par->Omega0_baryon)
           *(
                coffe_interp_spline(&bg->f, z1_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
                /* integrand */
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
               *(
                    (2*r21/3. + (costheta*costheta - 1)*lambda2*lambda2)
                   *coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                   -coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/3.
                )
               +
                coffe_interp_spline(&bg->f, z2_const)*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
                /* integrand */
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
               *(
                    (2*r22/3. + (costheta*costheta - 1)*lambda1*lambda1)
                   *coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                   -coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/3.
                )
            );
    }
    /* rsd-g5 + g5-rsd term */
    if (
        par->correlation_contrib.rsd &&
        par->correlation_contrib.g5
    ){
        result +=
            3*(par->Omega0_cdm + par->Omega0_baryon)
           *(
                chi2*coffe_interp_spline(&bg->f, z1_const)*coffe_interp_spline(&bg->G2, z2_const)*coffe_interp_spline(&bg->D1, z1_const)
                /* integrand */
               *coffe_interp_spline(&bg->conformal_Hz, z2)*(coffe_interp_spline(&bg->f, z2) - 1)
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
               *(
                    (2*r21/3. + (costheta*costheta - 1)*lambda2*lambda2)
                   *coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
                   -coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))/3.
                )
               +
                chi1*coffe_interp_spline(&bg->f, z2_const)*coffe_interp_spline(&bg->G1, z1_const)*coffe_interp_spline(&bg->D1, z2_const)
                /* integrand */
               *coffe_interp_spline(&bg->conformal_Hz, z1)*(coffe_interp_spline(&bg->f, z1) - 1)
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
               *(
                    (2*r22/3. + (costheta*costheta - 1)*lambda1*lambda1)
                   *coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
                   -coffe_interp_spline(&coffe_find_integral(integral, 2, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))/3.
                )
            );
    }
    /* d1-g4 + d1-g4 term */
    if (
        par->correlation_contrib.d1 &&
        par->correlation_contrib.g4
    ){
        result +=
            3*(par->Omega0_cdm + par->Omega0_baryon)
           *(
                coffe_interp_spline(&bg->conformal_Hz, z1_const)*coffe_interp_spline(&bg->f, z1_const)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)*(lambda2*costheta - chi1)
               *coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
               +
                coffe_interp_spline(&bg->conformal_Hz, z2_const)*coffe_interp_spline(&bg->f, z2_const)*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)*(lambda1*costheta - chi2)
               *coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
            );
    }
    /* d1-g5 + d1-g5 term */
    if (
        par->correlation_contrib.d1 &&
        par->correlation_contrib.g5
    ){
        result +=
            3*(par->Omega0_cdm + par->Omega0_baryon)
           *(
                chi2*coffe_interp_spline(&bg->conformal_Hz, z1_const)*coffe_interp_spline(&bg->f, z1_const)
               *coffe_interp_spline(&bg->G2, z2_const)*coffe_interp_spline(&bg->D1, z1_const)
               *coffe_interp_spline(&bg->conformal_Hz, z2)*(coffe_interp_spline(&bg->f, z2) - 1)
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)*(lambda2*costheta - chi1)
               *coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r21))
               +
                chi1*coffe_interp_spline(&bg->conformal_Hz, z2_const)*coffe_interp_spline(&bg->f, z2_const)
               *coffe_interp_spline(&bg->G1, z1_const)*coffe_interp_spline(&bg->D1, z2_const)
               *coffe_interp_spline(&bg->conformal_Hz, z1)*(coffe_interp_spline(&bg->f, z1) - 1)
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)*(lambda1*costheta - chi2)
               *coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r22))
            );
    }
    /* d2-g4 + g4-d2 term */
    if (
        par->correlation_contrib.d2 &&
        par->correlation_contrib.g4
    ){
        result +=
           -3*(par->Omega0_cdm + par->Omega0_baryon)
           *(
                (3 - coffe_interp_spline(&par->evolution_bias1, z1_const))*coffe_interp_spline(&bg->f, z1_const)
               *pow(coffe_interp_spline(&bg->conformal_Hz, z1_const), 2)*(2 - 5*s2)*coffe_interp_spline(&bg->D1, z1_const)
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
               *ren1
               +
                (3 - coffe_interp_spline(&par->evolution_bias2, z2_const))*coffe_interp_spline(&bg->f, z2_const)
               *pow(coffe_interp_spline(&bg->conformal_Hz, z2_const), 2)*(2 - 5*s1)*coffe_interp_spline(&bg->D1, z2_const)
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
               *ren2
            );
    }
    /* d2-g5 + g5-d2 term */
    if (
        par->correlation_contrib.d2 &&
        par->correlation_contrib.g5
    ){
        result +=
           -3*(par->Omega0_cdm + par->Omega0_baryon)
           *(
                chi2*(3 - coffe_interp_spline(&par->evolution_bias1, z1_const))*coffe_interp_spline(&bg->f, z1_const)
               *pow(coffe_interp_spline(&bg->conformal_Hz, z1_const), 2)*coffe_interp_spline(&bg->G2, z2_const)*coffe_interp_spline(&bg->D1, z1_const)
               *coffe_interp_spline(&bg->conformal_Hz, z2)*(coffe_interp_spline(&bg->f, z2) - 1)
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)
               *ren1
               +
                chi1*(3 - coffe_interp_spline(&par->evolution_bias2, z2_const))*coffe_interp_spline(&bg->f, z2_const)
               *pow(coffe_interp_spline(&bg->conformal_Hz, z2_const), 2)*coffe_interp_spline(&bg->G1, z1_const)*coffe_interp_spline(&bg->D1, z2_const)
               *coffe_interp_spline(&bg->conformal_Hz, z1)*(coffe_interp_spline(&bg->f, z1) - 1)
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)
               *ren2
            );
    }
    /* g1-g4 + g4-g1 term */
    if (
        par->correlation_contrib.g1 &&
        par->correlation_contrib.g4
    ){
        result +=
            9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/2.
           *(
                (1 + coffe_interp_spline(&bg->G1, z1_const))*(2 - 5*s2)
               *coffe_interp_spline(&bg->D1, z1_const)/coffe_interp_spline(&bg->a, z1_const)
                /* integrand */
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)*ren1
                +
                (1 + coffe_interp_spline(&bg->G2, z2_const))*(2 - 5*s1)
               *coffe_interp_spline(&bg->D1, z2_const)/coffe_interp_spline(&bg->a, z2_const)
                /* integrand */
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)*ren2
            );
    }
    /* g1-g5 + g5-g1 term */
    if (
        par->correlation_contrib.g1 &&
        par->correlation_contrib.g5
    ){
        result +=
            9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/2.
           *(
                chi2*(1 + coffe_interp_spline(&bg->G1, z1_const))*coffe_interp_spline(&bg->G2, z2_const)
               *coffe_interp_spline(&bg->D1, z1_const)/coffe_interp_spline(&bg->a, z1_const)
                /* integrand */
               *coffe_interp_spline(&bg->conformal_Hz, lambda2)
               *(coffe_interp_spline(&bg->f, lambda2) - 1)
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)*ren1
                +
                chi1*(1 + coffe_interp_spline(&bg->G2, z2_const))*coffe_interp_spline(&bg->G1, z1_const)
               *coffe_interp_spline(&bg->D1, z2_const)/coffe_interp_spline(&bg->a, z2_const)
                /* integrand */
                *coffe_interp_spline(&bg->conformal_Hz, lambda1)
               *(coffe_interp_spline(&bg->f, lambda1) - 1)
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)*ren2
            );
    }
    /* g2-g4 + g4-g2 term */
    if (
        par->correlation_contrib.g2 &&
        par->correlation_contrib.g4
    ){
        result +=
            9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/2.
           *(
                (5*s1 - 2)*(2 - 5*s2)
               *coffe_interp_spline(&bg->D1, z1_const)/coffe_interp_spline(&bg->a, z1_const)
                /* integrand */
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)*ren1
                +
                (5*s2 - 2)*(2 - 5*s1)
               *coffe_interp_spline(&bg->D1, z2_const)/coffe_interp_spline(&bg->a, z2_const)
                /* integrand */
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)*ren2
            );
    }
    /* g2-g5 + g5-g2 term */
    if (
        par->correlation_contrib.g2 &&
        par->correlation_contrib.g5
    ){
        result +=
            9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/2.
           *(
                chi2*(5*s1 - 2)*coffe_interp_spline(&bg->G2, z2_const)
               *coffe_interp_spline(&bg->D1, z1_const)/coffe_interp_spline(&bg->a, z1_const)
                /* integrand */
               *coffe_interp_spline(&bg->conformal_Hz, lambda2)
               *(coffe_interp_spline(&bg->f, lambda2) - 1)
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)*ren1
                +
                chi1*(5*s2 - 2)*coffe_interp_spline(&bg->G1, z1_const)
               *coffe_interp_spline(&bg->D1, z2_const)/coffe_interp_spline(&bg->a, z2_const)
                /* integrand */
                *coffe_interp_spline(&bg->conformal_Hz, lambda1)
               *(coffe_interp_spline(&bg->f, lambda1) - 1)
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)*ren2
            );
    }
    /* g3-g4 + g4-g3 term */
    if (
        par->correlation_contrib.g3 &&
        par->correlation_contrib.g4
    ){
        result +=
            9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/2.
           *(
                (coffe_interp_spline(&bg->f, z1_const) - 1)*(2 - 5*s2)
               *coffe_interp_spline(&bg->D1, z1_const)/coffe_interp_spline(&bg->a, z1_const)
                /* integrand */
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)*ren1
                +
                (coffe_interp_spline(&bg->f, z2_const) - 1)*(2 - 5*s1)
               *coffe_interp_spline(&bg->D1, z2_const)/coffe_interp_spline(&bg->a, z2_const)
                /* integrand */
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)*ren2
            );
    }
    /* g3-g5 + g5-g3 term */
    if (
        par->correlation_contrib.g3 &&
        par->correlation_contrib.g5
    ){
        result +=
            9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/2.
           *(
                chi2*(coffe_interp_spline(&bg->f, z1_const) - 1)*coffe_interp_spline(&bg->G2, z2_const)
               *coffe_interp_spline(&bg->D1, z1_const)/coffe_interp_spline(&bg->a, z1_const)
                /* integrand */
               *coffe_interp_spline(&bg->conformal_Hz, lambda2)
               *(coffe_interp_spline(&bg->f, lambda2) - 1)
               *coffe_interp_spline(&bg->D1, z2)/coffe_interp_spline(&bg->a, z2)*ren1
                +
                chi1*(coffe_interp_spline(&bg->f, z2_const) - 1)*coffe_interp_spline(&bg->G1, z1_const)
               *coffe_interp_spline(&bg->D1, z2_const)/coffe_interp_spline(&bg->a, z2_const)
                /* integrand */
                *coffe_interp_spline(&bg->conformal_Hz, lambda1)
               *(coffe_interp_spline(&bg->f, lambda1) - 1)
               *coffe_interp_spline(&bg->D1, z1)/coffe_interp_spline(&bg->a, z1)*ren2
            );
    }

    if (gsl_finite(result)){
        return result;
    }
    else{
        fprintf(stderr,
            "ERROR: in function %s, values:\n"
            "mu = %e\n"
            "z_mean = %e\n"
            "chi_mean = %e\n"
            "sep = %e\n"
            "z1 = %e\n"
            "z2 = %e\n"
            "chi1 = %e\n"
            "chi2 = %e\n"
            "costheta = %e\n"
            "r21 = %e\n"
            "r22 = %e\n"
            "z1_const = %e\n"
            "z2_const = %e\n"
            "s1, s2 = %e, %e\n"
            "b1, b2 = %e, %e\n"
            "x = %e\n",
            __func__,
            mu, z_mean, chi_mean, sep,
            z1, z2, chi1, chi2, costheta,
            r21, r22, z1_const, z2_const,
            s1, s2, b1, b2, x
        );
        exit(EXIT_FAILURE);
    }
}


double functions_double_integrated(
    const coffe_parameters_t *par,
    const coffe_background_t *bg,
    const coffe_integral_array_t integral[],
    const double z_mean,
    const double mu,
    const double sep,
    const double x1,
    const double x2
)
{
    double result = 0;

    const double chi_mean = coffe_interp_spline(&bg->comoving_distance, z_mean);
    const double chi1 = chi_mean - sep*mu/2.;
    const double chi2 = chi_mean + sep*mu/2.;
    const double costheta =
        (2*chi_mean*chi_mean - sep*sep + mu*mu*sep*sep/2.)
       /(2*chi_mean*chi_mean - mu*mu*sep*sep/2.);
    const double lambda1 = chi1*x1, lambda2 = chi2*x2;
    double r2 = lambda1*lambda1 + lambda2*lambda2 - 2*lambda1*lambda2*costheta;
    if (r2 < 0) r2 = 0;

    const double z1_const = coffe_interp_spline(&bg->z_as_chi, chi1);
    const double z2_const = coffe_interp_spline(&bg->z_as_chi, chi2);
    const double z1 = coffe_interp_spline(&bg->z_as_chi, lambda1);
    const double z2 = coffe_interp_spline(&bg->z_as_chi, lambda2);

    const double s1 = coffe_interp_spline(&par->magnification_bias1, z1_const);
    const double s2 = coffe_interp_spline(&par->magnification_bias2, z2_const);
    const double sz_mean1 = coffe_interp_spline(&par->magnification_bias1, z_mean);
    const double sz_mean2 = coffe_interp_spline(&par->magnification_bias2, z_mean);

    double ren = 0;
    if (par->divergent){
        if (r2 <= pow(0.000001*COFFE_H0, 2)){
            ren = coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->renormalization_zero_separation, lambda1);
        }
        else{
            ren = coffe_interp_spline(&coffe_find_integral(integral, 4, 0, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r2))
                    /* renormalization term */
                    -
                    coffe_interp_spline2d(
                        &coffe_find_integral(
                            integral,
                            4,
                            0,
                            COFFE_INTEGER,
                            COFFE_INTEGER
                        )->renormalization,
                        lambda1,
                        lambda2
                    );
        }
    }

    /* len-len term */
    if (
        par->correlation_contrib.len &&
        !par->only_cross_correlations
    ){
        if (par->flatsky_nonlocal && par->output_type != MULTIPOLES){
            /* len-len modified by flatsky */
            result +=
            /* constant in front */
            9 * pow(par->Omega0_cdm + par->Omega0_baryon, 2)
          * (2 - 5 * sz_mean1)
          * (2 - 5 * sz_mean2)
          * pow(chi_mean, 3)
          / 8. / M_PI
            /* integrand */
          * coffe_interp_spline(
                &coffe_find_integral(
                    integral,
                    1,
                    -1,
                    COFFE_HALF_INTEGER,
                    COFFE_HALF_INTEGER
                )->result,
                x1 * sep * sqrt(1 - mu * mu)
            )
            /* D1(z)^2 */
          * pow(
                coffe_interp_spline(
                    &bg->D1,
                    coffe_interp_spline(
                        &bg->z_as_chi,
                        x1 * chi_mean
                    )
                ),
                2
            )
            /* (1 + z)^2 */
          * pow(
                1 + coffe_interp_spline(
                    &bg->z_as_chi,
                    x1 * chi_mean
                ),
                2
            )
          * pow(x1 * (1 - x1), 2);
        }
        else if (!par->flatsky_nonlocal){
            if (r2 > 1e-20){
#ifdef HAVE_CLASS
                if (
                    par->pk_type != COFFE_PK_LINEAR &&
                    par->midpoint_approximation
                ){
                    result +=
                    /* constant in front */
                    9.
                   *(par->Omega0_cdm + par->Omega0_baryon)
                   *(par->Omega0_cdm + par->Omega0_baryon)
                   *(2 - 5 * s1)
                   *(2 - 5 * s2)
                   /4.
                   *chi1
                   *chi2
                    /* integrand */
                   /coffe_interp_spline(&bg->a, z1)
                   /coffe_interp_spline(&bg->a, z2)
                   *(1 - x1)
                   *(1 - x2)
                   *(
                        2 * (costheta * costheta - 1)
                       *lambda1
                       *lambda2
                       *coffe_interp_spline2d(
                            &coffe_find_integral(
                                integral,
                                0,
                                0,
                                COFFE_INTEGER,
                                COFFE_INTEGER
                            )->result2d,
                            (z1 + z2) / 2.,
                            sqrt(r2)
                        )
                       /5.
                       +
                        4 * costheta
                       *coffe_interp_spline2d(
                            &coffe_find_integral(
                                integral,
                                2,
                                0,
                                COFFE_INTEGER,
                                COFFE_INTEGER
                            )->result2d,
                            (z1 + z2) / 2.,
                            sqrt(r2)
                        )
                       /3.
                       +
                        4 * costheta
                       *(r2 + 6 * costheta * lambda1 * lambda2)
                       *coffe_interp_spline2d(
                            &coffe_find_integral(
                                integral,
                                1,
                                1,
                                COFFE_INTEGER,
                                COFFE_INTEGER
                            )->result2d,
                            (z1 + z2) / 2.,
                            sqrt(r2)
                        )
                       /15.
                        +
                        2 * (costheta * costheta - 1)
                       *lambda1
                       *lambda2
                       *(2 * r2 + 3 * costheta * lambda1 * lambda2)
                       *coffe_interp_spline2d(
                            &coffe_find_integral(
                                integral,
                                0,
                                2,
                                COFFE_INTEGER,
                                COFFE_INTEGER
                            )->result2d,
                            (z1 + z2) / 2.,
                            sqrt(r2)
                        )
                       /7.
                       /r2
                        +
                        2 * costheta
                       *(
                            2 * r2 * r2
                            +
                            12 * costheta * r2 * lambda1 * lambda2
                            +
                            15 * (costheta * costheta - 1)
                           *lambda1 * lambda1 * lambda2 * lambda2
                        )
                       *coffe_interp_spline2d(
                            &coffe_find_integral(
                                integral,
                                1,
                                3,
                                COFFE_INTEGER,
                                COFFE_INTEGER
                            )->result2d,
                            (z1 + z2) / 2.,
                            sqrt(r2)
                        )
                       /15.
                       /r2
                        +
                        (costheta * costheta - 1)
                       *lambda1 * lambda2
                       *(
                            6 * r2 * r2
                            +
                            30 * costheta * r2 * lambda1 * lambda2
                            +
                            35 * (costheta * costheta - 1)
                           *lambda1 * lambda1 * lambda2 * lambda2
                        )
                       *coffe_interp_spline2d(
                            &coffe_find_integral(
                                integral,
                                0,
                                4,
                                COFFE_INTEGER,
                                COFFE_INTEGER
                            )->result2d,
                            (z1 + z2) / 2.,
                            sqrt(r2)
                        )
                       /35.
                       /r2
                       /r2
                    );
                }
                else
#endif
                    result +=
                    /* constant in front */
                    9.
                   *(par->Omega0_cdm + par->Omega0_baryon)
                   *(par->Omega0_cdm + par->Omega0_baryon)
                    * (2 - 5 * s1)
                    * (2 - 5 * s2)
                   /4.
                   *chi1
                   *chi2
                    /* integrand */
                   *coffe_interp_spline(&bg->D1, z1)
                   *coffe_interp_spline(&bg->D1, z2)
                   /coffe_interp_spline(&bg->a, z1)
                   /coffe_interp_spline(&bg->a, z2)
                   *(1 - x1)
                   *(1 - x2)
                   *(
                        2 * (costheta * costheta - 1) * lambda1 * lambda2
                       *coffe_interp_spline(
                            &coffe_find_integral(
                                integral,
                                0,
                                0,
                                COFFE_INTEGER,
                                COFFE_INTEGER
                            )->result,
                            sqrt(r2)
                        )/5.
                       +
                        4 * costheta
                       *coffe_interp_spline(
                            &coffe_find_integral(
                                integral,
                                2,
                                0,
                                COFFE_INTEGER,
                                COFFE_INTEGER
                            )->result,
                            sqrt(r2)
                        )/3.
                       +
                        4 * costheta * (r2 + 6 * costheta * lambda1 * lambda2)
                       *coffe_interp_spline(
                            &coffe_find_integral(
                                integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER
                            )->result, sqrt(r2))/15.
                       +
                        2 * (costheta * costheta - 1) * lambda1 * lambda2
                        * (2 * r2 + 3 * costheta * lambda1 * lambda2)
                       *coffe_interp_spline(
                            &coffe_find_integral(
                                integral, 0, 2, COFFE_INTEGER, COFFE_INTEGER
                            )->result, sqrt(r2))/7./r2
                       +
                        2*costheta
                        * (2 * r2 * r2 + 12 * costheta * r2 * lambda1 * lambda2 + 15 * (costheta * costheta - 1) * lambda1 * lambda1 * lambda2 * lambda2)
                       *coffe_interp_spline(
                            &coffe_find_integral(
                                integral, 1, 3, COFFE_INTEGER, COFFE_INTEGER
                            )->result, sqrt(r2))/15./r2
                       +
                        (costheta * costheta - 1) * lambda1 * lambda2
                        * (6 * r2 * r2 + 30 * costheta * r2 * lambda1 * lambda2 + 35 * (costheta * costheta - 1) * lambda1 * lambda1 * lambda2 * lambda2)
                       *coffe_interp_spline(
                            &coffe_find_integral(
                                integral,
                                0,
                                4,
                                COFFE_INTEGER,
                                COFFE_INTEGER
                            )->result,
                            sqrt(r2)
                        )/ 35. / r2 / r2
                    );
            }
            else{
#ifdef HAVE_CLASS
                if (
                    par->pk_type != COFFE_PK_LINEAR &&
                    par->midpoint_approximation
                ){
                    result +=
                    /* constant in front */
                    9./4
                   *pow((par->Omega0_cdm + par->Omega0_baryon), 2)
                   *(2 - 5 * s1)
                   *(2 - 5 * s2)
                   *chi1
                   *chi2
                    /* integrand */
                   /coffe_interp_spline(&bg->a, z1)
                   /coffe_interp_spline(&bg->a, z2)
                   *(1 - x1)
                   *(1 - x2)
                   *(
                       4.
                       +
                        24. * lambda1 * lambda2
                       *coffe_interp_spline2d(
                            &coffe_find_integral(
                                integral,
                                2,
                                0,
                                COFFE_INTEGER,
                                COFFE_INTEGER
                            )->result2d,
                            (z1 + z2) / 2.,
                            sqrt(r2)
                        )
                       /15.
                    );
                }
                else
#endif
                    result +=
                    /* constant in front */
                    9. / 4*pow((par->Omega0_cdm + par->Omega0_baryon), 2)
                    * (2 - 5 * s1)
                    * (2 - 5 * s2)
                   *chi1
                   *chi2
                    /* integrand */
                   *coffe_interp_spline(&bg->D1, z1)
                   *coffe_interp_spline(&bg->D1, z2)
                   /coffe_interp_spline(&bg->a, z1)
                   /coffe_interp_spline(&bg->a, z2)
                   *(1 - x1)
                   *(1 - x2)
                   *(
                       4.
                       +
                        24. * lambda1 * lambda2
                       *coffe_interp_spline(&coffe_find_integral(integral, 1, 1, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)/15.
                    );
            }
        }
    }
    /* g4-g4 term */
    if (
        par->correlation_contrib.g4 &&
        !par->only_cross_correlations
    ){
        result +=
        /* constant in front */
        9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)*(2 - 5*s1)*(2 - 5*s2)
       *
            /* integrand */
            coffe_interp_spline(&bg->D1, z1)
           *coffe_interp_spline(&bg->D1, z2)
           /coffe_interp_spline(&bg->a, z1)
           /coffe_interp_spline(&bg->a, z2)
           *ren;
    }
    /* g5-g5 term */
    if (
        par->correlation_contrib.g5 &&
        !par->only_cross_correlations
    ){
        result +=
        /* constant in front */
        9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)
       *coffe_interp_spline(&bg->G1, z1_const)
       *coffe_interp_spline(&bg->G2, z2_const)
       *chi1*chi2
       *
        /* integrand */
            coffe_interp_spline(&bg->D1, z1)
           *coffe_interp_spline(&bg->D1, z2)
           /coffe_interp_spline(&bg->a, z1)
           /coffe_interp_spline(&bg->a, z2)
           *coffe_interp_spline(&bg->conformal_Hz, z1)
           *coffe_interp_spline(&bg->conformal_Hz, z2)
           *(coffe_interp_spline(&bg->f, z1) - 1)
           *(coffe_interp_spline(&bg->f, z2) - 1)
           *ren;
    }
    /* g4-len + len-g4 term */
    if (
        par->correlation_contrib.g4 &&
        par->correlation_contrib.len
    ){
        if (r2 != 0){
            result +=
                /* constant in front */
                9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    (2 - 5*s1)*(2 - 5*s2)
                   *(1 - x2)/x2*coffe_interp_spline(&bg->D1, z1)*coffe_interp_spline(&bg->D1, z2)
                   /coffe_interp_spline(&bg->a, z1)/coffe_interp_spline(&bg->a, z2)
                   *(
                        2*lambda1*lambda2*costheta*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r2))
                       -lambda1*lambda1*lambda2*lambda2*(1 - costheta*costheta)*coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r2))
                    )
                    +
                    (2 - 5*s1)*(2 - 5*s2)
                   *(1 - x1)/x1*coffe_interp_spline(&bg->D1, z1)*coffe_interp_spline(&bg->D1, z2)
                   /coffe_interp_spline(&bg->a, z1)/coffe_interp_spline(&bg->a, z2)
                   *(
                        2*lambda1*lambda2*costheta*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r2))
                       -lambda1*lambda1*lambda2*lambda2*(1 - costheta*costheta)*coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r2))
                    )
                );
        }
        else{
            result +=
                9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    (2 - 5*s1)*(2 - 5*s2)
                   *(1 - x2)/x2*coffe_interp_spline(&bg->D1, z1)*coffe_interp_spline(&bg->D1, z2)
                   /coffe_interp_spline(&bg->a, z1)/coffe_interp_spline(&bg->a, z2)
                   *2*lambda1*lambda2*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)
                    +
                    (2 - 5*s1)*(2 - 5*s2)
                   *(1 - x1)/x1*coffe_interp_spline(&bg->D1, z2)*coffe_interp_spline(&bg->D1, z1)
                   /coffe_interp_spline(&bg->a, z2)/coffe_interp_spline(&bg->a, z1)
                   *2*lambda1*lambda2*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)
                );
        }
    }
    /* g5-len + len-g5 term */
    if (
        par->correlation_contrib.g5 &&
        par->correlation_contrib.len
    ){
        if (r2 != 0){
            result +=
                /* constant in front */
                9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    (2 - 5*s2)*coffe_interp_spline(&bg->G1, z1_const)*chi1
                   *coffe_interp_spline(&bg->conformal_Hz, z1)*(coffe_interp_spline(&bg->conformal_Hz, z1) - 1)
                   *(1 - x2)/x2*coffe_interp_spline(&bg->D1, z1)*coffe_interp_spline(&bg->D1, z2)
                   /coffe_interp_spline(&bg->a, z1)/coffe_interp_spline(&bg->a, z2)
                   *(
                        2*lambda1*lambda2*costheta*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r2))
                       -lambda1*lambda1*lambda2*lambda2*(1 - costheta*costheta)*coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r2))
                    )
                    +
                    (2 - 5*s1)*coffe_interp_spline(&bg->G2, z2_const)*chi2
                   *coffe_interp_spline(&bg->conformal_Hz, z2)*(coffe_interp_spline(&bg->conformal_Hz, z2) - 1)
                   *(1 - x1)/x1*coffe_interp_spline(&bg->D1, z1)*coffe_interp_spline(&bg->D1, z2)
                   /coffe_interp_spline(&bg->a, z1)/coffe_interp_spline(&bg->a, z2)
                   *(
                        2*lambda1*lambda2*costheta*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r2))
                       -lambda1*lambda1*lambda2*lambda2*(1 - costheta*costheta)*coffe_interp_spline(&coffe_find_integral(integral, 2, 2, COFFE_INTEGER, COFFE_INTEGER)->result, sqrt(r2))
                    )
                );
        }
        else{
            result +=
                9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    (2 - 5*s2)*coffe_interp_spline(&bg->G1, z1_const)*chi1
                   *coffe_interp_spline(&bg->conformal_Hz, z1)*(coffe_interp_spline(&bg->conformal_Hz, z1) - 1)
                   *(1 - x2)/x2*coffe_interp_spline(&bg->D1, z1)*coffe_interp_spline(&bg->D1, z2)
                   /coffe_interp_spline(&bg->a, z1)/coffe_interp_spline(&bg->a, z2)
                   *2*lambda1*lambda2*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)
                    +
                    (2 - 5*s1)*coffe_interp_spline(&bg->G2, z2_const)*chi2
                   *coffe_interp_spline(&bg->conformal_Hz, z2)*(coffe_interp_spline(&bg->conformal_Hz, z2) - 1)
                   *(1 - x1)/x1*coffe_interp_spline(&bg->D1, z1)*coffe_interp_spline(&bg->D1, z2)
                   /coffe_interp_spline(&bg->a, z1)/coffe_interp_spline(&bg->a, z2)
                   *2*lambda1*lambda2*coffe_interp_spline(&coffe_find_integral(integral, 3, 1, COFFE_INTEGER, COFFE_INTEGER)->result, 0.0)
                );
        }
    }
    /* g4-g5 + g5-g4 term */
    if (
        par->correlation_contrib.g4 &&
        par->correlation_contrib.g5
    ){
        result +=
            /* constant in front */
            9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)
           *(
                coffe_interp_spline(&bg->G2, z2_const)*(2 - 5*s1)*chi2
               *coffe_interp_spline(&bg->conformal_Hz, z2)*(coffe_interp_spline(&bg->f, z2) - 1)
               *coffe_interp_spline(&bg->D1, z1)*coffe_interp_spline(&bg->D1, z2)
               /coffe_interp_spline(&bg->a, z1)/coffe_interp_spline(&bg->a, z2)
               *ren
               +
                coffe_interp_spline(&bg->G1, z1_const)*(2 - 5*s2)*chi1
               *coffe_interp_spline(&bg->conformal_Hz, z1)*(coffe_interp_spline(&bg->f, z1) - 1)
               *coffe_interp_spline(&bg->D1, z1)*coffe_interp_spline(&bg->D1, z2)
               /coffe_interp_spline(&bg->a, z1)/coffe_interp_spline(&bg->a, z2)
               *ren
            );
    }
    if (gsl_finite(result)){
    return
        result;
    }
    else{
        fprintf(stderr,
            "ERROR: in function %s, values:\n"
            "x1 = %e, x2 = %e\n"
            "r2 = %e\n"
            "mu = %e\n"
            "z_mean = %e\n"
            "chi_mean = %e\n"
            "sep = %e\n"
            "z1 = %e\n"
            "z2 = %e\n"
            "chi1 = %e\n"
            "chi2 = %e\n"
            "lambda1 = %e\n"
            "lambda2 = %e\n"
            "ren(sqrt(r2)) = %e\n",
            __func__, x1, x2, r2,
            mu, z_mean, chi_mean, sep,
            z1, z2, chi1, chi2, lambda1, lambda2, ren
        );
        exit(EXIT_FAILURE);
    }
}


/* lensing-lensing multipoles are special */
double functions_flatsky_lensing_lensing_multipoles(
    const coffe_parameters_t *par,
    const coffe_background_t *bg,
    const coffe_integral_array_t *integral,
    const double z_mean,
    const double sep,
    const int l,
    const double x
)
{
    const double chi_mean = coffe_interp_spline(
        &bg->comoving_distance,
        z_mean
    );
    const double lambda = chi_mean * x;

    const double sz_mean1 = coffe_interp_spline(
        &par->magnification_bias1,
        z_mean
    );
    const double sz_mean2 = coffe_interp_spline(
        &par->magnification_bias2,
        z_mean
    );

    double result = 0;
#ifdef HAVE_CLASS
    if (
        par->pk_type != COFFE_PK_LINEAR &&
        par->midpoint_approximation
    ){
        result =
            (2 * l + 1) * gsl_sf_fact(l) / pow(2, l) / pow(gsl_sf_fact(l / 2), 2)
          * 9 * pow(par->Omega0_cdm + par->Omega0_baryon, 2)
          * (2 - 5 * sz_mean1)
          * (2 - 5 * sz_mean2)
          * pow(chi_mean, 3)
          / 8. / M_PI
            /* integrand */
          * coffe_interp_spline2d(
                &coffe_find_integral(
                    integral,
                    1,
                    l,
                    COFFE_INTEGER,
                    COFFE_INTEGER
                )->result2d,
                /* must be interpolated at z(\lambda) */
                coffe_interp_spline(
                    &bg->z_as_chi,
                    lambda
                ),
                x * sep
            )
            /* (1 + z)^2 */
          * pow(
                1 + coffe_interp_spline(
                    &bg->z_as_chi,
                    lambda
                ),
                2
            )
          * pow(x * (1 - x), 2);
    }
    else{
#endif
        result =
            (2 * l + 1) * gsl_sf_fact(l) / pow(2, l) / pow(gsl_sf_fact(l / 2), 2)
          * 9 * pow(par->Omega0_cdm + par->Omega0_baryon, 2)
          * (2 - 5 * sz_mean1)
          * (2 - 5 * sz_mean2)
          * pow(chi_mean, 3)
          / 8. / M_PI
            /* integrand */
          * coffe_interp_spline(
                &coffe_find_integral(
                    integral,
                    1,
                    l,
                    COFFE_INTEGER,
                    COFFE_INTEGER
                )->result,
                x * sep
            )
            /* D1(z)^2 (only in linear theory) */
          * pow(
                coffe_interp_spline(
                    &bg->D1,
                    coffe_interp_spline(
                        &bg->z_as_chi,
                        lambda
                    )
                ),
                2
            )
            /* (1 + z)^2 */
          * pow(
                1 + coffe_interp_spline(
                    &bg->z_as_chi,
                    lambda
                ),
                2
            )
          * pow(x * (1 - x), 2);
#ifdef HAVE_CLASS
    }
#endif

    if (l > 0)
        return result * x * sep;
    else
        return result;
}


/* density-lensing multipoles are special */
double functions_flatsky_density_lensing_multipoles(
    const coffe_parameters_t *par,
    const coffe_background_t *bg,
    const coffe_integral_array_t *integral,
    const double z_mean,
    const double sep,
    const int l
)
{
    const double sz_mean1 = coffe_interp_spline(
        &par->magnification_bias1,
        z_mean
    );
    const double sz_mean2 = coffe_interp_spline(
        &par->magnification_bias2,
        z_mean
    );

    const double bz_mean1 = coffe_interp_spline(
        &par->galaxy_bias1,
        z_mean
    );
    const double bz_mean2 = coffe_interp_spline(
        &par->galaxy_bias2,
        z_mean
    );

    double result = 0;

#ifdef HAVE_CLASS
    if (
        par->pk_type != COFFE_PK_LINEAR &&
        par->midpoint_approximation
    ){
        /* part with even multipoles */
        if (l % 2 == 0){
            for (size_t k = 0; k <= (size_t)l / 2; ++k){
                result +=
                    pow(-1, k)
                   /pow(2, k)
                   *gsl_sf_choose(l, k)
                   *gsl_sf_choose(2 * l - 2 * k, l)
                   *gsl_sf_fact(l / 2 - k)
                   *coffe_interp_spline2d(
                        &coffe_find_integral(
                            integral,
                            l - 2 * k + 3,
                            l - 2 * k + 1,
                            COFFE_HALF_INTEGER,
                            COFFE_HALF_INTEGER
                        )->result2d,
                        z_mean,
                        sep
                    );
            }
            result *=
               -3 * (par->Omega0_cdm + par->Omega0_baryon) / M_PI / 8
               *(2 * l + 1)
               *(1 + z_mean)
               *sep
               *(
                    (2 - 5 * sz_mean1) * bz_mean2
                    +
                    (2 - 5 * sz_mean2) * bz_mean1
                )
               *pow(2, 1 - l / 2);
        }
        /* part with odd multipoles (TODO) */
        else{
            result += 0;
        }    }
    else{
#endif
        /* part with even multipoles */
        if (l % 2 == 0){
            for (size_t k = 0; k <= (size_t)l / 2; ++k){
                result +=
                    pow(-1, k)
                   /pow(2, k)
                   *gsl_sf_choose(l, k)
                   *gsl_sf_choose(2 * l - 2 * k, l)
                   *gsl_sf_fact(l / 2 - k)
                   *coffe_interp_spline(
                        &coffe_find_integral(
                            integral,
                            l - 2 * k + 3,
                            l - 2 * k + 1,
                            COFFE_HALF_INTEGER,
                            COFFE_HALF_INTEGER
                        )->result,
                        sep
                    );
            }
            result *=
               -3 * (par->Omega0_cdm + par->Omega0_baryon) / M_PI / 8
               *(2 * l + 1)
               *coffe_interp_spline(&bg->D1, z_mean)
               *coffe_interp_spline(&bg->D1, z_mean)
               *(1 + z_mean)
               *sep
               *(
                    (2 - 5 * sz_mean1) * bz_mean2
                    +
                    (2 - 5 * sz_mean2) * bz_mean1
                )
               *pow(2, 1 - l / 2);

        }
        /* part with odd multipoles (TODO) */
        else{
            result += 0;
        }
#ifdef HAVE_CLASS
    }
#endif

    return result / sqrt(2 * M_PI);
}
