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

/**
    all the nonintegrated terms in one place
**/

double functions_nonintegrated(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t integral[],
    double z_mean,
    double mu,
    double sep
)
{
    double chi_mean = interp_spline(&bg->comoving_distance, z_mean);
    double chi1 = chi_mean - sep*mu/2.;
    double chi2 = chi_mean + sep*mu/2.;
    double costheta =
        (2.*chi_mean*chi_mean - sep*sep + mu*mu*sep*sep/2.)
       /(2.*chi_mean*chi_mean - mu*mu*sep*sep/2.);

    double result = 0;
    double z1 = interp_spline(&bg->z_as_chi, chi1);
    double z2 = interp_spline(&bg->z_as_chi, chi2);
    double f1 = interp_spline(&bg->f, z1);
    double f2 = interp_spline(&bg->f, z2);
    double fmean = interp_spline(&bg->f, z_mean);
    double curlyH1 = interp_spline(&bg->conformal_Hz, z1); // dimensionless
    double curlyH2 = interp_spline(&bg->conformal_Hz, z2); // dimensionless
    double b1 = interp_spline(&par->matter_bias1, z1);
    double b2 = interp_spline(&par->matter_bias2, z2);
    double bz_mean1 = interp_spline(&par->matter_bias1, z_mean);
    double bz_mean2 = interp_spline(&par->matter_bias2, z_mean);
    double G1 = interp_spline(&bg->G1, z1);
    double G2 = interp_spline(&bg->G2, z2);
    double s1 = interp_spline(&par->magnification_bias1, z1);
    double s2 = interp_spline(&par->magnification_bias2, z2);
    double fevo1 = interp_spline(&par->evolution_bias1, z1);
    double fevo2 = interp_spline(&par->evolution_bias2, z2);
    double a1 = interp_spline(&bg->a, z1);
    double a2 = interp_spline(&bg->a, z2);
    /* den-den term */
    if (par->correlation_contrib.den){
        /* den-den modified by flatsky */
        if (par->flatsky){
            result += bz_mean1*bz_mean2
               *interp_spline(&integral[0].result, sep);
        }
        else{
        result += b1*b2
           *interp_spline(&integral[0].result, sep);
        }
    }
    /* rsd-rsd term */
    if (par->correlation_contrib.rsd){
        /* rsd-rsd modified by flatsky */
        if (par->flatsky){
            result +=
                fmean*fmean*interp_spline(&integral[0].result, sep)/5.
                -
                4*fmean*fmean*interp_spline(&integral[1].result, sep)/7.
               *gsl_sf_legendre_P2(mu)
                +
                8*fmean*fmean*interp_spline(&integral[2].result, sep)/35.
               *gsl_sf_legendre_Pl(4, mu);
        }
        else{
        result +=
            f1*f2*(1 + 2*pow(costheta, 2))/15
           *interp_spline(&integral[0].result, sep)
            -
            f1*f2/21.*(
                (1 + 11.*pow(costheta, 2)) + 18*costheta*(pow(costheta, 2) - 1)*chi1*chi2/sep/sep
            )
           *interp_spline(&integral[1].result, sep)
            +
            f1*f2*(
                4*(3*pow(costheta, 2) - 1)*(pow(chi1, 4) + pow(chi2, 4))/35./pow(sep, 4)
                +
                chi1*chi2*(3 + pow(costheta, 2))*(
                    3*(3 + pow(costheta, 2))*chi1*chi2 - 8*(pow(chi1, 2) + pow(chi2, 2))*costheta
                )/35./pow(sep, 4)
            )
           *interp_spline(&integral[2].result, sep);
        }
    }
    /* d1-d1 term */
    if (par->correlation_contrib.d1){
        result +=
            (
                curlyH1*curlyH2*f1*f2*G1*G2
               *costheta/3.
               *interp_spline(&integral[5].result, sep)
                +
                curlyH1*curlyH2*f1*f2*G1*G2
               *(
                    (chi2 - chi1*costheta)*(chi1 - chi2*costheta)
                    +
                    pow(sep, 2)*costheta/3.
                )
               *interp_spline(&integral[6].result, sep)
            );
    }
    /* d2-d2 term */
    if (par->correlation_contrib.d2){
        result +=
            (3 - fevo1)*(3 - fevo2)*pow(curlyH1, 2)*pow(curlyH2, 2)*f1*f2
           *(
                interp_spline(&integral[8].result, sep)
                /* renormalization term */
               -gsl_spline2d_eval(
                    integral[8].renormalization.spline,
                    chi1, chi2,
                    integral[8].renormalization.xaccel,
                    integral[8].renormalization.yaccel
                )
            );
    }
    /* g1-g1 term */
    if (par->correlation_contrib.g1){
        result += 9*pow((par->Omega0_cdm + par->Omega0_baryon), 2)
           *(1 + G1)*(1 + G2)/4/a1/a2
           *(
                interp_spline(&integral[8].result, sep)
                /* renormalization term */
               -gsl_spline2d_eval(
                    integral[8].renormalization.spline,
                    chi1, chi2,
                    integral[8].renormalization.xaccel,
                    integral[8].renormalization.yaccel
                )
            );
    }
    /* g2-g2 term */
    if (par->correlation_contrib.g2){
        result += 9*pow((par->Omega0_cdm + par->Omega0_baryon), 2)
           *(5*s1 - 2)*(5*s2 - 2)/4/a1/a2
           *(
                interp_spline(&integral[8].result, sep)
                /* renormalization term */
               -gsl_spline2d_eval(
                    integral[8].renormalization.spline,
                    chi1, chi2,
                    integral[8].renormalization.xaccel,
                    integral[8].renormalization.yaccel
                )
            );
    }
    /* g3-g3 term */
    if (par->correlation_contrib.g3){
        result += 9*pow((par->Omega0_cdm + par->Omega0_baryon), 2)
           *(f1 - 1)*(f2 - 1)/4/a1/a2
           *(
                interp_spline(&integral[8].result, sep)
                /* renormalization term */
               -gsl_spline2d_eval(
                    integral[8].renormalization.spline,
                    chi1, chi2,
                    integral[8].renormalization.xaccel,
                    integral[8].renormalization.yaccel
                )
            );
    }
    /* den-rsd + rsd-den term */
    if (
        par->correlation_contrib.den &&
        par->correlation_contrib.rsd
    ){
        /* den-rsd modified by flatsky */
        if (par->flatsky){
            result +=
                (bz_mean1*fmean/3. + bz_mean2*fmean/3.)
               *interp_spline(&integral[0].result, sep)
               -
                (2*bz_mean1*fmean/3. + 2*bz_mean2*fmean/3.)
               *interp_spline(&integral[1].result, sep)
               *gsl_sf_legendre_P2(mu);
        }
        else{
        result += (b1*f2/3. + b2*f1/3.)
           *interp_spline(&integral[0].result, sep)
           -
            (
                b1*f2*(2./3. - (1. - pow(costheta, 2))*pow(chi1/sep, 2))
                +
                b2*f1*(2./3. - (1. - pow(costheta, 2))*pow(chi2/sep, 2))
            )
           *interp_spline(&integral[1].result, sep);
        }
    }
    /* den-d1 + d1-den term */
    if (
        par->correlation_contrib.den &&
        par->correlation_contrib.d1
    ){
        result += -(
                b1*f2*curlyH2*G2*(chi1*costheta - chi2)
                +
                b2*f1*curlyH1*G1*(chi2*costheta - chi1)
            )
           *interp_spline(&integral[3].result, sep);
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
           *interp_spline(&integral[5].result, sep);
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
           *interp_spline(&integral[5].result, sep);
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
           *interp_spline(&integral[5].result, sep);
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
           *interp_spline(&integral[5].result, sep);
    }
    /* rsd-d1 + d1-rsd term */
    if (
        par->correlation_contrib.rsd &&
        par->correlation_contrib.d1
    ){
        result += (
            (
                f1*f2*curlyH2*G2*((1. + 2*pow(costheta, 2))*chi2 - 3*chi1*costheta)/5.
                +
                f2*f1*curlyH1*G1*((1. + 2*pow(costheta, 2))*chi1 - 3*chi2*costheta)/5.
            )
           *interp_spline(&integral[3].result, sep)
            + (
                f1*f2*curlyH2*G2*(
                    (1. - 3*costheta*costheta)*pow(chi2, 3)
                    +
                    costheta*(5. + pow(costheta, 2))*pow(chi2, 2)*chi1
                    -
                    2*(2. + pow(costheta, 2))*chi2*pow(chi1, 2)
                    +
                    2*pow(chi1, 3)*costheta
                )/5
                +
                f2*f1*curlyH1*G1*(
                    (1. - 3*costheta*costheta)*pow(chi1, 3)
                    +
                    costheta*(5. + pow(costheta, 2))*pow(chi1, 2)*chi2
                    -
                    2*(2. + pow(costheta, 2))*chi1*pow(chi2, 2)
                    +
                    2*pow(chi2, 3)*costheta
                )/5
            )
           *interp_spline(&integral[4].result, sep)/pow(sep, 2)
        );
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
           *interp_spline(&integral[5].result, sep)
            - (
                (3 - fevo2)*f1*f2*pow(curlyH2, 2)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi2, 2))
                +
                (3 - fevo1)*f2*f1*pow(curlyH1, 2)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi1, 2))
            )
           *interp_spline(&integral[6].result, sep)
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
           *interp_spline(&integral[5].result, sep)
            + (
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a2*f1*(1 + G2)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi2, 2))
                +
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a1*f2*(1 + G1)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi1, 2))
            )
           *interp_spline(&integral[6].result, sep);
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
           *interp_spline(&integral[5].result, sep)
            + (
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a2*f1*(5*s2 - 2)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi2, 2))
                +
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a1*f2*(5*s1 - 2)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi1, 2))
            )
           *interp_spline(&integral[6].result, sep);
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
           *interp_spline(&integral[5].result, sep)
            + (
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a2*f1*(f2 - 1)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi2, 2))
                +
                3*(par->Omega0_cdm + par->Omega0_baryon)/2./a1*f2*(f1 - 1)*(2./3*pow(sep, 2) - (1 - pow(costheta, 2))*pow(chi1, 2))
            )
           *interp_spline(&integral[6].result, sep);

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
           *interp_spline(&integral[7].result, sep);
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
           *interp_spline(&integral[7].result, sep);
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
           *interp_spline(&integral[7].result, sep);
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
           *interp_spline(&integral[7].result, sep);
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
                interp_spline(&integral[8].result, sep)
                /* renormalization term */
               -gsl_spline2d_eval(
                    integral[8].renormalization.spline,
                    chi1, chi2,
                    integral[8].renormalization.xaccel,
                    integral[8].renormalization.yaccel
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
                interp_spline(&integral[8].result, sep)
                /* renormalization term */
               -gsl_spline2d_eval(
                    integral[8].renormalization.spline,
                    chi1, chi2,
                    integral[8].renormalization.xaccel,
                    integral[8].renormalization.yaccel
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
                interp_spline(&integral[8].result, sep)
                /* renormalization term */
               -gsl_spline2d_eval(
                    integral[8].renormalization.spline,
                    chi1, chi2,
                    integral[8].renormalization.xaccel,
                    integral[8].renormalization.yaccel
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
                interp_spline(&integral[8].result, sep)
                /* renormalization term */
               -gsl_spline2d_eval(
                    integral[8].renormalization.spline,
                    chi1, chi2,
                    integral[8].renormalization.xaccel,
                    integral[8].renormalization.yaccel
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
                interp_spline(&integral[8].result, sep)
                /* renormalization term */
               -gsl_spline2d_eval(
                    integral[8].renormalization.spline,
                    chi1, chi2,
                    integral[8].renormalization.xaccel,
                    integral[8].renormalization.yaccel
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
                interp_spline(&integral[8].result, sep)
            /* renormalization term */
               -gsl_spline2d_eval(
                    integral[8].renormalization.spline,
                    chi1, chi2,
                    integral[8].renormalization.xaccel,
                    integral[8].renormalization.yaccel
                )
            );
    }
    if (gsl_finite(result)){
    return
        result*interp_spline(&bg->D1, z1)*interp_spline(&bg->D1, z2);
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
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t integral[],
    double z_mean,
    double mu,
    double sep,
    double x
)
{
    double result = 0;

    double chi_mean = interp_spline(&bg->comoving_distance, z_mean);
    double chi1 = chi_mean - sep*mu/2.;
    double chi2 = chi_mean + sep*mu/2.;
    double costheta =
        (2*chi_mean*chi_mean - sep*sep + mu*mu*sep*sep/2.)
       /(2*chi_mean*chi_mean - mu*mu*sep*sep/2.);
    double lambda1 = chi1*x, lambda2 = chi2*x;

    double r21 = lambda2*lambda2 + chi1*chi1 - 2*chi1*lambda2*costheta;
    double r22 = lambda1*lambda1 + chi2*chi2 - 2*chi2*lambda1*costheta;
    if (r21 < 0) r21 = 0;
    if (r22 < 0) r22 = 0;
    double z1_const = interp_spline(&bg->z_as_chi, chi1);
    double z2_const = interp_spline(&bg->z_as_chi, chi2);
    double z1 = interp_spline(&bg->z_as_chi, lambda1);
    double z2 = interp_spline(&bg->z_as_chi, lambda2);

    double s1 = interp_spline(&par->magnification_bias1, z1_const);
    double s2 = interp_spline(&par->magnification_bias2, z2_const);
    double sz_mean1 = interp_spline(&par->magnification_bias1, z_mean);
    double sz_mean2 = interp_spline(&par->magnification_bias2, z_mean);
    double b1 = interp_spline(&par->matter_bias1, z1_const);
    double b2 = interp_spline(&par->matter_bias2, z2_const);
    double bz_mean1 = interp_spline(&par->matter_bias1, z_mean);
    double bz_mean2 = interp_spline(&par->matter_bias2, z_mean);

    double ren1 = 0, ren2 = 0;
    if (par->divergent){
        if (r21 == 0.0) ren1 = interp_spline(&integral[8].renormalization0, lambda2);
        else ren1 = interp_spline(&integral[8].result, sqrt(r21))
                    /* renormalization term */
                   -gsl_spline2d_eval(
                        integral[8].renormalization.spline,
                        lambda2, chi1,
                        integral[8].renormalization.xaccel,
                        integral[8].renormalization.yaccel
                );
        if (r22 == 0.0) ren2 = interp_spline(&integral[8].renormalization0, lambda1);
        else ren2 = interp_spline(&integral[8].result, sqrt(r22))
                    /* renormalization term */
                   -gsl_spline2d_eval(
                        integral[8].renormalization.spline,
                        lambda1, chi2,
                        integral[8].renormalization.xaccel,
                        integral[8].renormalization.yaccel
                );
    }

    /* den-len + len-den term */
    if (
        par->correlation_contrib.den &&
        par->correlation_contrib.len
    ){
        /* den-len + len-den modified in flatsky */
        if (par->flatsky){
            result +=
               -3*(par->Omega0_cdm + par->Omega0_baryon)/M_PI/4
               *interp_spline(&bg->D1, z_mean)*(1 + z_mean)*fabs(mu)*sep
               *((2 - 5*sz_mean1)*bz_mean2 + (2 - 5*sz_mean2)*bz_mean1)
               *interp_spline(&integral[9].result, sep*sqrt(1 - mu*mu));
        }
        else{
        if (r21 != 0.0 && r22 != 0.0){
            result +=
               -3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    b1*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)*chi2
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                        2*chi1*costheta*interp_spline(&integral[3].result, sqrt(r21))
                       -chi1*chi1*lambda2*(1 - costheta*costheta)
                       *interp_spline(&integral[1].result, sqrt(r21))
                       /r21
                    )
                    +
                    b2*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)*chi1
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                        2*chi2*costheta*interp_spline(&integral[3].result, sqrt(r22))
                       -chi2*chi2*lambda1*(1 - costheta*costheta)
                       *interp_spline(&integral[1].result, sqrt(r22))
                       /r22
                    )
                );
        }
        else if (r21 == 0.0 && r22 != 0){
            result +=
               -3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    b1*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)*chi2
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                        2*chi1*interp_spline(&integral[3].result, 0.0)
                    )
                    +
                    b2*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)*chi1
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                        2*chi2*costheta*interp_spline(&integral[3].result, sqrt(r22))
                       -chi2*chi2*lambda1*(1 - costheta*costheta)
                       *interp_spline(&integral[1].result, sqrt(r22))
                       /r22
                    )
                );
        }
        else if (r21 != 0 && r22 == 0.0){
            result +=
               -3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    b1*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)*chi2
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                        2*chi1*costheta*interp_spline(&integral[3].result, sqrt(r21))
                       -chi1*chi1*lambda2*(1 - costheta*costheta)
                       *interp_spline(&integral[1].result, sqrt(r21))
                       /r21
                    )
                    +
                    b2*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)*chi1
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                        2*chi2*interp_spline(&integral[3].result, 0.0)
                    )
                );
        }
        else{
            result +=
               -3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    b1*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                        2*chi1*chi2*interp_spline(&integral[3].result, 0.0)
                    )
                    +
                    b2*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                        2*chi2*chi1*interp_spline(&integral[3].result, 0.0)
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
        if (r21 != 0 && r22 != 0){
            result +=
                /* constant in front */
                3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    chi2*interp_spline(&bg->f, z1_const)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                        (lambda2 - 6*chi1*costheta + 3*lambda2*(2*costheta*costheta - 1))
                       *interp_spline(&integral[0].result, sqrt(r21))/15.
                       -(
                            6*chi1*chi1*chi1*costheta - chi1*chi1*lambda2*(9*costheta*costheta + 11)
                           +chi1*lambda2*lambda2*costheta*(3*(2*costheta*costheta - 1) + 19)
                           -2*lambda2*lambda2*lambda2*(3*(2*costheta*costheta - 1) + 1)
                        )
                       *interp_spline(&integral[1].result, sqrt(r21))/r21/21.
                    )
                    +
                    chi1*interp_spline(&bg->f, z2_const)*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                        (lambda1 - 6*chi2*costheta + 3*lambda1*(2*costheta*costheta - 1))
                       *interp_spline(&integral[0].result, sqrt(r22))/15.
                       -(
                            6*chi2*chi2*chi2*costheta - chi2*chi2*lambda1*(9*costheta*costheta + 11)
                           +chi2*lambda1*lambda1*costheta*(3*(2*costheta*costheta - 1) + 19)
                           -2*lambda1*lambda1*lambda1*(3*(2*costheta*costheta - 1) + 1)
                        )
                       *interp_spline(&integral[1].result, sqrt(r22))/r22/21.
                    )
                );
            if (fabs(mu) < 0.999){
                result +=
                3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    chi2*interp_spline(&bg->f, z1_const)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                       -(
                           -4*pow(chi1, 5)*costheta
                           -pow(chi1, 3)*pow(lambda2, 2)*costheta*((2*costheta*costheta - 1) + 7)
                           +pow(chi1, 2)*pow(lambda2, 3)*(pow(costheta, 4) + 12*costheta*costheta - 21)
                           -3*chi1*pow(lambda2, 4)*costheta*((2*costheta*costheta - 1) - 5)
                           -pow(lambda2, 5)*(3*(2*costheta*costheta - 1) + 1)
                           +12*pow(chi1, 4)*lambda2
                        )
                       *interp_spline(&integral[2].result, sqrt(r21))/r21/r21/35.
                    )
                    +
                    chi1*interp_spline(&bg->f, z2_const)*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                       -(
                           -4*pow(chi2, 5)*costheta
                           -pow(chi2, 3)*pow(lambda1, 2)*costheta*((2*costheta*costheta - 1) + 7)
                           +pow(chi2, 2)*pow(lambda1, 3)*(pow(costheta, 4) + 12*costheta*costheta - 21)
                           -3*chi2*pow(lambda1, 4)*costheta*((2*costheta*costheta - 1) - 5)
                           -pow(lambda1, 5)*(3*(2*costheta*costheta - 1) + 1)
                           +12*pow(chi2, 4)*lambda1
                        )
                       *interp_spline(&integral[2].result, sqrt(r22))/r22/r22/35.
                    )
                );
            }
            else{
                result +=
                3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    chi2*interp_spline(&bg->f, z1_const)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *4.*(lambda2 + chi1)*interp_spline(&integral[2].result, sqrt(r21))/35.
                    +
                    chi1*interp_spline(&bg->f, z2_const)*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *4.*(lambda1 + chi2)*interp_spline(&integral[2].result, sqrt(r22))/35.
                );
            }
        }
        else if (r21 == 0 && r22 != 0){
            result +=
                3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    chi2*interp_spline(&bg->f, z1_const)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                    /* integrand */
                   *(1 - chi1/chi2)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                       -2*chi1*interp_spline(&integral[0].result, 0.0)/15.
                    )
                    +
                    chi1*interp_spline(&bg->f, z2_const)*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                        (lambda1 - 6*chi2*costheta + 3*lambda1*(2*costheta*costheta - 1))
                       *interp_spline(&integral[0].result, sqrt(r22))/15.
                       -(
                            6*chi2*chi2*chi2*costheta - chi2*chi2*lambda1*(9*costheta*costheta + 11)
                           +chi2*lambda1*lambda1*costheta*(3*(2*costheta*costheta - 1) + 19)
                           -2*lambda1*lambda1*lambda1*(3*(2*costheta*costheta - 1) + 1)
                        )
                       *interp_spline(&integral[1].result, sqrt(r22))/r22/21.
                       -(
                           -4*pow(chi2, 5)*costheta
                           -pow(chi2, 3)*pow(lambda1, 2)*costheta*((2*costheta*costheta - 1) + 7)
                           +pow(chi2, 2)*pow(lambda1, 3)*(pow(costheta, 4) + 12*costheta*costheta - 21)
                           -3*chi2*pow(lambda1, 4)*costheta*((2*costheta*costheta - 1) - 5)
                           -pow(lambda1, 5)*(3*(2*costheta*costheta - 1) + 1)
                           +12*pow(chi2, 4)*lambda1
                        )
                       *interp_spline(&integral[2].result, sqrt(r22))/r22/r22/35.
                    )
                );
        }
        else if (r21 != 0 && r22 == 0){
            result +=
                /* constant in front */
                3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    chi2*interp_spline(&bg->f, z1_const)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                        (lambda2 - 6*chi1*costheta + 3*lambda2*(2*costheta*costheta - 1))
                       *interp_spline(&integral[0].result, sqrt(r21))/15.
                       -(
                            6*chi1*chi1*chi1*costheta - chi1*chi1*lambda2*(9*costheta*costheta + 11)
                           +chi1*lambda2*lambda2*costheta*(3*(2*costheta*costheta - 1) + 19)
                           -2*lambda2*lambda2*lambda2*(3*(2*costheta*costheta - 1) + 1)
                        )
                       *interp_spline(&integral[1].result, sqrt(r21))/r21/21.
                       -(
                           -4*pow(chi1, 5)*costheta
                           -pow(chi1, 3)*pow(lambda2, 2)*costheta*((2*costheta*costheta - 1) + 7)
                           +pow(chi1, 2)*pow(lambda2, 3)*(pow(costheta, 4) + 12*costheta*costheta - 21)
                           -3*chi1*pow(lambda2, 4)*costheta*((2*costheta*costheta - 1) - 5)
                           -pow(lambda2, 5)*(3*(2*costheta*costheta - 1) + 1)
                           +12*pow(chi1, 4)*lambda2
                        )
                       *interp_spline(&integral[2].result, sqrt(r21))/r21/r21/35.
                    )
                    +
                    chi1*interp_spline(&bg->f, z2_const)*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
                    /* integrand */
                   *(1 - chi2/chi1)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                      -2*chi2*interp_spline(&integral[0].result, 0.0)/15.
                    )
                );
        }
        else{
            result +=
                3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    chi2*interp_spline(&bg->f, z1_const)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                    /* integrand */
                   *(1 - chi1/chi2)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                      -2*chi1*interp_spline(&integral[0].result, 0.0)/15.
                    )
                    +
                    chi1*interp_spline(&bg->f, z2_const)*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
                    /* integrand */
                   *(1 - chi2/chi1)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                       -2*chi2*interp_spline(&integral[0].result, 0.0)/15.
                    )
                );
        }
    }
    /* d1-len + len-d1 term */
    if (
        par->correlation_contrib.d1 &&
        par->correlation_contrib.len
    ){
        if (r21 != 0 && r22 != 0){
            result +=
                /* constant in front */
                3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    chi2*interp_spline(&bg->conformal_Hz, z1_const)*interp_spline(&bg->f, z1_const)
                   *interp_spline(&bg->G1, z1_const)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                        2*(costheta*(lambda2*lambda2 - 2*chi1*chi1) + chi1*lambda2*(2*(2*costheta*costheta - 1) - 1))
                       *interp_spline(&integral[3].result, sqrt(r21))/15.
                       +2*costheta*interp_spline(&integral[5].result, sqrt(r21))/3.
                       -(
                            4*pow(chi1, 4)*costheta
                           -pow(chi1, 3)*lambda2*(costheta*costheta + 9)
                           +chi1*chi1*lambda2*lambda2*costheta*(costheta*costheta + 5)
                           -2*chi1*pow(lambda2, 3)*((2*costheta*costheta - 1) - 2)
                           -2*pow(lambda2, 4)*costheta
                        )*interp_spline(&integral[4].result, sqrt(r21))/r21/15.
                    )
                    +
                    chi1*interp_spline(&bg->conformal_Hz, z2_const)*interp_spline(&bg->f, z2_const)
                   *interp_spline(&bg->G1, z2_const)*(2 - 5*s2)*interp_spline(&bg->D1, z2_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                        2*(costheta*(lambda1*lambda1 - 2*chi2*chi2) + chi2*lambda1*(2*(2*costheta*costheta - 1) - 1))
                       *interp_spline(&integral[3].result, sqrt(r22))/15.
                       +2*costheta*interp_spline(&integral[5].result, sqrt(r22))/3.
                       -(
                            4*pow(chi2, 4)*costheta
                           -pow(chi2, 3)*lambda1*(costheta*costheta + 9)
                           +chi2*chi2*lambda1*lambda1*costheta*(costheta*costheta + 5)
                           -2*chi2*pow(lambda1, 3)*((2*costheta*costheta - 1) - 2)
                           -2*pow(lambda1, 4)*costheta
                        )*interp_spline(&integral[4].result, sqrt(r22))/r22/15.
                    )
                );
        }
        else if (r21 == 0 && r22 != 0){
            result +=
                /* constant in front */
                3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    chi2*interp_spline(&bg->conformal_Hz, z1_const)*interp_spline(&bg->f, z1_const)
                   *interp_spline(&bg->G1, z1_const)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                       2*interp_spline(&integral[5].result, 0.0)/3.
                    )
                    +
                    chi1*interp_spline(&bg->conformal_Hz, z2_const)*interp_spline(&bg->f, z2_const)
                   *interp_spline(&bg->G1, z2_const)*(2 - 5*s2)*interp_spline(&bg->D1, z2_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                        2*(costheta*(lambda1*lambda1 - 2*chi2*chi2) + chi2*lambda1*(2*(2*costheta*costheta - 1) - 1))
                       *interp_spline(&integral[3].result, sqrt(r22))/15.
                       +2*costheta*interp_spline(&integral[5].result, sqrt(r22))/3.
                       -(
                            4*pow(chi2, 4)*costheta
                           -pow(chi2, 3)*lambda1*(costheta*costheta + 9)
                           +chi2*chi2*lambda1*lambda1*costheta*(costheta*costheta + 5)
                           -2*chi2*pow(lambda1, 3)*((2*costheta*costheta - 1) - 2)
                           -2*pow(lambda1, 4)*costheta
                        )*interp_spline(&integral[4].result, sqrt(r22))/r22/15.
                    )
                );
        }
        else if (r21 != 0 && r22 == 0){
            result +=
                /* constant in front */
                3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    chi2*interp_spline(&bg->conformal_Hz, z1_const)*interp_spline(&bg->f, z1_const)
                   *interp_spline(&bg->G1, z1_const)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                        2*(costheta*(lambda2*lambda2 - 2*chi1*chi1) + chi1*lambda2*(2*(2*costheta*costheta - 1) - 1))
                       *interp_spline(&integral[3].result, sqrt(r21))/15.
                       +2*costheta*interp_spline(&integral[5].result, sqrt(r21))/3.
                       -(
                            4*pow(chi1, 4)*costheta
                           -pow(chi1, 3)*lambda2*(costheta*costheta + 9)
                           +chi1*chi1*lambda2*lambda2*costheta*(costheta*costheta + 5)
                           -2*chi1*pow(lambda2, 3)*((2*costheta*costheta - 1) - 2)
                           -2*pow(lambda2, 4)*costheta
                        )*interp_spline(&integral[4].result, sqrt(r21))/r21/15.
                    )
                    +
                    chi1*interp_spline(&bg->conformal_Hz, z2_const)*interp_spline(&bg->f, z2_const)
                   *interp_spline(&bg->G1, z2_const)*(2 - 5*s2)*interp_spline(&bg->D1, z2_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                       2*interp_spline(&integral[5].result, 0.0)/3.
                    )
                );
        }
        else{
            result +=
                /* constant in front */
                3*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    chi2*interp_spline(&bg->conformal_Hz, z1_const)*interp_spline(&bg->f, z1_const)
                   *interp_spline(&bg->G1, z1_const)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                       2*interp_spline(&integral[5].result, 0.0)/3.
                    )
                    +
                    chi1*interp_spline(&bg->conformal_Hz, z2_const)*interp_spline(&bg->f, z2_const)
                   *interp_spline(&bg->G1, z2_const)*(2 - 5*s2)*interp_spline(&bg->D1, z2_const)
                    /* integrand */
                   *(1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                       2*interp_spline(&integral[5].result, 0.0)/3.
                    )
                );
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
                chi2*(3 - interp_spline(&par->evolution_bias1, z1_const))*interp_spline(&bg->f, z1_const)
               *pow(interp_spline(&bg->conformal_Hz, z1_const), 2)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
               *(
                    /* integrand */
                    (1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                        2*chi1*costheta*interp_spline(&integral[7].result, sqrt(r21))
                       -chi1*chi1*lambda2*(1 - costheta*costheta)*interp_spline(&integral[6].result, sqrt(r21))
                    )
                )
                +
                chi1*(3 - interp_spline(&par->evolution_bias2, z2_const))*interp_spline(&bg->f, z2_const)
               *pow(interp_spline(&bg->conformal_Hz, z2_const), 2)*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
               *(
                    /* integrand */
                    (1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                        2*chi2*costheta*interp_spline(&integral[7].result, sqrt(r22))
                       -chi2*chi2*lambda1*(1 - costheta*costheta)*interp_spline(&integral[6].result, sqrt(r22))
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
                chi2*(1 + interp_spline(&bg->G1, z1_const))*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
               *(
                    /* integrand */
                    (1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                        2*chi1*costheta*interp_spline(&integral[7].result, sqrt(r21))
                       -chi1*chi1*lambda2*(1 - costheta*costheta)*interp_spline(&integral[6].result, sqrt(r21))
                    )
                )
                +
                chi1*(1 + interp_spline(&bg->G2, z2_const))*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
               *(
                    /* integrand */
                    (1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                        2*chi2*costheta*interp_spline(&integral[7].result, sqrt(r22))
                       -chi2*chi2*lambda1*(1 - costheta*costheta)*interp_spline(&integral[6].result, sqrt(r22))
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
                chi2*(5*s1 - 2)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
               *(
                    /* integrand */
                    (1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                        2*chi1*costheta*interp_spline(&integral[7].result, sqrt(r21))
                       -chi1*chi1*lambda2*(1 - costheta*costheta)*interp_spline(&integral[6].result, sqrt(r21))
                    )
                )
                +
                chi1*(5*s2 - 2)*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
               *(
                    /* integrand */
                    (1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                        2*chi2*costheta*interp_spline(&integral[7].result, sqrt(r22))
                       -chi2*chi2*lambda1*(1 - costheta*costheta)*interp_spline(&integral[6].result, sqrt(r22))
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
                chi2*(interp_spline(&bg->f, z1_const) - 1)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
               *(
                    /* integrand */
                    (1 - x)*interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
                   *(
                        2*chi1*costheta*interp_spline(&integral[7].result, sqrt(r21))
                       -chi1*chi1*lambda2*(1 - costheta*costheta)*interp_spline(&integral[6].result, sqrt(r21))
                    )
                )
                +
                chi1*(interp_spline(&bg->f, z2_const) - 1)*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
               *(
                    /* integrand */
                    (1 - x)*interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
                   *(
                        2*chi2*costheta*interp_spline(&integral[7].result, sqrt(r22))
                       -chi2*chi2*lambda1*(1 - costheta*costheta)*interp_spline(&integral[6].result, sqrt(r22))
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
                b1*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                /* integrand */
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
               *interp_spline(&integral[5].result, sqrt(r21))
               +
                b2*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
                /* integrand */
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
               *interp_spline(&integral[5].result, sqrt(r22))
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
                chi2*b1*interp_spline(&bg->G2, z2_const)*interp_spline(&bg->D1, z1_const)
                /* integrand */
               *interp_spline(&bg->conformal_Hz, z2)*(interp_spline(&bg->f, z2) - 1)
               *interp_spline(&bg->D1, z2)*interp_spline(&bg->a, z2)
               *interp_spline(&integral[5].result, sqrt(r21))
               +
                chi1*b2*interp_spline(&bg->G1, z1_const)*interp_spline(&bg->D1, z2_const)
                /* integrand */
               *interp_spline(&bg->conformal_Hz, z1)*(interp_spline(&bg->f, z1) - 1)
               *interp_spline(&bg->D1, z1)*interp_spline(&bg->a, z1)
               *interp_spline(&integral[5].result, sqrt(r22))
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
                interp_spline(&bg->f, z1_const)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
                /* integrand */
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
               *(
                    (2*r21/3. + (costheta*costheta - 1)*lambda2*lambda2)
                   *interp_spline(&integral[6].result, sqrt(r21))
                   -interp_spline(&integral[5].result, sqrt(r21))/3.
                )
               +
                interp_spline(&bg->f, z2_const)*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
                /* integrand */
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
               *(
                    (2*r22/3. + (costheta*costheta - 1)*lambda1*lambda1)
                   *interp_spline(&integral[6].result, sqrt(r22))
                   -interp_spline(&integral[5].result, sqrt(r22))/3.
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
                chi2*interp_spline(&bg->f, z1_const)*interp_spline(&bg->G2, z2_const)*interp_spline(&bg->D1, z1_const)
                /* integrand */
               *interp_spline(&bg->conformal_Hz, z2)*(interp_spline(&bg->f, z2) - 1)
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
               *(
                    (2*r21/3. + (costheta*costheta - 1)*lambda2*lambda2)
                   *interp_spline(&integral[6].result, sqrt(r21))
                   -interp_spline(&integral[5].result, sqrt(r21))/3.
                )
               +
                chi1*interp_spline(&bg->f, z2_const)*interp_spline(&bg->G1, z1_const)*interp_spline(&bg->D1, z2_const)
                /* integrand */
               *interp_spline(&bg->conformal_Hz, z1)*(interp_spline(&bg->f, z1) - 1)
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
               *(
                    (2*r22/3. + (costheta*costheta - 1)*lambda1*lambda1)
                   *interp_spline(&integral[6].result, sqrt(r22))
                   -interp_spline(&integral[5].result, sqrt(r22))/3.
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
                interp_spline(&bg->conformal_Hz, z1_const)*interp_spline(&bg->f, z1_const)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)*(lambda2*costheta - chi1)
               *interp_spline(&integral[7].result, sqrt(r21))
               +
                interp_spline(&bg->conformal_Hz, z2_const)*interp_spline(&bg->f, z2_const)*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)*(lambda1*costheta - chi2)
               *interp_spline(&integral[7].result, sqrt(r22))
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
                chi2*interp_spline(&bg->conformal_Hz, z1_const)*interp_spline(&bg->f, z1_const)
               *interp_spline(&bg->G2, z2_const)*interp_spline(&bg->D1, z1_const)
               *interp_spline(&bg->conformal_Hz, z2)*(interp_spline(&bg->f, z2) - 1)
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)*(lambda2*costheta - chi1)
               *interp_spline(&integral[7].result, sqrt(r21))
               +
                chi1*interp_spline(&bg->conformal_Hz, z2_const)*interp_spline(&bg->f, z2_const)
               *interp_spline(&bg->G1, z1_const)*interp_spline(&bg->D1, z2_const)
               *interp_spline(&bg->conformal_Hz, z1)*(interp_spline(&bg->f, z1) - 1)
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)*(lambda1*costheta - chi2)
               *interp_spline(&integral[7].result, sqrt(r22))
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
                (3 - interp_spline(&par->evolution_bias1, z1_const))*interp_spline(&bg->f, z1_const)
               *pow(interp_spline(&bg->conformal_Hz, z1_const), 2)*(2 - 5*s2)*interp_spline(&bg->D1, z1_const)
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
               *ren1
               +
                (3 - interp_spline(&par->evolution_bias2, z2_const))*interp_spline(&bg->f, z2_const)
               *pow(interp_spline(&bg->conformal_Hz, z2_const), 2)*(2 - 5*s1)*interp_spline(&bg->D1, z2_const)
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
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
                chi2*(3 - interp_spline(&par->evolution_bias1, z1_const))*interp_spline(&bg->f, z1_const)
               *pow(interp_spline(&bg->conformal_Hz, z1_const), 2)*interp_spline(&bg->G2, z2_const)*interp_spline(&bg->D1, z1_const)
               *interp_spline(&bg->conformal_Hz, z2)*(interp_spline(&bg->f, z2) - 1)
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)
               *ren1
               +
                chi1*(3 - interp_spline(&par->evolution_bias2, z2_const))*interp_spline(&bg->f, z2_const)
               *pow(interp_spline(&bg->conformal_Hz, z2_const), 2)*interp_spline(&bg->G1, z1_const)*interp_spline(&bg->D1, z2_const)
               *interp_spline(&bg->conformal_Hz, z1)*(interp_spline(&bg->f, z1) - 1)
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)
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
                (1 + interp_spline(&bg->G1, z1_const))*(2 - 5*s2)
               *interp_spline(&bg->D1, z1_const)/interp_spline(&bg->a, z1_const)
                /* integrand */
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)*ren1
                +
                (1 + interp_spline(&bg->G2, z2_const))*(2 - 5*s1)
               *interp_spline(&bg->D1, z2_const)/interp_spline(&bg->a, z2_const)
                /* integrand */
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)*ren2
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
                chi2*(1 + interp_spline(&bg->G1, z1_const))*interp_spline(&bg->G2, z2_const)
               *interp_spline(&bg->D1, z1_const)/interp_spline(&bg->a, z1_const)
                /* integrand */
               *interp_spline(&bg->conformal_Hz, lambda2)
               *(interp_spline(&bg->f, lambda2) - 1)
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)*ren1
                +
                chi1*(1 + interp_spline(&bg->G2, z2_const))*interp_spline(&bg->G1, z1_const)
               *interp_spline(&bg->D1, z2_const)/interp_spline(&bg->a, z2_const)
                /* integrand */
                *interp_spline(&bg->conformal_Hz, lambda1)
               *(interp_spline(&bg->f, lambda1) - 1)
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)*ren2
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
               *interp_spline(&bg->D1, z1_const)/interp_spline(&bg->a, z1_const)
                /* integrand */
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)*ren1
                +
                (5*s2 - 2)*(2 - 5*s1)
               *interp_spline(&bg->D1, z2_const)/interp_spline(&bg->a, z2_const)
                /* integrand */
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)*ren2
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
                chi2*(5*s1 - 2)*interp_spline(&bg->G2, z2_const)
               *interp_spline(&bg->D1, z1_const)/interp_spline(&bg->a, z1_const)
                /* integrand */
               *interp_spline(&bg->conformal_Hz, lambda2)
               *(interp_spline(&bg->f, lambda2) - 1)
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)*ren1
                +
                chi1*(5*s2 - 2)*interp_spline(&bg->G1, z1_const)
               *interp_spline(&bg->D1, z2_const)/interp_spline(&bg->a, z2_const)
                /* integrand */
                *interp_spline(&bg->conformal_Hz, lambda1)
               *(interp_spline(&bg->f, lambda1) - 1)
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)*ren2
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
                (interp_spline(&bg->f, z1_const) - 1)*(2 - 5*s2)
               *interp_spline(&bg->D1, z1_const)/interp_spline(&bg->a, z1_const)
                /* integrand */
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)*ren1
                +
                (interp_spline(&bg->f, z2_const) - 1)*(2 - 5*s1)
               *interp_spline(&bg->D1, z2_const)/interp_spline(&bg->a, z2_const)
                /* integrand */
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)*ren2
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
                chi2*(interp_spline(&bg->f, z1_const) - 1)*interp_spline(&bg->G2, z2_const)
               *interp_spline(&bg->D1, z1_const)/interp_spline(&bg->a, z1_const)
                /* integrand */
               *interp_spline(&bg->conformal_Hz, lambda2)
               *(interp_spline(&bg->f, lambda2) - 1)
               *interp_spline(&bg->D1, z2)/interp_spline(&bg->a, z2)*ren1
                +
                chi1*(interp_spline(&bg->f, z2_const) - 1)*interp_spline(&bg->G1, z1_const)
               *interp_spline(&bg->D1, z2_const)/interp_spline(&bg->a, z2_const)
                /* integrand */
                *interp_spline(&bg->conformal_Hz, lambda1)
               *(interp_spline(&bg->f, lambda1) - 1)
               *interp_spline(&bg->D1, z1)/interp_spline(&bg->a, z1)*ren2
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
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t integral[],
    double z_mean,
    double mu,
    double sep,
    double x1,
    double x2
)
{
    double result = 0;

    double chi_mean = interp_spline(&bg->comoving_distance, z_mean);
    double chi1 = chi_mean - sep*mu/2.;
    double chi2 = chi_mean + sep*mu/2.;
    double costheta =
        (2*chi_mean*chi_mean - sep*sep + mu*mu*sep*sep/2.)
       /(2*chi_mean*chi_mean - mu*mu*sep*sep/2.);
    double lambda1 = chi1*x1, lambda2 = chi2*x2;
    double r2 = lambda1*lambda1 + lambda2*lambda2 - 2*lambda1*lambda2*costheta;
    if (r2 < 0) r2 = 0;

    double z1_const = interp_spline(&bg->z_as_chi, chi1);
    double z2_const = interp_spline(&bg->z_as_chi, chi2);
    double z1 = interp_spline(&bg->z_as_chi, lambda1);
    double z2 = interp_spline(&bg->z_as_chi, lambda2);

    double s1 = interp_spline(&par->magnification_bias1, z1_const);
    double s2 = interp_spline(&par->magnification_bias2, z2_const);
    double sz_mean1 = interp_spline(&par->magnification_bias1, z_mean);
    double sz_mean2 = interp_spline(&par->magnification_bias2, z_mean);

    double ren = 0;
    if (par->divergent){
        if (r2 <= pow(0.000001*COFFE_H0, 2)){
            ren = interp_spline(&integral[8].renormalization0, lambda1);
        }
        else{
            ren = interp_spline(&integral[8].result, sqrt(r2))
                    /* renormalization term */
                   -gsl_spline2d_eval(
                        integral[8].renormalization.spline,
                        lambda1, lambda2,
                        integral[8].renormalization.xaccel,
                        integral[8].renormalization.yaccel
                    );
        }
    }

    /* len-len term */
    if (par->correlation_contrib.len){
        if (par->flatsky){
            double temp_sep = chi_mean*sqrt(2.*(1. - costheta));
            /* len-len modified by flatsky */
            result +=
            /* constant in front */
            9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)*(2 - 5*sz_mean1)*(2 - 5*sz_mean2)/8./M_PI*pow(chi_mean, 3)
            /* integrand */
           *interp_spline(&integral[9].result, x1*temp_sep*sqrt(1 - mu*mu))
           *pow(interp_spline(&bg->D1, interp_spline(&bg->z_as_chi, x1*chi_mean)), 2) // D(lambda)^2
           *pow(1 + interp_spline(&bg->z_as_chi, x1*chi_mean), 2) // (1 + z(lambda))^2
           *pow(x1*(1 - x1), 2);
        }
        else{
        if (r2 > 1e-20){
            result +=
            /* constant in front */
            9.*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)*(2 - 5*s1)*(2 - 5*s2)/4.*chi1*chi2
           *
            /* integrand */
            interp_spline(&bg->D1, z1)
           *interp_spline(&bg->D1, z2)
           /interp_spline(&bg->a, z1)
           /interp_spline(&bg->a, z2)
           *(1 - x1)*(1 - x2)
           *(
                2*(costheta*costheta - 1)*lambda1*lambda2
               *interp_spline(&integral[0].result, sqrt(r2))/5.
               +
                4*costheta
               *interp_spline(&integral[5].result, sqrt(r2))/3.
               +
                4*costheta*(r2 + 6*costheta*lambda1*lambda2)
               *interp_spline(&integral[3].result, sqrt(r2))/15.
               +
                2*(costheta*costheta - 1)*lambda1*lambda2
               *(2*r2 + 3*costheta*lambda1*lambda2)
               *interp_spline(&integral[1].result, sqrt(r2))/7./r2
               +
                2*costheta
               *(2*r2*r2 + 12*costheta*r2*lambda1*lambda2 + 15*(costheta*costheta - 1)*lambda1*lambda1*lambda2*lambda2)
               *interp_spline(&integral[4].result, sqrt(r2))/15./r2
               +
                (costheta*costheta - 1)*lambda1*lambda2
               *(6*r2*r2 + 30*costheta*r2*lambda1*lambda2 + 35*(costheta*costheta - 1)*lambda1*lambda1*lambda2*lambda2)
               *interp_spline(&integral[2].result, sqrt(r2))/35./r2/r2
            );
        }
        else{
            result +=
            /* constant in front */
            9./4*pow((par->Omega0_cdm + par->Omega0_baryon), 2)*(2 - 5*s1)*(2 - 5*s2)*chi1*chi2
           *
            /* integrand */
            interp_spline(&bg->D1, z1)
           *interp_spline(&bg->D1, z2)
           /interp_spline(&bg->a, z1)
           /interp_spline(&bg->a, z2)
           *(1 - x1)*(1 - x2)
           *(
               4*interp_spline(&integral[5].result, 0.0)/3.
               +
                24.*lambda1*lambda2
               *interp_spline(&integral[3].result, 0.0)/15.
            );
        }
        }
    }
    /* g4-g4 term */
    if (par->correlation_contrib.g4){
        result +=
        /* constant in front */
        9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)*(2 - 5*s1)*(2 - 5*s2)
       *
            /* integrand */
            interp_spline(&bg->D1, z1)
           *interp_spline(&bg->D1, z2)
           /interp_spline(&bg->a, z1)
           /interp_spline(&bg->a, z2)
           *ren;
    }
    /* g5-g5 term */
    if (par->correlation_contrib.g5){
        result +=
        /* constant in front */
        9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)
       *interp_spline(&bg->G1, z1_const)
       *interp_spline(&bg->G2, z2_const)
       *chi1*chi2
       *
        /* integrand */
            interp_spline(&bg->D1, z1)
           *interp_spline(&bg->D1, z2)
           /interp_spline(&bg->a, z1)
           /interp_spline(&bg->a, z2)
           *interp_spline(&bg->conformal_Hz, z1)
           *interp_spline(&bg->conformal_Hz, z2)
           *(interp_spline(&bg->f, z1) - 1)
           *(interp_spline(&bg->f, z2) - 1)
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
                   *(1 - x2)/x2*interp_spline(&bg->D1, z1)*interp_spline(&bg->D1, z2)
                   /interp_spline(&bg->a, z1)/interp_spline(&bg->a, z2)
                   *(
                        2*lambda1*lambda2*costheta*interp_spline(&integral[7].result, sqrt(r2))
                       -lambda1*lambda1*lambda2*lambda2*(1 - costheta*costheta)*interp_spline(&integral[6].result, sqrt(r2))
                    )
                    +
                    (2 - 5*s1)*(2 - 5*s2)
                   *(1 - x1)/x1*interp_spline(&bg->D1, z1)*interp_spline(&bg->D1, z2)
                   /interp_spline(&bg->a, z1)/interp_spline(&bg->a, z2)
                   *(
                        2*lambda1*lambda2*costheta*interp_spline(&integral[7].result, sqrt(r2))
                       -lambda1*lambda1*lambda2*lambda2*(1 - costheta*costheta)*interp_spline(&integral[6].result, sqrt(r2))
                    )
                );
        }
        else{
            result +=
                9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    (2 - 5*s1)*(2 - 5*s2)
                   *(1 - x2)/x2*interp_spline(&bg->D1, z1)*interp_spline(&bg->D1, z2)
                   /interp_spline(&bg->a, z1)/interp_spline(&bg->a, z2)
                   *2*lambda1*lambda2*interp_spline(&integral[7].result, 0.0)
                    +
                    (2 - 5*s1)*(2 - 5*s2)
                   *(1 - x1)/x1*interp_spline(&bg->D1, z2)*interp_spline(&bg->D1, z1)
                   /interp_spline(&bg->a, z2)/interp_spline(&bg->a, z1)
                   *2*lambda1*lambda2*interp_spline(&integral[7].result, 0.0)
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
                    (2 - 5*s2)*interp_spline(&bg->G1, z1_const)*chi1
                   *interp_spline(&bg->conformal_Hz, z1)*(interp_spline(&bg->conformal_Hz, z1) - 1)
                   *(1 - x2)/x2*interp_spline(&bg->D1, z1)*interp_spline(&bg->D1, z2)
                   /interp_spline(&bg->a, z1)/interp_spline(&bg->a, z2)
                   *(
                        2*lambda1*lambda2*costheta*interp_spline(&integral[7].result, sqrt(r2))
                       -lambda1*lambda1*lambda2*lambda2*(1 - costheta*costheta)*interp_spline(&integral[6].result, sqrt(r2))
                    )
                    +
                    (2 - 5*s1)*interp_spline(&bg->G2, z2_const)*chi2
                   *interp_spline(&bg->conformal_Hz, z2)*(interp_spline(&bg->conformal_Hz, z2) - 1)
                   *(1 - x1)/x1*interp_spline(&bg->D1, z1)*interp_spline(&bg->D1, z2)
                   /interp_spline(&bg->a, z1)/interp_spline(&bg->a, z2)
                   *(
                        2*lambda1*lambda2*costheta*interp_spline(&integral[7].result, sqrt(r2))
                       -lambda1*lambda1*lambda2*lambda2*(1 - costheta*costheta)*interp_spline(&integral[6].result, sqrt(r2))
                    )
                );
        }
        else{
            result +=
                9*(par->Omega0_cdm + par->Omega0_baryon)*(par->Omega0_cdm + par->Omega0_baryon)/2.
               *(
                    (2 - 5*s2)*interp_spline(&bg->G1, z1_const)*chi1
                   *interp_spline(&bg->conformal_Hz, z1)*(interp_spline(&bg->conformal_Hz, z1) - 1)
                   *(1 - x2)/x2*interp_spline(&bg->D1, z1)*interp_spline(&bg->D1, z2)
                   /interp_spline(&bg->a, z1)/interp_spline(&bg->a, z2)
                   *2*lambda1*lambda2*interp_spline(&integral[7].result, 0.0)
                    +
                    (2 - 5*s1)*interp_spline(&bg->G2, z2_const)*chi2
                   *interp_spline(&bg->conformal_Hz, z2)*(interp_spline(&bg->conformal_Hz, z2) - 1)
                   *(1 - x1)/x1*interp_spline(&bg->D1, z1)*interp_spline(&bg->D1, z2)
                   /interp_spline(&bg->a, z1)/interp_spline(&bg->a, z2)
                   *2*lambda1*lambda2*interp_spline(&integral[7].result, 0.0)
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
                interp_spline(&bg->G2, z2_const)*(2 - 5*s1)*chi2
               *interp_spline(&bg->conformal_Hz, z2)*(interp_spline(&bg->f, z2) - 1)
               *interp_spline(&bg->D1, z1)*interp_spline(&bg->D1, z2)
               /interp_spline(&bg->a, z1)/interp_spline(&bg->a, z2)
               *ren
               +
                interp_spline(&bg->G1, z1_const)*(2 - 5*s2)*chi1
               *interp_spline(&bg->conformal_Hz, z1)*(interp_spline(&bg->f, z1) - 1)
               *interp_spline(&bg->D1, z1)*interp_spline(&bg->D1, z2)
               /interp_spline(&bg->a, z1)/interp_spline(&bg->a, z2)
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
