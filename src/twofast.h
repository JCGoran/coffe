/*
 * This file is part of COFFE
 * Copyright (C) 2018 Goran Jelic-Cizmek
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

#ifndef TWOFAST_H
#define TWOFAST_H

#ifdef __cplusplus
extern "C" {
#endif


/*****
    computes the integral of the form k^2 P(k) j_l(k r)/(k r)^nu/(2 pi^2)
    INPUT:
        output_x - pointer to output x-array (MUST BE ALLOCATED BEFOREHAND)
        output_y - pointer to output y-array (MUST BE ALLOCATED BEFOREHAND USING FFTW_MALLOC)
        output_len - length of previous 2 arrays
        input_x - pointer to input x-array (k in above nomenclature)
        input_y - pointer to input y-array (P(k) in above nomenclature)
        input_len - length of previous 2 arrays
        l - degree of spherical bessel function
        nu - real number
        r0 - smallest separation for the output
        k0 - smallest value of input x-array
        kmin - smallest value of input x-array
        kmax - largest value of input y-array
        flag - FFTW transformation flag (0 is usually sufficient)
*****/
void twofast_1bessel(
    double *output_x,
    double *output_y,
    size_t output_len,
    double *input_x,
    double *input_y,
    size_t input_len,
    int l,
    double nu,
    double r0,
    double k0,
    double kmin,
    double kmax,
    unsigned flag
);

/**** DEPRACATED ****/
/*****
    computes the integral of the form 2/pi k^2 P(k) j_l1(k chi1) j_l2(k chi2)
    INPUT:
        output_x - pointer to output x-array (MUST BE ALLOCATED BEFOREHAND)
        output_y - pointer to output y-array (MUST BE ALLOCATED BEFOREHAND USING FFTW_MALLOC)
        output_len - length of previous 2 arrays
        input_x - pointer to input x-array (k in above nomenclature)
        input_y - pointer to input y-array (P(k) in above nomenclature)
        input_len - length of previous 2 arrays
        q - biasing parameter; a value of 1.1 usually works best
        l1 - degree of first spherical bessel function
        l2 - degree of second spherical bessel function
        R - ratio between chi1 and chi2 := R chi1
        r0 - smallest separation for the output
        k0 - smallest value of input x-array
        kmin - smallest value of input x-array
        kmax - largest value of input y-array
        flag - FFTW transformation flag (0 is usually sufficient)
    NOTE: this function can only be used if the user links the following:
          1) the FLINT library (http://www.flintlib.org/)
          2) the ARB library (http://arblib.org/)
          3) the C99 wrapper for ARB (available at https://github.com/fredrik-johansson/arbcmath)
*****/
#ifdef HAVE_ARB
void twofast_2bessel(
    double *output_x,
    double *output_y,
    size_t output_len,
    double *input_x,
    double *input_y,
    size_t input_len,
    double q,
    int l1,
    int l2,
    double R,
    double r0,
    double k0,
    double kmin,
    double kmax,
    unsigned flag
);
#endif
#ifdef __cplusplus
}
#endif

#endif
