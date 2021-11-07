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

#ifdef __cplusplus
#include <complex>
typedef std::complex<double> complex_t;
using namespace std::complex_literals;
#define I 1i
#include "compatibility.h"
#else
#include <complex.h>
typedef double complex complex_t;
#endif

#include <math.h>
#include <fftw3.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#include "twofast.h"

/* for C99 compatibility */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PIl
#define M_PIl 3.141592653589793238462643383279502884L
#endif

#ifndef SQRT_PI
#define SQRT_PI 1.7724538509055160272981674833411451827975494561224
#endif

static double twofast_window(
    double value,
    double xmin,
    double xmax,
    double xleft,
    double xright
)
{
    if (xmin <= xleft && xleft <= xright && xright <= xmax){
        double result;
        if (value > xleft && value < xright && value > xmin && value < xmax){
            result = 1.0;
        }
        else if (value <= xmin || value >= xmax){
            result = 0.0;
        }
        else{
            if (value < xleft && value > xmin){
                result = (value - xmin)/(xleft - xmin);
            }
            else if (value > xright && value < xmax){
                result = (xmax - value)/(xmax - xright);
            }
            result = result - sin(2*M_PI*result)/2./M_PI;
        }
        return result;
    }
    else{
        fprintf(stderr, "ERROR: window function received undefined limits!\n");
        fprintf(stderr,
            "input values are: "
            "xmin = %e, xleft = %e "
            "xright = %e, xmax = %e\n",
            xmin, xleft, xright, xmax);
        exit(EXIT_FAILURE);
    }
}

/**
    Copyright: Wikipedia, Creative Commons Attribution-ShareAlike 3.0 Unported License
    https://en.wikipedia.org/wiki/Lanczos_approximation
**/
// TODO check if this lgamma has the same accuracy as ac_lgamma
static complex_t twofast_lgamma(complex_t z)
{
    static const long double gamma_coeff[] = {
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };
    static const size_t gamma_coeff_size =
        sizeof(gamma_coeff) / sizeof(gamma_coeff[0]);
    complex_t result;
    if (creal(z) < 0.5){
        result = log(M_PI) - clog(csin(M_PI * z)) - twofast_lgamma(1. - z);
    }
    else{
        complex_t x = 0.99999999999980993;
        z -= 1;
        for (size_t i = 0; i < gamma_coeff_size; ++i)
            x += (complex_t)gamma_coeff[i] / (z + (double)i + 1.);
        complex_t t = z + (double)gamma_coeff_size - 0.5;
        result = log(2 * M_PI) / 2. + (z + 0.5) * clog(t) - t + clog(x);
        if (!(gsl_finite(creal(result)) && gsl_finite(cimag(result)))){
            fprintf(
                stderr,
                "ERROR: file %s, function %s\n", __FILE__, __func__
            );
            fprintf(
                stderr,
                "Values: "
                "z = %.3e + %.3eI, "
                "x = %.3e + %.3eI, "
                "t = %.3e + %.3eI, "
                "result = %.3e + %.3eI\n",
                creal(z), cimag(z),
                creal(x), cimag(x),
                creal(t), cimag(t),
                creal(result), cimag(result)
            );
            exit(EXIT_FAILURE);
        }
    }

    return result;
}

static complex_t twofast_gamma(complex_t z)
{
    static const long double gamma_coeff[] = {
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };
    static const size_t gamma_coeff_size =
        sizeof(gamma_coeff) / sizeof(gamma_coeff[0]);
    complex_t result;
    if (creal(z) < 0.5){
        result = M_PI / csin(M_PI * z) / twofast_gamma(1. - z);
    }
    else{
        complex_t x = 0.99999999999980993;
        z -= 1;
        for (size_t i = 0; i < gamma_coeff_size; ++i)
            x += (complex_t)gamma_coeff[i] / (z + (double)i + 1.);
        complex_t t = z + (double)gamma_coeff_size - 0.5;
        result = csqrt(2 * M_PI) * cpow(t, z + 0.5) / cexp(t) * x;
        if (!(gsl_finite(creal(result)) && gsl_finite(cimag(result)))){
            fprintf(
                stderr,
                "ERROR: file %s, function %s\n", __FILE__, __func__
            );
            fprintf(
                stderr,
                "Values: "
                "z = %.3e + %.3eI, "
                "x = %.3e + %.3eI, "
                "t = %.3e + %.3eI, "
                "result = %.3e + %.3eI\n",
                creal(z), cimag(z),
                creal(x), cimag(x),
                creal(t), cimag(t),
                creal(result), cimag(result)
            );
            exit(EXIT_FAILURE);
        }
    }

    return result;
}


static complex_t twofast_mql(
    double t,
    double q,
    double l,
    double alpha
)
{
    const complex_t n = q - 1. - t * I;
    complex_t unl =
        cpow(2, n - 1.)
       *SQRT_PI
       *cexp(
            twofast_lgamma((1. + (double)l + n) / 2.) - twofast_lgamma((2. + (double)l - n) / 2.)
        );
    if (gsl_finite(creal(unl)) && gsl_finite(cimag(unl))){
        return cpow(alpha, t * I - q) * unl;
    }
    else{
        fprintf(
            stderr,
            "ERROR: "
            "file %s, function %s\n"
            "inputs: t = %e, q = %e, l = %e, alpha = %e\n",
            __FILE__, __func__,
            t, q, l, alpha
        );
        exit(EXIT_FAILURE);
    }
}

static inline double twofast_min(double a, double b)
{
    if (a > b)
        return b;
    return a;
}

static inline double twofast_max(double a, double b)
{
    if (a < b)
        return b;
    return a;
}

static double twofast_selectqnu(double l, double nu)
{
    static const double n1 = 0.9, n2 = 0.9999;
    double qmin = twofast_max(n2 - 1.0 - nu, -l);
    double qmax = twofast_min(n1 + 3.0 - nu, 2.0);
    double qbest = (2. + n1 + n2)/2. - nu;
    double q = qbest;
    if (!(qmin < q && q < qmax))
        q = (qmin + 2.*qmax)/3.;
    return q;
}

static void twofast_fft_input(
    fftw_complex *output_y,
    size_t output_len,
    gsl_spline *spline,
    gsl_interp_accel *accel,
    double q,
    double k0,
    double kmin,
    double kmax,
    unsigned flag
)
{
    const size_t N2 = output_len/2 + 1;
    const double L = 2*M_PI*output_len/log(kmax/kmin);
    double *input_x_mod = (double *)malloc(sizeof(double)*output_len);
    double *input_y_mod = (double *)fftw_malloc(sizeof(double)*output_len);
    fftw_complex *input_y_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N2);
    fftw_plan p;
    #pragma omp critical
    {
        p = fftw_plan_dft_r2c_1d(output_len, input_y_mod, input_y_fft, flag);
    }

    for (size_t i = 0; i<output_len; ++i){
        input_x_mod[i] = k0*pow(kmax/kmin, (double)i/output_len);
    }

    for (size_t i = 0; i<output_len; ++i){
        input_y_mod[i] =
            pow(kmax/kmin, (3. - q)*i/output_len)
           *gsl_spline_eval(spline, k0*pow(kmax/kmin, (double)i/output_len), accel)
           *twofast_window(
                k0*pow(kmax/kmin, (double)i/output_len),
                kmin,
                k0*pow(kmax/kmin, (double)(output_len - 1)/output_len),
                exp(0.46)*kmin,
                exp(-0.46)*k0*pow(kmax/kmin, (double)(output_len - 1)/output_len)
            );
    }

    fftw_execute(p);
    #pragma omp critical
    {
        fftw_destroy_plan(p);
    }

    for (size_t i = 0; i<N2; ++i){
#ifdef __cplusplus
        output_y[i][0] =
            creal(
                twofast_window(
                    input_x_mod[N2 - 2 + i],
                    kmin,
                    k0*pow(kmax/kmin, (double)(output_len - 1)/output_len),
                    exp(0.46)*kmin,
                    exp(-0.46)*k0*pow(kmax/kmin, (double)(output_len - 1)/output_len)
                )
            )
           *input_y_fft[i][0]/L
            +
            cimag(
                twofast_window(
                    input_x_mod[N2 - 2 + i],
                    kmin,
                    k0*pow(kmax/kmin, (double)(output_len - 1)/output_len),
                    exp(0.46)*kmin,
                    exp(-0.46)*k0*pow(kmax/kmin, (double)(output_len - 1)/output_len)
                )
            )
           *input_y_fft[i][1]/L;
        output_y[i][1] =
            -creal(
                twofast_window(
                    input_x_mod[N2 - 2 + i],
                    kmin,
                    k0*pow(kmax/kmin, (double)(output_len - 1)/output_len),
                    exp(0.46)*kmin,
                    exp(-0.46)*k0*pow(kmax/kmin, (double)(output_len - 1)/output_len)
                )
            )
           *input_y_fft[i][1]/L
            +
            cimag(
                twofast_window(
                    input_x_mod[N2 - 2 + i],
                    kmin,
                    k0*pow(kmax/kmin, (double)(output_len - 1)/output_len),
                    exp(0.46)*kmin,
                    exp(-0.46)*k0*pow(kmax/kmin, (double)(output_len - 1)/output_len)
                )
            )
           *input_y_fft[i][0]/L;
#else
        output_y[i] =
            twofast_window(
                input_x_mod[N2 - 2 + i],
                kmin,
                k0*pow(kmax/kmin, (double)(output_len - 1)/output_len),
                exp(0.46)*kmin,
                exp(-0.46)*k0*pow(kmax/kmin, (double)(output_len - 1)/output_len)
            )
           *conj(input_y_fft[i])/L;
#endif
    }

    free(input_x_mod);
    fftw_free(input_y_mod);
    fftw_free(input_y_fft);
}

void twofast_1bessel(
    double *output_x,
    double *output_y,
    size_t output_len,
    double *input_x,
    double *input_y,
    size_t input_len,
    double l,
    double nu,
    double r0,
    double k0,
    double kmin,
    double kmax,
    unsigned flag
)
{
    switch(flag){
        case 0:
            flag = FFTW_ESTIMATE;
            break;
        case 1:
            flag = FFTW_MEASURE;
            break;
        case 2:
            flag = FFTW_PATIENT;
            break;
        case 3:
            flag = FFTW_EXHAUSTIVE;
            break;
        default:
            flag = FFTW_ESTIMATE;
            break;
    }
    const size_t N2 = output_len/2 + 1;
    const double qnu = twofast_selectqnu(l, nu);
    const double G = log(kmax/kmin);

    gsl_spline *input_spline = gsl_spline_alloc(gsl_interp_cspline, input_len);
    gsl_interp_accel *input_accel = gsl_interp_accel_alloc();
    gsl_spline_init(input_spline, input_x, input_y, input_len);
    fftw_complex *input_y_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N2);

    twofast_fft_input(
        input_y_fft, output_len,
        input_spline, input_accel,
        qnu + nu, k0, kmin, kmax,
        flag
    );
    gsl_spline_free(input_spline);
    gsl_interp_accel_free(input_accel);

    double *prefactors = (double *)fftw_malloc(sizeof(double)*output_len);

    for (size_t i = 0; i<output_len; ++i){
        output_x[i] = r0*pow(kmax/kmin, (double)i/output_len);
    }

    for (size_t i = 0; i<output_len; ++i){
        prefactors[i] =
            k0*k0*k0*pow(kmax/kmin, -((qnu + nu)*i/output_len))/M_PI/pow(r0*k0, nu)/G;
    }

    fftw_complex *temp_input = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N2);

    for (size_t i = 0; i<N2; ++i){
        // copying the input as it's destroyed by the planner
#ifdef __cplusplus
        temp_input[i][0] =
            input_y_fft[i][0]
           *creal(twofast_mql(2*M_PI*i/G, qnu, l, k0*r0))
            -
            input_y_fft[i][1]
           *cimag(twofast_mql(2*M_PI*i/G, qnu, l, k0*r0));
        temp_input[i][1] =
            input_y_fft[i][0]
           *cimag(twofast_mql(2*M_PI*i/G, qnu, l, k0*r0))
            +
            input_y_fft[i][1]
           *creal(twofast_mql(2*M_PI*i/G, qnu, l, k0*r0));
#else
        temp_input[i] =
            input_y_fft[i]
           *twofast_mql(2*M_PI*i/G, qnu, l, k0*r0);
#endif
    }

    double *temp_output_y = (double *)fftw_malloc(sizeof(double)*output_len);

    fftw_plan p;
    #pragma omp critical
    {
        p = fftw_plan_dft_c2r_1d(output_len, input_y_fft, temp_output_y, flag);
    }

    for (size_t i = 0; i<N2; ++i){
#ifdef __cplusplus
        input_y_fft[i][0] = temp_input[i][0];
        input_y_fft[i][1] = temp_input[i][1];
#else
        input_y_fft[i] = temp_input[i]; // putting it back
#endif
    }

    fftw_execute(p);
    #pragma omp critical
    {
        fftw_destroy_plan(p);
    }

    for (size_t i = 0; i<output_len; ++i){
        temp_output_y[i] *= prefactors[i];
    }

    for (size_t i = 0; i<output_len; ++i){
        output_y[i] = temp_output_y[i];
    }

    fftw_free(prefactors);
    fftw_free(input_y_fft);
    fftw_free(temp_input);
    fftw_free(temp_output_y);
}


#ifdef HAVE_ARB


static int minus_one(
    const int exponent
)
{
    const int value = exponent % 2 ? 1 : -1;
    return value;
}

// only 00 for now
static complex_t power_bessel2(
    const complex_t z,
    const double a,
    const double b,
    const int l1,
    const int l2
)
{
    const complex_t mu = z - l1 - l2 - 1;
    return 2 * minus_one(l1 + l2)
        * twofast_gamma(mu)
        * ccos(M_PI / 2. * mu)
        * (1. / a / b * (1. / cpow(a + b, mu) - 1. / cpow(a - b, mu)));
}

static complex_t twofast_mql1l2(double t, double r, double q, int l1, int l2, double alpha)
{
    const complex_t n = q - 1. - t*I;
    complex_t ul1l2;
    static const double limit = 1E-8;
    if (r < 1 - limit){
        ul1l2 =
            cpow(alpha, t*I - q) * power_bessel2(n, 1., r, l1, l2);
    }
    else if (r > 1 + limit){
        // TODO
    }
    else{
        ul1l2 =
            cpow(alpha, t*I - q)
           *cpow(2, n - 2)*M_PI
           *cexp(
                twofast_lgamma((1. + l1 + l2 + n)/2.)
              - twofast_lgamma((2. + l1 - l2 - n)/2.)
              + twofast_lgamma(1. - n)
              - twofast_lgamma((3. + l2 + l1 - n)/2.)
              - twofast_lgamma((2. + l2 - l1 - n)/2.)
            );
    }
    if (gsl_finite(ul1l2)) return ul1l2;
    else{
        fprintf(stderr, "ERROR: file %s, function %s\n", __FILE__, __func__);
        exit(EXIT_FAILURE);
    }
}

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
    double r,
    double r0,
    double k0,
    double kmin,
    double kmax,
    unsigned flag
)
{
    switch(flag){
        case 0:
            flag = FFTW_ESTIMATE;
            break;
        case 1:
            flag = FFTW_MEASURE;
            break;
        case 2:
            flag = FFTW_PATIENT;
            break;
        case 3:
            flag = FFTW_EXHAUSTIVE;
            break;
        default:
            flag = FFTW_ESTIMATE;
            break;
    }

    const size_t N2 = output_len/2 + 1;
    const double G = log(kmax/kmin);

    gsl_spline *input_spline = gsl_spline_alloc(gsl_interp_cspline, input_len);
    gsl_interp_accel *input_accel = gsl_interp_accel_alloc();
    gsl_spline_init(input_spline, input_x, input_y, input_len);
    fftw_complex *input_y_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N2);

    twofast_fft_input(
        input_y_fft, output_len,
        input_spline, input_accel,
        q, k0, kmin, kmax,
        flag
    );
    gsl_spline_free(input_spline);
    gsl_interp_accel_free(input_accel);

    double *prefactors = (double *)malloc(sizeof(double)*output_len);

    for (size_t i = 0; i<output_len; ++i){
        output_x[i] = r0*pow(kmax/kmin, (double)i/output_len);
    }

    for (size_t i = 0; i<output_len; ++i){
        prefactors[i] =
            4*k0*k0*k0*pow(kmax/kmin, -q*i/output_len)/G;
    }

    fftw_complex *temp_input = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N2);
    double *t_array = (double *)malloc(sizeof(double)*N2);

    for (size_t i = 0; i<N2; ++i){
        t_array[i] = 2*M_PI*i/G;
    }

    #pragma omp parallel for
    for (size_t i = 0; i<N2; ++i){
        // copying the input as it's destroyed by the planner
        temp_input[i] =
            input_y_fft[i]
           *twofast_mql1l2(2*M_PI*i/G, r, q, l1, l2, k0*r0);
    }

    double *temp_output_y = (double *)fftw_malloc(sizeof(double)*output_len);

    fftw_plan p =
        fftw_plan_dft_c2r_1d(output_len, input_y_fft, temp_output_y, flag);

    for (size_t i = 0; i<N2; ++i){
        input_y_fft[i] = temp_input[i]; // putting it back
    }

    fftw_execute(p);
    fftw_destroy_plan(p);

    for (size_t i = 0; i<output_len; ++i){
        temp_output_y[i] *= prefactors[i];
    }
    for (size_t i = 0; i<output_len; ++i){
        output_y[i] = temp_output_y[i];
    }

    fftw_free(prefactors);
    fftw_free(input_y_fft);
    fftw_free(temp_input);
    free(t_array);
    fftw_free(temp_output_y);
}
#endif
