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

#ifndef COFFE_COMMON_H
#define COFFE_COMMON_H

#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

#ifndef COFFE_CYTHON
#include <libconfig.h>
#endif

#ifdef HAVE_CUBA
#include "cuba.h"
#endif

#include "tanhsinh.h"

#ifndef COFFE_COPYRIGHT
#define COFFE_COPYRIGHT  \
                        "Copyright (C) 2019 Goran Jelic-Cizmek\n" \
                        "This is free software; see the source for copying conditions.\n" \
                        "There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n" \
                        "Report bugs to: goran.jelic-cizmek@unige.ch\n"
#endif


#ifndef COFFE_MAX_ALLOCABLE
#define COFFE_MAX_ALLOCABLE 1024 // largest allowed number of columns
#endif

#ifndef COFFE_MAX_STRLEN
#define COFFE_MAX_STRLEN 256 // largest allowed length of the filename
#endif

#ifndef COFFE_TRUE
#define COFFE_TRUE 1
#endif

#ifndef COFFE_FALSE
#define COFFE_FALSE 0
#endif

#ifndef COFFE_LOGO
#define COFFE_LOGO \
        "   _____ ____  ______ ______ ______ \n" \
        "  / ____/ __ \\|  ____|  ____|  ____|\n" \
        " | |   | |  | | |__  | |__  | |__   \n" \
        " | |   | |  | |  __| |  __| |  __|  \n" \
        " | |___| |__| | |    | |    | |____ \n" \
        "  \\_____\\____/|_|    |_|    |______|\n"
#endif


#ifndef COFFE_MAX_INTSPACE
#define COFFE_MAX_INTSPACE 50000 // for 1D integration, the size of workspace
#endif

#ifndef COFFE_H0
#define COFFE_H0 (1./(2997.92458)) // H0 in units h/Mpc
#endif

/* useful macro for determining the size of an array */
#ifndef COFFE_ARRAY_SIZE
#define COFFE_ARRAY_SIZE(arr) ((size_t) (sizeof (arr) / sizeof (arr)[0]))
#endif



typedef struct coffe_class_struct_t
{
    void *file_content,
         *background,
         *thermodynamics,
         *perturb,
         *primordial,
         *nonlinear,
         *transfer,
         *spectra;
} coffe_class_struct_t;


int coffe_new_class_struct(
    coffe_class_struct_t *
);

int coffe_free_class_struct(
    coffe_class_struct_t *
);

typedef struct coffe_class_precision_parameters_t
{
    double k_per_decade_for_bao;
    double k_per_decade_for_pk;
    double start_large_k_at_tau_h_over_tau_k;
    int l_max_g;
    int l_max_ur;
    double tol_perturb_integration;
    double radiation_streaming_trigger_tau_over_tau_k;
    double ur_fluid_trigger_tau_over_tau_k;
} coffe_class_precision_parameters_t;


/**
    finds the max of an array of doubles
**/

double coffe_max_array_double(
    const double *array,
    const size_t size
);



/**
    contains the coefficients for the bias fitting at a given redshift
**/
typedef struct coffe_fit_coefficients_t
{
    double *coefficients;
    size_t size;
    /**
        I'm not sure how to parametrize the convergence; for now, I'm just
        reading the maximum degree and using that to fit, which outputs the chi squared;
        whether this is any good, I don't know
    **/
    double chisq;
    int degree;
    double z_min;
    double z_max;
} coffe_fit_coefficients_t;


typedef struct coffe_fit_coefficients_array_t
{
    coffe_fit_coefficients_t *array;
    size_t size;
} coffe_fit_coefficients_array_t;


int coffe_new_fit_coefficients(
    coffe_fit_coefficients_t *input
);


/**
    contains the coordinates for the correlation function at (z, r, mu)
**/

typedef struct coffe_corrfunc_coords_t
{
    double z_mean;
    double separation;
    double mu;
} coffe_corrfunc_coords_t;


/**
    same as above for multipoles
**/

typedef struct coffe_multipoles_coords_t
{
    double z_mean;
    double separation;
    int l;
} coffe_multipoles_coords_t;


/**
    -||- for RAMPs
**/

typedef struct coffe_average_multipoles_coords_t
{
    double z_mean;
    double separation;
    int l;
} coffe_average_multipoles_coords_t;


/**
    -||- for covariance of multipoles
**/

typedef struct coffe_covariance_coords_t
{
    int l1, l2;
    double separation1, separation2;
    double z_mean;
} coffe_covariance_coords_t;


/**
    simple wrapper with failsafe for malloc
**/

void *coffe_malloc(size_t size);


/**
    wrapper for 1D GSL interpolation structure
**/

typedef struct coffe_interpolation
{
    gsl_spline *spline;
    gsl_interp_accel *accel;
} coffe_interpolation;


/**
    wrapper for 2D GSL interpolation structure
**/

typedef struct coffe_interpolation2d
{
    gsl_spline2d *spline;
    gsl_interp_accel *xaccel;
    gsl_interp_accel *yaccel;
} coffe_interpolation2d;


/**
    simple enum to specify which correlation terms we need to compute
**/

enum coffe_integral_type
{
    NONINTEGRATED, SINGLE_INTEGRATED, DOUBLE_INTEGRATED
};


/**
    same as above, but for specifing the output type
**/

enum coffe_output_type
{
    CORRFUNC = 1,
    MULTIPOLES = 2,
    AVERAGE_MULTIPOLES = 3,
    COVARIANCE_MULTIPOLES = 4,
    COVARIANCE_AVERAGE_MULTIPOLES = 5
};


enum coffe_interp1d_type
{
    COFFE_INTERP_LINEAR = 1,
    COFFE_INTERP_POLYNOMIAL = 2,
    COFFE_INTERP_CSPLINE = 3,
    COFFE_INTERP_CSPLINE_PERIODIC = 4,
    COFFE_INTERP_AKIMA = 5,
    COFFE_INTERP_AKIMA_PERIODIC = 6,
    COFFE_INTERP_STEFFEN = 7
};


enum coffe_interp2d_type
{
    COFFE_INTERP2D_BILINEAR = 1,
    COFFE_INTERP2D_BICUBIC = 2
};


enum coffe_pk_type
{
    COFFE_PK_LINEAR = 0,
    COFFE_PK_LINEAR_CLASS = 1,
    COFFE_PK_NONLINEAR_HALOFIT = 2,
    COFFE_PK_NONLINEAR_HMCODE = 3
};


enum coffe_covariance_integration_method
{
    COFFE_COVARIANCE_INTEGRATION_GSL = 1,
    COFFE_COVARIANCE_INTEGRATION_FFTLOG = 2
};


enum coffe_integration_1d_type
{
    COFFE_INTEGRATION_GSL = 1,
    COFFE_INTEGRATION_DOUBLE_EXPONENTIAL = 2
};


/**
    contains all the values for n and l
    for which we need to compute the I^n_l;
    others are zero
**/

typedef struct nl_terms
{
    int n, l;
} nl_terms;


/**
    structure with all of the correlation contributions
**/

typedef struct coffe_correlation_contributions
{
    int den, rsd, len, d1, d2, g1, g2, g3, g4, g5;
} coffe_correlation_contributions;


/**
    contains all the parameters necessary to carry
    out the computation
**/

typedef struct coffe_parameters_t
{
    int output_type; /* correlation function (angular or full sky), or multipoles, or redshift averaged multipoles */

    double *mu; /* value for the angles */

    size_t mu_len; /* number of angles */

    coffe_correlation_contributions correlation_contrib; /* all of the correlation contributions (internally) */

    char **type_bg; /* background values to output */

    size_t type_bg_len; /* number of background values */

    int background_bins; /* number of bins for the background */

    int bessel_bins; /* number of bins for the bessel integrals */

    double inv_k_window; /* inverse of the exponential scale for the power spectrum */

    double Omega0_cdm; /* omega parameter for cold dark matter */

    double Omega0_m;

    int b_derivative, f_derivative, b_tilde_derivative, f_tilde_derivative;

    double Omega0_baryon; /* omega parameter for baryons */

    double Omega0_gamma; /* omega parameter for photons */

    double w0; /* parameter for the equation of state of dark energy */

    double wa; /* parameter for the equation of state of dark energy */

    double Omega0_de; /* present omega parameter for dark energy-like component */

    double *z_mean; /* mean redshift */

    size_t z_mean_len;

    double *deltaz;

    size_t deltaz_len;

    double *sep; /* all the separations */

    size_t sep_len; /* number of separations */

    int integration_method;

    int integration_bins;

    int nthreads; /* how many threads are used for the computation */

    char file_power_spectrum[COFFE_MAX_STRLEN]; /* file containing the PS */

    coffe_interpolation power_spectrum;

    coffe_interpolation power_spectrum_norm;

    coffe_interpolation2d power_spectrum2d;

    coffe_interpolation2d power_spectrum2d_norm;

    double k_min; /* min value to be taken from PS file */

    double k_max; /* max value to be taken from PS file */

    double k_min_norm; /* min value to be taken from PS file */

    double k_max_norm; /* max value to be taken from PS file */

    /* different biases */
    int read_galaxy_bias1, read_galaxy_bias2;

    char file_galaxy_bias1[COFFE_MAX_STRLEN], file_galaxy_bias2[COFFE_MAX_STRLEN];

    coffe_interpolation galaxy_bias1, galaxy_bias2, galaxy_bias3, galaxy_bias4;

    int read_magnification_bias1, read_magnification_bias2;

    int galaxy_bias_analytic;

    char file_magnification_bias1[COFFE_MAX_STRLEN], file_magnification_bias2[COFFE_MAX_STRLEN];

    coffe_interpolation magnification_bias1, magnification_bias2;

    int read_evolution_bias1, read_evolution_bias2;

    char file_evolution_bias1[COFFE_MAX_STRLEN], file_evolution_bias2[COFFE_MAX_STRLEN];

    coffe_interpolation evolution_bias1, evolution_bias2;

    int divergent; /* flag for divergent integrals */

    #ifndef COFFE_CYTHON
    config_t *conf; /* contains all the settings */
    #endif

    char timestamp[COFFE_MAX_STRLEN]; /* when the calculation started */

    char output_path[COFFE_MAX_STRLEN]; /* output directory for all the files */

    char output_prefix[COFFE_MAX_STRLEN]; /* output prefix for all the files */

    enum coffe_interp1d_type interp_method; /* method used for interpolation (linear, poly, etc.) */

    int *multipole_values; /* the multipoles to calculate */

    size_t multipole_values_len;

    int has_class;

    int has_cuba;

    double *zmin;

    size_t zmin_len;

    double *zmax;

    size_t zmax_len;

    double *fsky;

    size_t fsky_len;

    double *density1, *density2;

    size_t density1_len, density2_len;

    double *pixelsize;

    size_t pixelsize_len;

    /* the covariance is a 4-point function, hence we can in general
    specify 4 different populations */
    int covariance_pop1, covariance_pop2, covariance_pop3, covariance_pop4;

    /* which contributions should be considered for the covariance (CV-CV, CV-Poisson, and Poisson-Poisson */
    int covariance_cosmic, covariance_mixed, covariance_poisson;

    int covariance_integration_method;

    int covariance_integration_bins;

    enum coffe_interp2d_type covariance_interpolation_method;

    int covariance_window;

    int have_window;

    double window_size;

    /* for redshift averaged multipoles */
    double z_min, z_max;

    int flatsky_local;

    int flatsky_local_nonlocal;

    int flatsky_nonlocal;

    int verbose;

    /* stuff for CLASS only */

    int have_class; /* should CLASS be used or not? */

    double n_s; /* spectral index n_s */

    double sigma8; /* amplitude of primordial power spectrum sigma8 */

    double A_s; /* amplitude of primordial power spectrum A_s */

    double T_cmb, N_ur, m_ncdm, Omega0_nu; /* last one = neutrino energy density fraction */

    double YHe; /* helium fraction */

    int N_ncdm;

    double h; /* Damn you, little h! */

    double k_pivot; /* for CLASS */

    /* in C we can use void for anything that's a pointer */
    coffe_class_struct_t class_struct;

    coffe_class_precision_parameters_t class_precision;

    int flag;

    enum coffe_pk_type pk_type;

    int midpoint_approximation;

    int only_cross_correlations;

    enum coffe_integration_1d_type integration_1d_type;

    double integration_1d_prec;

    coffe_fit_coefficients_array_t galaxy_bias1_coefficients;
    coffe_fit_coefficients_array_t galaxy_bias2_coefficients;
    coffe_fit_coefficients_array_t magnification_bias1_coefficients;
    coffe_fit_coefficients_array_t magnification_bias2_coefficients;
    coffe_fit_coefficients_array_t evolution_bias1_coefficients;
    coffe_fit_coefficients_array_t evolution_bias2_coefficients;

    int degree_galaxy_bias1, degree_galaxy_bias2;
    int degree_magnification_bias1, degree_magnification_bias2;
    int degree_evolution_bias1, degree_evolution_bias2;

} coffe_parameters_t;

char *coffe_get_time(void);


/**
    auxiliary functions that may or may not be useful
**/

int read_1col(
    const char *filename,
    double **values,
    size_t *len
);

int read_2col(
    const char *filename,
    double **values1,
    double **values2,
    size_t *len
);

int write_1col(
    const char *filename,
    double *values,
    size_t len,
    const char *header
);

int write_2col(
    const char *filename,
    const double *values1,
    const double *values2,
    size_t len,
    const char *header,
    const char *sep
);

int write_ncol(
    size_t ncolumns,
    const char *filename,
    size_t len,
    const char *header,
    const char *sep,
    const double *values,
    ...
);

int write_ncol_null(
    const char *filename,
    size_t len,
    const char *header,
    const char *sep,
    const double *values,
    ...
);

int alloc_double_matrix(
    double ***values,
    size_t len1,
    size_t len2
);

int free_double_matrix(
    double ***values,
    size_t len
);

int write_matrix(
    const char *filename,
    const double **values,
    size_t len1,
    size_t len2,
    const char *header,
    const char *sep
);

int copy_matrix_array(
    double **destination,
    const double **source,
    size_t rows,
    size_t columns,
    size_t index,
    const char *type
);

int coffe_init_spline(
    coffe_interpolation *interp,
    const double *xi,
    const double *yi,
    const size_t bins,
    const enum coffe_interp1d_type interpolation_type
);

int coffe_init_spline2d(
    coffe_interpolation2d *interp,
    const double *xi,
    const double *yi,
    const double *zi,
    const size_t binsx,
    const size_t binsy,
    const enum coffe_interp2d_type interpolation_type
);

double coffe_interp_spline(
    const coffe_interpolation *interp,
    double value
);

double coffe_interp_spline2d(
    const coffe_interpolation2d *interp,
    const double value1,
    const double value2
);

int coffe_new_spline(
    coffe_interpolation *interp
);

int coffe_free_spline(
    coffe_interpolation *interp
);

int coffe_new_spline2d(
    coffe_interpolation2d *interp
);

int coffe_free_spline2d(
    coffe_interpolation2d *interp
);

int coffe_compare_ascending(
    const void *a,
    const void *b
);

int coffe_compare_descending(
    const void *a,
    const void *b
);

double coffe_dark_energy_eos(
    double z,
    const void *p
);

int coffe_parameters_free(
    coffe_parameters_t *par
);

double coffe_resolution_window(
    double x
);

int coffe_read_ncol(
    const char *filename,
    const size_t N,
    size_t *len,
    double **values,
    ...
);

double coffe_galaxy_bias(
    const double z
);


/**
    multiplies every item of a double `array` of `size` by a double `factor`
    MUST be alloc'd beforehand
**/
void coffe_rescale_array(
    double *array,
    const size_t size,
    const double factor
);


/**
    multiplies every item of a double `array` of `size` by a double `factor`
    MUST be alloc'd beforehand
**/
void coffe_rescale_array(
    double *array,
    const size_t size,
    const double factor
);


/**
    computes `array1` * `array2`^`power`, of same `size`,
    and stores the result into `array_out`
    MUST be alloc'd beforehand
**/
void coffe_multiply_power_array(
    double *array_out,
    const double *array1,
    const double *array2,
    const size_t size,
    const double power
);


/**
    integrates any 1D function `func` with arbitrary parameters `parameters`
    between `a` and `b`, with relative precision `prec`, and returns `result`,
    using either GSL (top) or double exponential quadrature (bottom)
**/

double coffe_integrate_1d_prec_gsl(
    double (*func)(
        double,
        void*
    ),
    const void *parameters,
    const double a,
    const double b,
    const double precision
);

double coffe_integrate_1d_prec_double_exponential(
    double (*func)(
        double,
        const void*
    ),
    const void *parameters,
    const double a,
    const double b,
    const double precision
);


double coffe_integrate_multidimensional(
#ifdef HAVE_CUBA
    int (*func)(
        const int *,
        const cubareal *,
        const int *,
        cubareal *,
        void *
    ),
#else
    double (*func)(
        double *,
        size_t,
        void *
    ),
#endif
    const void *parameters,
    const int integration_method,
    const int dims,
    const int integration_bins
);

int write_single_line(
    const char *filename,
    const double *values,
    ...
);

double *coffe_generate_range(
    const double xmin,
    const double xmax,
    const size_t steps
);

int coffe_approx_equal(
    const double a,
    const double b,
    const double rel_epsilon,
    const double abs_epsilon
);

int coffe_free_fit_coefficients(
    coffe_fit_coefficients_t *input
);

int coffe_new_fit_coefficients_array(
    coffe_fit_coefficients_array_t *input
);

int coffe_free_fit_coefficients_array(
    coffe_fit_coefficients_array_t *input
);


double coffe_legendre_integral(
    int n,
    int m,
    int a,
    int b
);


int coffe_kronecker_delta(
    int i,
    int j
);


int coffe_sign(int m);

void print_parameters(coffe_parameters_t *);

#endif
