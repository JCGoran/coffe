# Bunch of declarations from C to python. The idea here is to
# define only the quantities that will be used, for input, output
# or intermediate manipulation, by the python wrapper.  If, for
# whatever reason, you need another, existing parameter from COFFE,
# remember to add it inside this cdef.


cdef extern from "<gsl/gsl_errno.h>":

    cdef struct gsl_error_handler_t:
        pass

    cdef gsl_error_handler_t *gsl_set_error_handler_off()

    cdef void gsl_set_error_handler(gsl_error_handler_t *)


cdef extern from "<gsl/gsl_spline.h>":

    cdef struct gsl_spline:
        double *x
        double *y
        size_t size

    cdef struct gsl_interp_accel:
        pass


cdef extern from "<gsl/gsl_spline2d.h>":

    cdef struct gsl_spline2d:
        pass


DEF COFFE_MAX_STRLEN = 256


cdef extern from "common.h":

    cdef enum coffe_pk_type:
        COFFE_PK_LINEAR = 0,
        COFFE_PK_LINEAR_CLASS = 1,
        COFFE_PK_NONLINEAR_HALOFIT = 2,
        COFFE_PK_NONLINEAR_HMCODE = 3

    cdef enum coffe_interp1d_type:
        COFFE_INTERP_LINEAR = 1,
        COFFE_INTERP_POLYNOMIAL = 2,
        COFFE_INTERP_CSPLINE = 3,
        COFFE_INTERP_CSPLINE_PERIODIC = 4,
        COFFE_INTERP_AKIMA = 5,
        COFFE_INTERP_AKIMA_PERIODIC = 6,
        COFFE_INTERP_STEFFEN = 7

    cdef enum coffe_interp2d_type:
        COFFE_INTERP2D_BILINEAR = 1,
        COFFE_INTERP2D_BICUBIC = 2

    cdef enum coffe_covariance_integration_method:
        COFFE_COVARIANCE_INTEGRATION_STANDARD = 1,
        COFFE_COVARIANCE_INTEGRATION_FFTLOG = 2

    cdef struct coffe_corrfunc_coords_t:
        double z_mean
        double separation
        double mu

    cdef struct coffe_multipoles_coords_t:
        double z_mean
        double separation
        int l

    cdef struct coffe_covariance_coords_t:
        double separation1
        double separation2
        double z_mean
        int l1
        int l2

    cdef struct coffe_class_precision_parameters_t:
        double k_per_decade_for_bao
        double k_per_decade_for_pk
        double start_large_k_at_tau_h_over_tau_k
        int l_max_g
        int l_max_ur
        double tol_perturb_integration
        double radiation_streaming_trigger_tau_over_tau_k
        double ur_fluid_trigger_tau_over_tau_k

    cdef struct coffe_interpolation:
        gsl_spline *spline
        gsl_interp_accel *accel

    cdef struct coffe_interpolation2d:
        gsl_spline2d *spline
        gsl_interp_accel *xaccel
        gsl_interp_accel *yaccel

    cdef enum coffe_integral_type:
        NONINTEGRATED, SINGLE_INTEGRATED, DOUBLE_INTEGRATED

    cdef enum coffe_output_type:
        CORRFUNC = 1, MULTIPOLES = 2, AVERAGE_MULTIPOLES = 3, COVARIANCE_MULTIPOLES = 4

    cdef enum coffe_integration_1d_type:
        COFFE_INTEGRATION_GSL = 1,
        COFFE_INTEGRATION_DOUBLE_EXPONENTIAL = 2

    cdef struct coffe_correlation_contributions:
        int den
        int rsd
        int len
        int d1
        int d2
        int g1
        int g2
        int g3
        int g4
        int g5

    cdef struct coffe_parameters_t:

        coffe_class_precision_parameters_t class_precision

        int verbose

        int output_type

        double *mu

        size_t mu_len

        coffe_correlation_contributions correlation_contrib

        int background_bins

        int bessel_bins

        double Omega0_m

        double Omega0_cdm

        double Omega0_baryon

        double Omega0_gamma

        double w0

        double wa

        double Omega0_de

        double *z_mean

        size_t z_mean_len

        double *deltaz

        size_t deltaz_len

        double *sep

        size_t sep_len

        int integration_method

        int integration_bins

        int nthreads

        coffe_interpolation power_spectrum

        coffe_pk_type pk_type

        double k_min

        double k_max

        coffe_interpolation power_spectrum_norm

        double k_min_norm

        double k_max_norm

        coffe_interpolation2d power_spectrum2d
        coffe_interpolation2d power_spectrum2d_norm

        coffe_interpolation galaxy_bias1
        coffe_interpolation galaxy_bias2
        coffe_interpolation galaxy_bias3
        coffe_interpolation galaxy_bias4
        int degree_galaxy_bias1
        int degree_galaxy_bias2

        coffe_interpolation magnification_bias1
        coffe_interpolation magnification_bias2
        int degree_magnification_bias1
        int degree_magnification_bias2

        coffe_interpolation evolution_bias1
        coffe_interpolation evolution_bias2

        int divergent

        int covariance_window

        coffe_interp1d_type interp_method

        int *multipole_values

        size_t multipole_values_len

        double *zmin

        size_t zmin_len

        double *zmax

        size_t zmax_len

        double *fsky

        size_t fsky_len

        int has_class

        int has_cuba

        double *density1

        size_t density1_len

        double *density2

        size_t density2_len

        int covariance_pop1

        int covariance_pop2

        int covariance_pop3

        int covariance_pop4

        int covariance_cosmic

        int covariance_mixed

        int covariance_poisson

        double *pixelsize

        size_t pixelsize_len

        int covariance_integration_method

        int covariance_integration_bins

        int covariance_interpolation_method

        int have_window

        double window_size

        int flatsky_local

        int flatsky_local_nonlocal

        int flatsky_nonlocal

        double inv_k_window

        double n_s

        double sigma8

        double A_s

        double h

        double k_pivot

        double T_cmb

        double N_ur

        double m_ncdm

        int N_ncdm

        double Omega0_nu

        double YHe

        int only_cross_correlations

        coffe_integration_1d_type integration_1d_type

        double integration_1d_prec


    int coffe_parameters_free(
        coffe_parameters_t *
    )

    int coffe_init_spline(
        coffe_interpolation *interp,
        const double *xi,
        const double *yi,
        const size_t bins,
        const int interpolation_type
    )

    int coffe_init_spline2d(
        coffe_interpolation2d *interp,
        const double *xi,
        const double *yi,
        const double *zi,
        const size_t binsx,
        const size_t binsy,
        const int interpolation_type
    )

    double coffe_interp_spline(
        const coffe_interpolation *,
        const double
    )

    double coffe_interp_spline2d(
        const coffe_interpolation2d *,
        const double,
        const double
    )

    int coffe_free_spline(
        coffe_interpolation *
    )

    int coffe_free_spline2d(
        coffe_interpolation2d *
    )

    void *coffe_malloc(size_t)


cdef extern from "parser.h":

    int coffe_parse_default_parameters(
        coffe_parameters_t *
    )

    int parse_external_power_spectrum(
        coffe_parameters_t *
    )


cdef extern from "background.h":

    cdef struct coffe_background_t:

        coffe_interpolation z_as_chi
        coffe_interpolation a
        coffe_interpolation Hz
        coffe_interpolation conformal_Hz
        coffe_interpolation conformal_Hz_prime
        coffe_interpolation D1
        coffe_interpolation D1_prime
        coffe_interpolation f
        coffe_interpolation G1
        coffe_interpolation G2
        coffe_interpolation comoving_distance

        int flag

    int coffe_background_init(
        coffe_parameters_t *,
        coffe_background_t *
    )

    int coffe_background_free(
        coffe_background_t *
    )


cdef extern from "integrals.h":

    cdef enum coffe_integer_state:
        COFFE_INTEGER, COFFE_HALF_INTEGER

    cdef struct coffe_integral_t:
        coffe_interpolation result
        coffe_interpolation2d result2d
        coffe_interpolation2d renormalization
        coffe_interpolation renormalization_zero_separation
        int n
        int l
        coffe_integer_state state_n
        coffe_integer_state state_l

    cdef struct coffe_integral_array_t:
        coffe_integral_t *array
        size_t size

    int coffe_integrals_init(
        coffe_parameters_t *,
        coffe_background_t *,
        coffe_integral_array_t *
    )

    int coffe_integrals_free(
        coffe_integral_array_t *
    )

    coffe_integral_t *coffe_find_integral(
        const coffe_integral_array_t *array,
        const int n,
        const int l,
        const coffe_integer_state state_n,
        const coffe_integer_state state_l
    )




cdef extern from "signal.h":

    double coffe_integrate(
        coffe_parameters_t *par,
        coffe_background_t *bg,
        coffe_integral_array_t *integral,
        double z_mean,
        double sep,
        double mu,
        int l,
        coffe_integral_type flag_integral,
        coffe_output_type flag_output
    )


cdef extern from "corrfunc.h":

    cdef struct coffe_corrfunc_t:
        coffe_corrfunc_coords_t coords
        double value

    cdef struct coffe_corrfunc_array_t:
        coffe_corrfunc_t *array
        size_t size

    int coffe_corrfunc_init(
        coffe_parameters_t *,
        coffe_background_t *,
        coffe_integral_array_t *,
        coffe_corrfunc_array_t *
    )

    int coffe_corrfunc_free(
        coffe_corrfunc_array_t *
    )

cdef extern from "multipoles.h":
    cdef struct coffe_multipoles_t:
        coffe_multipoles_coords_t coords
        double value

    cdef struct coffe_multipoles_array_t:
        coffe_multipoles_t *array
        size_t size

    int coffe_multipoles_init(
        coffe_parameters_t *,
        coffe_background_t *,
        coffe_integral_array_t *,
        coffe_multipoles_array_t *
    )

    int coffe_multipoles_free(
        coffe_multipoles_array_t *
    )

cdef extern from "covariance.h":

    cdef struct coffe_covariance_t:
        coffe_covariance_coords_t coords
        double value

    cdef struct coffe_covariance_array_t:
        coffe_covariance_t *array
        size_t size

    int coffe_covariance_init(
        coffe_parameters_t *,
        coffe_background_t *,
        coffe_covariance_array_t *,
        coffe_covariance_array_t *
    )

    int coffe_covariance_free(
        coffe_covariance_array_t *
    )

    coffe_covariance_t coffe_covariance_find(
        const coffe_covariance_array_t *,
        const double z_mean,
        const int l1,
        const int l2,
        const double sep1,
        const double sep2
    )

#cdef from extern "average_multipoles.h":
#
#    cdef struct coffe_average_multipoles_t:
#        double **result
#        double *sep
#        size_t sep_len
#        int *l
#        size_t l_len
#        int flag
