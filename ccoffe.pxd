# Bunch of declarations from C to python. The idea here is to define only the
# quantities that will be used, for input, output or intermediate manipulation,
# by the python wrapper. For instance, in the precision structure, the only
# item used here is its error message. That is why nothing more is defined from
# this structure. The rest is internal in COFFE.
# If, for whatever reason, you need another, existing parameter from COFFE,
# remember to add it inside this cdef.

cdef extern from "<libconfig.h>":

    cdef struct config_t:
        pass


cdef extern from "<gsl/gsl_spline.h>":

    cdef struct gsl_spline:
        pass

    cdef struct gsl_interp_accel:
        pass


cdef extern from "<gsl/gsl_spline2d.h>":

    cdef struct gsl_spline2d:
        pass


DEF COFFE_MAX_STRLEN = 256


cdef extern from "error.h":
    print_error(int code)


cdef extern from "common.h":

    cdef struct coffe_corrfunc_coords_t:
        double z_mean
        double separation
        double mu

    cdef struct coffe_multipoles_coords_t:
        double z_mean
        double separation
        int l

    cdef struct coffe_corrfunc_t:
        coffe_corrfunc_coords_t coords
        double value

    cdef struct coffe_corrfunc_array_t:
        coffe_corrfunc_t *array
        size_t size

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
        CORRFUNC, MULTIPOLES, AVERAGE_MULTIPOLES

    cdef struct nl_terms:
        int n
        int l

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

        int output_type

        double *mu

        size_t mu_len

        coffe_correlation_contributions correlation_contrib

        nl_terms nonzero_terms[10]

        char **type_bg

        size_t type_bg_len

        int background_bins

        int bessel_bins

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

        char file_power_spectrum[COFFE_MAX_STRLEN]

        coffe_interpolation power_spectrum

        coffe_interpolation power_spectrum_norm

        double k_min

        double k_max

        double k_min_norm

        double k_max_norm

        int read_matter_bias1
        int read_matter_bias2

        char file_matter_bias1[COFFE_MAX_STRLEN]
        char file_matter_bias2[COFFE_MAX_STRLEN]

        coffe_interpolation matter_bias1
        coffe_interpolation matter_bias2

        int read_magnification_bias1
        int read_magnification_bias2

        int matter_bias_analytic

        char file_magnification_bias1[COFFE_MAX_STRLEN]
        char file_magnification_bias2[COFFE_MAX_STRLEN]

        coffe_interpolation magnification_bias1
        coffe_interpolation magnification_bias2

        int read_evolution_bias1
        int read_evolution_bias2

        char file_evolution_bias1[COFFE_MAX_STRLEN]
        char file_evolution_bias2[COFFE_MAX_STRLEN]

        coffe_interpolation evolution_bias1
        coffe_interpolation evolution_bias2

        int divergent

        config_t *conf

        char timestamp[COFFE_MAX_STRLEN]

        char output_path[COFFE_MAX_STRLEN]

        char output_prefix[COFFE_MAX_STRLEN]

        int interp_method

        int *multipole_values

        size_t multipole_values_len

        double *covariance_z_mean

        size_t covariance_z_mean_len

        double *covariance_deltaz

        size_t covariance_deltaz_len

        double *covariance_zmin

        size_t covariance_zmin_len

        double *covariance_zmax

        size_t covariance_zmax_len

        double *covariance_fsky

        size_t covariance_fsky_len

        double *covariance_density

        size_t covariance_density_len

        double covariance_pixelsize

        double covariance_minimum_separation

        int covariance_integration_method

        int covariance_integration_bins

        int covariance_interpolation_method

        int have_window

        double window_size

        double z_min
        double z_max

        int theta_len

        int flatsky

        int verbose

    #ifdef HAVE_CLASS

        int have_class

        double n_s

        double ln_10_pow_10_A_s

        double h

        double k_pivot

    #endif

    double *coffe_generate_range(
        double xmin,
        double xmax,
        size_t size
    )

    int coffe_parameters_free(
        coffe_parameters_t *
    )

cdef extern from "parser.h":

    int coffe_parse_default_parameters(
        coffe_parameters_t *
    )

    int coffe_parser_init(
        char *filename,
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


cdef extern from "corrfunc.h":
    cdef struct coffe_corrfunc_t:
        coffe_corrfunc_coords_t coords
        double value

    cdef struct coffe_corrfunc_array_t:
        coffe_corrfunc_t *array
        size_t size

cdef extern from "multipoles.h":
    cdef struct coffe_multipoles_t:
        coffe_multipoles_coords_t coords
        double value

    cdef struct coffe_multipoles_array_t:
        coffe_multipoles_t *array
        size_t size


#cdef from extern "multipoles.h":
#
#    cdef struct coffe_multipoles_t:
#        double **result
#        int *l
#        double *sep
#        size_t l_len, sep_len
#        int flag
#
#cdef from extern "average_multipoles.h":
#
#    cdef struct coffe_average_multipoles_t:
#        double **result
#        double *sep
#        size_t sep_len
#        int *l
#        size_t l_len
#        int flag
#
#cdef from extern "covariance.h":
#
#    cdef struct coffe_covariance_t:
#        double *z_mean, *deltaz, *density, *fsky
#        double *zmin, *zmax
#        size_t list_len
#        double pixelsize
#        double **sep
#        int *l
#        size_t *sep_len, l_len
#        double ***result
#        int flag
#
#    int coffe_covariance_init(
#        struct coffe_parameters_t *,
#        struct coffe_background_t *,
#        struct coffe_covariance_t *,
#        struct coffe_covariance_t *
#    )
#
#    int coffe_covariance_free(
#        struct coffe_covariance_t *
#    )
#
