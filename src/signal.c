#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_legendre.h>
/* for the factorial function */
#include <gsl/gsl_sf_gamma.h>
#include "common.h"
#include "background.h"
#include "integrals.h"
#include "functions.h"
#include "signal.h"

#ifdef HAVE_CUBA
#include "cuba.h"
#else
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#endif


static double corrfunc_single_integrated_integrand(
    double x,
    void *p
)
{
    const struct coffe_integration_parameters_t *test =
        (const struct coffe_integration_parameters_t *) p;
    const coffe_background_t *bg = test->bg;
    const coffe_parameters_t *par = test->par;
    const coffe_integral_array_t *integral = test->integral;
    const double z_mean = test->z_mean;
    const double mu = test->mu;
    const double sep = test->sep;
    return
        functions_single_integrated(
            par, bg, integral,
            z_mean, mu, sep, x
        );
}


#ifdef HAVE_CUBA
static int corrfunc_double_integrated_integrand(
    const int *ndim, const cubareal var[],
    const int *ncomp, cubareal value[],
    void *p
)
#else
static double corrfunc_double_integrated_integrand(
    double *var, size_t dim, void *p
)
#endif
{
    const struct coffe_integration_parameters_t *test =
        (const struct coffe_integration_parameters_t *) p;
    const coffe_background_t *bg = test->bg;
    const coffe_parameters_t *par = test->par;
    const coffe_integral_array_t *integral = test->integral;
    const double z_mean = test->z_mean;
    const double mu = test->mu;
    const double sep = test->sep;
    const double x1 = var[0], x2 = var[1];
#ifdef HAVE_CUBA
    value[0] =
#else
    return
#endif
        functions_double_integrated(
            par, bg, integral,
            z_mean, mu, sep, x1, x2
        );
#ifdef HAVE_CUBA
    return EXIT_SUCCESS;
#endif
}


/**
    calculates all the nonintegrated terms
**/

static double multipoles_nonintegrated_integrand(
    double x,
    void *p
)
{
    const struct coffe_integration_parameters_t *all_params =
        (const struct coffe_integration_parameters_t *) p;
    const coffe_parameters_t *par = all_params->par;
    const coffe_background_t *bg = all_params->bg;
    const coffe_integral_array_t *integral = all_params->integral;
    const double z_mean = all_params->z_mean;
    const double sep = all_params->sep;
    const int l = all_params->l;

    const double mu = 2*x - 1;
    if (l == 0){
        return functions_nonintegrated(
            par, bg, integral,
            z_mean, mu, sep
        );
    }
    else{
        return functions_nonintegrated(
            par, bg, integral,
            z_mean, mu, sep
            )
           *gsl_sf_legendre_Pl(l, mu);
    }
}


static double multipoles_flatsky_integrand(
    double x,
    void *p
)
{
    const struct coffe_integration_parameters_t *all_params =
        (const struct coffe_integration_parameters_t *) p;
    const coffe_parameters_t *par = all_params->par;
    const coffe_background_t *bg = all_params->bg;
    const coffe_integral_array_t *integral = all_params->integral;
    const double z_mean = all_params->z_mean;
    const double sep = all_params->sep;
    const int l = all_params->l;

    return functions_flatsky_lensing_lensing_multipoles(
        par,
        bg,
        integral,
        z_mean,
        sep,
        l,
        x
    );
}


#ifdef HAVE_CUBA
static int multipoles_single_integrated_integrand(
    const int *ndim, const cubareal var[],
    const int *ncomp, cubareal value[],
    void *p
)
#else
static double multipoles_single_integrated_integrand(
    double *var, size_t dim, void *p
)
#endif
{
    const struct coffe_integration_parameters_t *params = (const struct coffe_integration_parameters_t *) p;
    const coffe_parameters_t *par = params->par;
    const coffe_background_t *bg = params->bg;
    const coffe_integral_array_t *integral = params->integral;
    const double sep = params->sep;
    const double z_mean = params->z_mean;

    const double mu = 2*var[0] - 1, x = var[1];
#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_single_integrated(
            par, bg, integral, z_mean, mu, sep, x
        );
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_single_integrated(
            par, bg, integral, z_mean, mu, sep, x
        )
       *gsl_sf_legendre_Pl(params->l, mu);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_single_integrated(
            par, bg, integral, z_mean, mu, sep, x
        );
    }
    else{
        return functions_single_integrated(
            par, bg, integral, z_mean, mu, sep, x
        )
       *gsl_sf_legendre_Pl(params->l, mu);
    }
#endif
}


#ifdef HAVE_CUBA
static int multipoles_double_integrated_integrand(
    const int *ndim, const cubareal var[],
    const int *ncomp, cubareal value[],
    void *p
)
#else
static double multipoles_double_integrated_integrand(
    double *var, size_t dim, void *p
)
#endif
{
    const struct coffe_integration_parameters_t *params = (struct coffe_integration_parameters_t *) p;
    const coffe_parameters_t *par = params->par;
    const coffe_background_t *bg = params->bg;
    const coffe_integral_array_t *integral = params->integral;
    const double sep = params->sep;
    const double z_mean = params->z_mean;

    const double mu = 2*var[0] - 1, x1 = var[1], x2 = var[2];
#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_double_integrated(
            par, bg, integral, z_mean, mu, sep, x1, x2
        );
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_double_integrated(
            par, bg, integral, z_mean, mu, sep, x1, x2
        )
       *gsl_sf_legendre_Pl(params->l, mu);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_double_integrated(
            par, bg, integral, z_mean, mu, sep, x1, x2
        );
    }
    else{
        return functions_double_integrated(
            par, bg, integral, z_mean, mu, sep, x1, x2
        )
       *gsl_sf_legendre_Pl(params->l, mu);
    }
#endif
}


/* integrand of nonintegrated terms for redshift averaged multipoles */

#ifdef HAVE_CUBA
static int average_multipoles_nonintegrated_integrand(
    const int *ndim, const cubareal var[],
    const int *ncomp, cubareal value[],
    void *p
)
#else
static double average_multipoles_nonintegrated_integrand(
    double *var, size_t dim, void *p
)
#endif
{
    const struct coffe_integration_parameters_t *params = (struct coffe_integration_parameters_t *) p;
    const coffe_parameters_t *par = params->par;
    const coffe_background_t *bg = params->bg;
    const coffe_integral_array_t *integral = params->integral;
    const double sep = params->sep;

    const double z1 =
        coffe_interp_spline(
            &bg->z_as_chi,
            coffe_interp_spline(&bg->comoving_distance, par->z_min) + sep/2.
        );
    const double z2 =
        coffe_interp_spline(
            &bg->z_as_chi,
            coffe_interp_spline(&bg->comoving_distance, par->z_max) - sep/2.
        );

    const double z = (z2 - z1)*var[0] + z1;
    const double mu = 2*var[1] - 1;

#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_nonintegrated(
            par, bg, integral, z, mu, sep
        )
       /coffe_interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_nonintegrated(
            par, bg, integral, z, mu, sep
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /coffe_interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_nonintegrated(
            par, bg, integral, z, mu, sep
        )
       /coffe_interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
    else{
        return functions_nonintegrated(
            par, bg, integral, z, mu, sep
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /coffe_interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
#endif
}




/* integrand of single integrated terms for redshift averaged multipoles */

#ifdef HAVE_CUBA
static int average_multipoles_single_integrated_integrand(
    const int *ndim, const cubareal var[],
    const int *ncomp, cubareal value[],
    void *p
)
#else
static double average_multipoles_single_integrated_integrand(
    double *var, size_t dim, void *p
)
#endif
{
    const struct coffe_integration_parameters_t *params = (struct coffe_integration_parameters_t *) p;
    const coffe_parameters_t *par = params->par;
    const coffe_background_t *bg = params->bg;
    const coffe_integral_array_t *integral = params->integral;
    const double sep = params->sep;

    const double z1 =
        coffe_interp_spline(
            &bg->z_as_chi,
            coffe_interp_spline(&bg->comoving_distance, par->z_min) + sep/2.
        );
    const double z2 =
        coffe_interp_spline(
            &bg->z_as_chi,
            coffe_interp_spline(&bg->comoving_distance, par->z_max) - sep/2.
        );

    const double z = (z2 - z1)*var[0] + z1;
    const double mu = 2*var[1] - 1;
    const double x = var[2];

#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_single_integrated(
            par, bg, integral, z, mu, sep, x
        )
       /coffe_interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_single_integrated(
            par, bg, integral, z, mu, sep, x
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /coffe_interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_single_integrated(
            par, bg, integral, z, mu, sep, x
        )
       /coffe_interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
    else{
        return functions_single_integrated(
            par, bg, integral, z, mu, sep, x
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /coffe_interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
#endif
}


/* integrand of double integrated terms for redshift averaged multipoles */

#ifdef HAVE_CUBA
static int average_multipoles_double_integrated_integrand(
    const int *ndim, const cubareal var[],
    const int *ncomp, cubareal value[],
    void *p
)
#else
static double average_multipoles_double_integrated_integrand(
    double *var, size_t dim, void *p
)
#endif
{
    const struct coffe_integration_parameters_t *params = (struct coffe_integration_parameters_t *) p;
    const coffe_parameters_t *par = params->par;
    const coffe_background_t *bg = params->bg;
    const coffe_integral_array_t *integral = params->integral;
    const double sep = params->sep;

    const double z1 =
        coffe_interp_spline(
            &bg->z_as_chi,
            coffe_interp_spline(&bg->comoving_distance, par->z_min) + sep/2.
        );
    const double z2 =
        coffe_interp_spline(
            &bg->z_as_chi,
            coffe_interp_spline(&bg->comoving_distance, par->z_max) - sep/2.
        );

    const double z = (z2 - z1)*var[0] + z1;
    const double mu = 2*var[1] - 1;
    const double x1 = var[2], x2 = var[3];

#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_double_integrated(
            par, bg, integral, z, mu, sep, x1, x2
        )
       /coffe_interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_double_integrated(
            par, bg, integral, z, mu, sep, x1, x2
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /coffe_interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_double_integrated(
            par, bg, integral, z, mu, sep, x1, x2
        )
       /coffe_interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
    else{
        return functions_double_integrated(
            par, bg, integral, z, mu, sep, x1, x2
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /coffe_interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
#endif
}


/**
    Computes all of the integrals for the output.
    Sadly, C doesn't allow function overloading,
    so it's a bit clumsy
**/

double coffe_integrate(
    coffe_parameters_t *par,
    coffe_background_t *bg,
    coffe_integral_array_t *integral,
    const double z_mean,
    const double sep,
    const double mu,
    const int l,
    enum coffe_integral_type flag_integral,
    enum coffe_output_type flag_output
)
{
    struct coffe_integration_parameters_t test = {
    .par = par,
    .bg = bg,
    .integral = integral,
    .z_mean = z_mean,
    .sep = sep,
    .mu = mu,
    .l = l
    };

    switch (flag_integral){
        case NONINTEGRATED:{
            /* first check if we really need it */
            int flag = COFFE_FALSE;
            if (
                par->correlation_contrib.den ||
                par->correlation_contrib.rsd ||
                par->correlation_contrib.d1 ||
                par->correlation_contrib.d2 ||
                par->correlation_contrib.g1 ||
                par->correlation_contrib.g2 ||
                par->correlation_contrib.g3
            ) flag = COFFE_TRUE;
            if (flag != COFFE_TRUE) return 0.0;

            /* if we do, check which one */
            switch (flag_output){
                case CORRFUNC:{
                    return functions_nonintegrated(
                        par, bg, integral,
                        z_mean, mu, sep
                    );
                }
                case MULTIPOLES:{
                    double result = 0;
                    if (par->integration_1d_type == COFFE_INTEGRATION_GSL){
                        result = coffe_integrate_1d_prec_gsl(
                            &multipoles_nonintegrated_integrand,
                            &test,
                            0,
                            1,
                            par->integration_1d_prec
                        );
                    }
                    else{
                        result = coffe_integrate_1d_prec_double_exponential(
                            &multipoles_nonintegrated_integrand,
                            &test,
                            0,
                            1,
                            par->integration_1d_prec
                        );
                    }
                    return (2 * l + 1) * result;
                }
                case AVERAGE_MULTIPOLES:{
                    return (2 * l + 1) * coffe_integrate_multidimensional(
                        &average_multipoles_nonintegrated_integrand,
                        (void *)&test,
                        par->integration_method,
                        2,
                        par->integration_bins
                    );
                }
                default:
                    return 0.0;
                }
        }
        case SINGLE_INTEGRATED:{
            int flag = COFFE_FALSE;
            if (
                (par->correlation_contrib.len && par->correlation_contrib.den) ||
                (par->correlation_contrib.len && par->correlation_contrib.rsd) ||
                (par->correlation_contrib.len && par->correlation_contrib.d1) ||
                (par->correlation_contrib.len && par->correlation_contrib.d2) ||
                (par->correlation_contrib.len && par->correlation_contrib.g1) ||
                (par->correlation_contrib.len && par->correlation_contrib.g2) ||
                (par->correlation_contrib.len && par->correlation_contrib.g3) ||
                (par->correlation_contrib.g4 && par->correlation_contrib.den) ||
                (par->correlation_contrib.g4 && par->correlation_contrib.rsd) ||
                (par->correlation_contrib.g4 && par->correlation_contrib.d1) ||
                (par->correlation_contrib.g4 && par->correlation_contrib.d2) ||
                (par->correlation_contrib.g4 && par->correlation_contrib.g1) ||
                (par->correlation_contrib.g4 && par->correlation_contrib.g2) ||
                (par->correlation_contrib.g4 && par->correlation_contrib.g3) ||
                (par->correlation_contrib.g5 && par->correlation_contrib.den) ||
                (par->correlation_contrib.g5 && par->correlation_contrib.rsd) ||
                (par->correlation_contrib.g5 && par->correlation_contrib.d1) ||
                (par->correlation_contrib.g5 && par->correlation_contrib.d2) ||
                (par->correlation_contrib.g5 && par->correlation_contrib.g1) ||
                (par->correlation_contrib.g5 && par->correlation_contrib.g2) ||
                (par->correlation_contrib.g5 && par->correlation_contrib.g3)
            ) flag = COFFE_TRUE;
            if (flag != COFFE_TRUE) return 0.0;

            switch(flag_output){
                case CORRFUNC:{
                    double result = 0;
                    if (par->integration_1d_type == COFFE_INTEGRATION_GSL){
                        result = coffe_integrate_1d_prec_gsl(
                            &corrfunc_single_integrated_integrand,
                            &test,
                            0,
                            1,
                            par->integration_1d_prec
                        );
                    }
                    else{
                        result = coffe_integrate_1d_prec_double_exponential(
                            &corrfunc_single_integrated_integrand,
                            &test,
                            0,
                            1,
                            par->integration_1d_prec
                        );
                    }

                    return result;
                }
                case MULTIPOLES:{
                    double final_result = 0;
                    /* flatsky density|RSD|d1-lensing are special */
                    if (
                        (
                            par->correlation_contrib.den ||
                            par->correlation_contrib.rsd ||
                            par->correlation_contrib.d1
                        ) &&
                        par->correlation_contrib.len &&
                        par->flatsky_local_nonlocal
                    ){
                        /* only density-lensing is non-zero in flat-sky */
                        if (par->correlation_contrib.den){
                            const double result = functions_flatsky_density_lensing_multipoles(
                                par,
                                bg,
                                integral,
                                z_mean,
                                sep,
                                l
                            );
                            final_result += 2 * M_PI * M_PI * result;
                        }
                    }
                    if (
                        (par->correlation_contrib.len && par->correlation_contrib.d2) ||
                        (par->correlation_contrib.len && par->correlation_contrib.g1) ||
                        (par->correlation_contrib.len && par->correlation_contrib.g2) ||
                        (par->correlation_contrib.len && par->correlation_contrib.g3) ||
                        (par->correlation_contrib.g4 && par->correlation_contrib.den) ||
                        (par->correlation_contrib.g4 && par->correlation_contrib.rsd) ||
                        (par->correlation_contrib.g4 && par->correlation_contrib.d1) ||
                        (par->correlation_contrib.g4 && par->correlation_contrib.d2) ||
                        (par->correlation_contrib.g4 && par->correlation_contrib.g1) ||
                        (par->correlation_contrib.g4 && par->correlation_contrib.g2) ||
                        (par->correlation_contrib.g4 && par->correlation_contrib.g3) ||
                        (par->correlation_contrib.g5 && par->correlation_contrib.den) ||
                        (par->correlation_contrib.g5 && par->correlation_contrib.rsd) ||
                        (par->correlation_contrib.g5 && par->correlation_contrib.d1) ||
                        (par->correlation_contrib.g5 && par->correlation_contrib.d2) ||
                        (par->correlation_contrib.g5 && par->correlation_contrib.g1) ||
                        (par->correlation_contrib.g5 && par->correlation_contrib.g2) ||
                        (par->correlation_contrib.g5 && par->correlation_contrib.g3) ||
                        (
                            (
                                par->correlation_contrib.den ||
                                par->correlation_contrib.rsd ||
                                par->correlation_contrib.d1
                            ) &&
                            par->correlation_contrib.len &&
                            !par->flatsky_local_nonlocal
                        )
                    ){

                        final_result += (2 * l + 1) * coffe_integrate_multidimensional(
                            &multipoles_single_integrated_integrand,
                            (void *)&test,
                            par->integration_method,
                            2,
                            par->integration_bins
                        );
                    }
                    return final_result;
                }
                case AVERAGE_MULTIPOLES:{
                    return (2 * l + 1) * coffe_integrate_multidimensional(
                        &average_multipoles_single_integrated_integrand,
                        (void *)&test,
                        par->integration_method,
                        3,
                        par->integration_bins
                    );
                }
                default:
                    return 0.0;
            }
        }
        case DOUBLE_INTEGRATED:{
            int flag = COFFE_FALSE;
            if (
                par->correlation_contrib.len ||
                par->correlation_contrib.g4 ||
                par->correlation_contrib.g5
            ) flag = COFFE_TRUE;
            if (flag != COFFE_TRUE) return 0.0;

            switch (flag_output){
                case CORRFUNC:{
                    return coffe_integrate_multidimensional(
                        &corrfunc_double_integrated_integrand,
                        (void *)&test,
                        par->integration_method,
                        2,
                        par->integration_bins
                    );
                }
                case MULTIPOLES:{
                    double final_result = 0;
                    /* lensing-lensing is special */
                    if (
                        par->correlation_contrib.len &&
                        par->flatsky_nonlocal &&
                        !par->only_cross_correlations &&
                        /* the odd ones are zero by construction, so we care only about the even */
                        l % 2 == 0
                    ){
                        double result = 0;
                        if (par->integration_1d_type == COFFE_INTEGRATION_GSL){
                            result = coffe_integrate_1d_prec_gsl(
                                &multipoles_flatsky_integrand,
                                &test,
                                0,
                                1,
                                par->integration_1d_prec
                            );
                        }
                        else{
                            result = coffe_integrate_1d_prec_double_exponential(
                                &multipoles_flatsky_integrand,
                                &test,
                                0,
                                1,
                                par->integration_1d_prec
                            );
                        }
                        /* the factor of 2 pi^2 is because of the definition of the I^n_l integrals */
                        final_result += 2 * M_PI * M_PI * result;
                    }
                    /* this part can run now that we have lensing-lensing multipoles */
                    if (
                        par->correlation_contrib.g4 ||
                        par->correlation_contrib.g5 ||
                        (
                            par->correlation_contrib.len &&
                            !par->flatsky_nonlocal
                        )
                    ){
                        final_result += (2 * l + 1) * coffe_integrate_multidimensional(
                            &multipoles_double_integrated_integrand,
                            (void *)&test,
                            par->integration_method,
                            3,
                            par->integration_bins
                        );
                    }
                    return final_result;
                }
                case AVERAGE_MULTIPOLES:{
                    return (2 * l + 1) * coffe_integrate_multidimensional(
                        &average_multipoles_double_integrated_integrand,
                        (void *)&test,
                        par->integration_method,
                        4,
                        par->integration_bins
                    );
                }
                default:
                    return 0.0;
            }
        }
        default:
            return 0.0;
    }
}
