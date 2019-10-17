#include <stdlib.h>
#include <gsl/gsl_sf_legendre.h>
#include "common.h"
#include "background.h"
#include "integrals.h"
#include "functions.h"
#include "signal.h"

#ifdef HAVE_DOUBLE_EXPONENTIAL
#include "tanhsinh.h"
#else
#include <gsl/gsl_integration.h>
#endif

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
#ifdef HAVE_DOUBLE_EXPONENTIAL
    const void *p
#else
    void *p
#endif
)
{
    struct coffe_integration_parameters_t *test =
        (struct coffe_integration_parameters_t *) p;
    struct coffe_background_t *bg = test->bg;
    struct coffe_parameters_t *par = test->par;
    struct coffe_integrals_t *integral = test->integral;
    double mu = test->mu;
    double sep = test->sep;
    return
        functions_single_integrated(
            par, bg, integral,
            par->z_mean, mu, sep, x
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
    struct coffe_integration_parameters_t *test =
        (struct coffe_integration_parameters_t *) p;
    struct coffe_background_t *bg = test->bg;
    struct coffe_parameters_t *par = test->par;
    struct coffe_integrals_t *integral = test->integral;
    double mu = test->mu;
    double sep = test->sep;
    double x1 = var[0], x2 = var[1];
#ifdef HAVE_CUBA
    value[0] =
#else
    return
#endif
        functions_double_integrated(
            par, bg, integral,
            par->z_mean, mu, sep, x1, x2
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
#ifdef HAVE_DOUBLE_EXPONENTIAL
    const void *p
#else
    void *p
#endif
)
{
    struct coffe_integration_parameters_t *all_params =
        (struct coffe_integration_parameters_t *) p;
    struct coffe_parameters_t *par = all_params->par;
    struct coffe_background_t *bg = all_params->bg;
    struct coffe_integrals_t *integral = all_params->integral;
    double sep = all_params->sep;
    int l = all_params->l;

    double mu = 2*x - 1;
    if (l == 0){
        return functions_nonintegrated(
            par, bg, integral,
            par->z_mean, mu, sep
        );
    }
    else{
        return functions_nonintegrated(
            par, bg, integral,
            par->z_mean, mu, sep
            )
           *gsl_sf_legendre_Pl(l, mu);
    }
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
    struct coffe_integration_parameters_t *params = (struct coffe_integration_parameters_t *) p;
    struct coffe_parameters_t *par = params->par;
    struct coffe_background_t *bg = params->bg;
    struct coffe_integrals_t *integral = params->integral;
    double sep = params->sep;

    double mu = 2*var[0] - 1, x = var[1];
#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_single_integrated(
            par, bg, integral, par->z_mean, mu, sep, x
        );
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_single_integrated(
            par, bg, integral, par->z_mean, mu, sep, x
        )
       *gsl_sf_legendre_Pl(params->l, mu);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_single_integrated(
            par, bg, integral, par->z_mean, mu, sep, x
        );
    }
    else{
        return functions_single_integrated(
            par, bg, integral, par->z_mean, mu, sep, x
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
    struct coffe_integration_parameters_t *params = (struct coffe_integration_parameters_t *) p;
    struct coffe_parameters_t *par = params->par;
    struct coffe_background_t *bg = params->bg;
    struct coffe_integrals_t *integral = params->integral;
    double sep = params->sep;

    double mu = 2*var[0] - 1, x1 = var[1], x2 = var[2];
#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_double_integrated(
            par, bg, integral, par->z_mean, mu, sep, x1, x2
        );
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_double_integrated(
            par, bg, integral, par->z_mean, mu, sep, x1, x2
        )
       *gsl_sf_legendre_Pl(params->l, mu);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_double_integrated(
            par, bg, integral, par->z_mean, mu, sep, x1, x2
        );
    }
    else{
        return functions_double_integrated(
            par, bg, integral, par->z_mean, mu, sep, x1, x2
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
    struct coffe_integration_parameters_t *params = (struct coffe_integration_parameters_t *) p;
    struct coffe_parameters_t *par = params->par;
    struct coffe_background_t *bg = params->bg;
    struct coffe_integrals_t *integral = params->integral;
    double sep = params->sep;

    double z1 =
        interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, par->z_min) + sep/2.
        );
    double z2 =
        interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, par->z_max) - sep/2.
        );

    double z = (z2 - z1)*var[0] + z1;
    double mu = 2*var[1] - 1;

#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_nonintegrated(
            par, bg, integral, z, mu, sep
        )
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_nonintegrated(
            par, bg, integral, z, mu, sep
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_nonintegrated(
            par, bg, integral, z, mu, sep
        )
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
    else{
        return functions_nonintegrated(
            par, bg, integral, z, mu, sep
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
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
    struct coffe_integration_parameters_t *params = (struct coffe_integration_parameters_t *) p;
    struct coffe_parameters_t *par = params->par;
    struct coffe_background_t *bg = params->bg;
    struct coffe_integrals_t *integral = params->integral;
    double sep = params->sep;

    double z1 =
        interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, par->z_min) + sep/2.
        );
    double z2 =
        interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, par->z_max) - sep/2.
        );

    double z = (z2 - z1)*var[0] + z1;
    double mu = 2*var[1] - 1;
    double x = var[2];

#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_single_integrated(
            par, bg, integral, z, mu, sep, x
        )
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_single_integrated(
            par, bg, integral, z, mu, sep, x
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_single_integrated(
            par, bg, integral, z, mu, sep, x
        )
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
    else{
        return functions_single_integrated(
            par, bg, integral, z, mu, sep, x
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
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
    struct coffe_integration_parameters_t *params = (struct coffe_integration_parameters_t *) p;
    struct coffe_parameters_t *par = params->par;
    struct coffe_background_t *bg = params->bg;
    struct coffe_integrals_t *integral = params->integral;
    double sep = params->sep;

    double z1 =
        interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, par->z_min) + sep/2.
        );
    double z2 =
        interp_spline(
            &bg->z_as_chi,
            interp_spline(&bg->comoving_distance, par->z_max) - sep/2.
        );

    double z = (z2 - z1)*var[0] + z1;
    double mu = 2*var[1] - 1;
    double x1 = var[2], x2 = var[3];

#ifdef HAVE_CUBA
    if (params->l == 0){
        value[0] = functions_double_integrated(
            par, bg, integral, z, mu, sep, x1, x2
        )
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
    else{
        value[0] = functions_double_integrated(
            par, bg, integral, z, mu, sep, x1, x2
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
        return EXIT_SUCCESS;
    }
#else
    if (params->l == 0){
        return functions_double_integrated(
            par, bg, integral, z, mu, sep, x1, x2
        )
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
    else{
        return functions_double_integrated(
            par, bg, integral, z, mu, sep, x1, x2
        )
       *gsl_sf_legendre_Pl(params->l, mu)
       /interp_spline(&bg->conformal_Hz, z)/(1 + z);
    }
#endif
}

#ifndef HAVE_CUBA
static int signal_integrate_gsl(
    gsl_monte_function integrand,
    int integration_method,
    int dims,
    int integration_bins,
    double *result
)
{
    gsl_rng *random;
    gsl_rng_env_setup();
    const gsl_rng_type *rng = gsl_rng_default;
    random = gsl_rng_alloc(rng);
    double lower[dims];
    double upper[dims];
    double error;
    for (int i = 0; i<dims; ++i){
        lower[i] = 0.0;
        upper[i] = 1.0;
    }
    switch (integration_method){
        case 0:{
            gsl_monte_plain_state *state =
                gsl_monte_plain_alloc(dims);
            gsl_monte_plain_integrate(
                &integrand, lower, upper,
                dims, integration_bins, random,
                state,
                result, &error
            );
            gsl_monte_plain_free(state);
            break;
        }
        case 1:{
            gsl_monte_miser_state *state =
                gsl_monte_miser_alloc(dims);
            gsl_monte_miser_integrate(
                &integrand, lower, upper,
                dims, integration_bins, random,
                state,
                result, &error
            );
            gsl_monte_miser_free(state);
            break;
        }
        case 2:{
            gsl_monte_vegas_state *state =
                gsl_monte_vegas_alloc(dims);
            gsl_monte_vegas_integrate(
                &integrand, lower, upper,
                dims, integration_bins, random,
                state,
                result, &error
            );
            gsl_monte_vegas_free(state);
            break;
        }
    }
    gsl_rng_free(random);

    return EXIT_SUCCESS;
}
#endif

/**
    Computes all of the integrals for the output.
    Sadly, C doesn't allow function overloading,
    so it's a bit clumsy
**/

double coffe_integrate(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t *integral,
    double sep,
    double mu,
    int l,
    enum coffe_integral_type flag_integral,
    enum coffe_output_type flag_output
)
{
    struct coffe_integration_parameters_t test = {
    .par = par,
    .bg = bg,
    .integral = integral,
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
                        par->z_mean, mu, sep
                    )/interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                }
                case MULTIPOLES:{
                    double result, error;
                    #ifdef HAVE_DOUBLE_EXPONENTIAL
                    result = tanhsinh_quad(
                        &multipoles_nonintegrated_integrand,
                        &test,
                        0., 1., 0.,
                        &error, NULL
                    );
                    #else
                    double prec = 1E-5;
                    gsl_function integrand;
                    integrand.function = &multipoles_nonintegrated_integrand;
                    integrand.params = &test;

                    gsl_integration_workspace *wspace =
                        gsl_integration_workspace_alloc(COFFE_MAX_INTSPACE);
                    gsl_integration_qag(
                        &integrand, 0., 1., 0,
                        prec, COFFE_MAX_INTSPACE,
                        GSL_INTEG_GAUSS61, wspace,
                        &result, &error
                    );
                    gsl_integration_workspace_free(wspace);
                    #endif
                    return (2*l + 1)*result
                        /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                }
                case AVERAGE_MULTIPOLES:{
                    const int dims = 2;

                    #ifdef HAVE_CUBA
                    int nregions, neval, fail;
                    double result[1], error[1], prob[1];
                    Cuhre(dims, 1,
                        average_multipoles_nonintegrated_integrand,
                        (void *)&test, 1,
                        1e-3, 0, 0,
                        1, par->integration_bins, 7,
                        NULL, NULL,
                        &nregions, &neval, &fail, result, error, prob
                    );
                    return (2*l + 1)*result[0]
                    /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                    #else
                    double result;
                    gsl_monte_function integrand;
                    integrand.f = &average_multipoles_nonintegrated_integrand;
                    integrand.dim = dims;
                    integrand.params = &test;
                    signal_integrate_gsl(
                        integrand,
                        par->integration_method,
                        dims,
                        par->integration_bins,
                        &result
                    );

                    return (2*l + 1)*result
                    /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                    #endif

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
                    double result, error;

                    #ifdef HAVE_DOUBLE_EXPONENTIAL
                    result = tanhsinh_quad(
                        &corrfunc_single_integrated_integrand,
                        &test,
                        0., 1., 0.,
                        &error, NULL
                    );
                    #else
                    double prec = 1E-5;
                    gsl_function integrand;
                    integrand.function = &corrfunc_single_integrated_integrand;
                    integrand.params = &test;

                    gsl_integration_workspace *wspace =
                        gsl_integration_workspace_alloc(COFFE_MAX_INTSPACE);
                    gsl_integration_qag(
                        &integrand, 0., 1., 0,
                        prec, COFFE_MAX_INTSPACE,
                        GSL_INTEG_GAUSS61, wspace,
                        &result, &error
                    );
                    gsl_integration_workspace_free(wspace);
                    #endif

                    return result/interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                }
                case MULTIPOLES:{
                    const int dims = 2;
                    #ifdef HAVE_CUBA
                    int nregions, neval, fail;
                    double result[1], error[1], prob[1];

                    Cuhre(dims, 1,
                        multipoles_single_integrated_integrand,
                        (void *)&test, 1,
                        5e-4, 0, 0,
                        1, par->integration_bins, 7,
                        NULL, NULL,
                        &nregions, &neval, &fail, result, error, prob
                    );
                    return (2*l + 1)*result[0]/interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                    #else
                    double result;
                    gsl_monte_function integrand;
                    integrand.dim = dims;
                    integrand.params = &test;
                    integrand.f = &multipoles_single_integrated_integrand;
                    signal_integrate_gsl(
                        integrand,
                        par->integration_method,
                        dims,
                        par->integration_bins,
                        &result
                    );

                    return (2*l + 1)*result
                        /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                    #endif
                }
                case AVERAGE_MULTIPOLES:{
                    const int dims = 3;
                    #ifdef HAVE_CUBA
                    int nregions, neval, fail;
                    double result[1], error[1], prob[1];
                    Cuhre(dims, 1,
                        average_multipoles_single_integrated_integrand,
                        (void *)&test, 1,
                        5e-4, 0, 0,
                        1, par->integration_bins, 7,
                        NULL, NULL,
                        &nregions, &neval, &fail, result, error, prob
                    );
                    return (2*l + 1)*result[0]
                    /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                    #else
                    double result;
                    gsl_monte_function integrand;
                    integrand.f = &average_multipoles_single_integrated_integrand;
                    integrand.dim = dims;
                    integrand.params = &test;
                    signal_integrate_gsl(
                        integrand,
                        par->integration_method,
                        dims,
                        par->integration_bins,
                        &result
                    );

                    return (2*l + 1)*result
                    /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                    #endif
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
                    const int dims = 2;
                    #ifdef HAVE_CUBA
                    int nregions, neval, fail;
                    double result[1], error[1], prob[1];

                    Cuhre(dims, 1,
                        corrfunc_double_integrated_integrand,
                        (void *)&test, 1,
                        5e-4, 0, 0,
                        1, par->integration_bins, 7,
                        NULL, NULL,
                        &nregions, &neval, &fail, result, error, prob
                    );
                    return result[0]/interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                    #else
                    double result;
                    gsl_monte_function integrand;
                    integrand.dim = dims;
                    integrand.params = &test;
                    integrand.f = &corrfunc_double_integrated_integrand;
                    signal_integrate_gsl(
                        integrand,
                        par->integration_method,
                        dims,
                        par->integration_bins,
                        &result
                    );

                    return result/interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                    #endif
                }
                case MULTIPOLES:{
                    const int dims = 3;
                    #ifdef HAVE_CUBA
                    int nregions, neval, fail;
                    double result[1], error[1], prob[1];

                    Cuhre(dims, 1,
                        multipoles_double_integrated_integrand,
                        (void *)&test, 1,
                        5e-4, 0, 0,
                        1, par->integration_bins, 7,
                        NULL, NULL,
                        &nregions, &neval, &fail, result, error, prob
                    );
                    return (2*l + 1)*result[0]/interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                    #else
                    double result;
                    gsl_monte_function integrand;
                    integrand.dim = dims;
                    integrand.params = &test;
                    integrand.f = &multipoles_double_integrated_integrand;
                    signal_integrate_gsl(
                        integrand,
                        par->integration_method,
                        dims,
                        par->integration_bins,
                        &result
                    );

                    return (2*l + 1)*result
                        /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                    #endif
                }
                case AVERAGE_MULTIPOLES:{
                    const int dims = 4;
                    #ifdef HAVE_CUBA
                    int nregions, neval, fail;
                    double result[1], error[1], prob[1];
                    Cuhre(dims, 1,
                        average_multipoles_double_integrated_integrand,
                        (void *)&test, 1,
                        5e-4, 0, 0,
                        1, par->integration_bins, 7,
                        NULL, NULL,
                        &nregions, &neval, &fail, result, error, prob
                    );
                    return (2*l + 1)*result[0]
                    /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                    #else
                    double result;
                    gsl_monte_function integrand;
                    integrand.f = &average_multipoles_double_integrated_integrand;
                    integrand.dim = dims;
                    integrand.params = &test;
                    signal_integrate_gsl(
                        integrand,
                        par->integration_method,
                        dims,
                        par->integration_bins,
                        &result
                    );

                    return (2*l + 1)*result
                    /interp_spline(&bg->D1, 0)/interp_spline(&bg->D1, 0);
                    #endif
                }
                default:
                    return 0.0;
            }
        }
        default:
            return 0.0;
    }
}


