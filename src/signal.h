#ifndef COFFE_SIGNAL_H
#define COFFE_SIGNAL_H

#include "common.h"
#include "background.h"
#include "integrals.h"

struct coffe_integration_parameters_t
{
    struct coffe_parameters_t *par;
    struct coffe_background_t *bg;
    struct coffe_integral_array_t *integral;
    double z_mean;
    double z_min;
    double z_max;
    double sep;
    double mu;
    int l;
};

double coffe_integrate(
    coffe_parameters_t *par,
    coffe_background_t *bg,
    coffe_integral_array_t *integral,
    double z_mean,
    double z_min,
    double z_max,
    double sep,
    double mu,
    int l,
    enum coffe_integral_type flag_integral,
    enum coffe_output_type flag_output
);

#endif
