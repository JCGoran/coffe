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

#ifndef COFFE_FUNCTIONS_H
#define COFFE_FUNCTIONS_H

/**
    contains all the parameters needed for evaluation of non-integrated terms
**/
struct functions_params
{
    struct coffe_parameters_t *par;
    struct coffe_background_t *bg;
    struct coffe_integrals_t *integral;
    double sep, mu;
    int l;
};

double functions_nonintegrated(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t integral[],
    double z_mean,
    double mu,
    double r
);

double functions_single_integrated(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t integral[],
    double z_mean,
    double mu,
    double r,
    double x
);

double functions_double_integrated(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t integral[],
    double z_mean,
    double mu,
    double r,
    double x1,
    double x2
);

#endif

