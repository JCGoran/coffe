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

#ifndef COFFE_FUNCTIONS_H
#define COFFE_FUNCTIONS_H

double functions_nonintegrated(
    const struct coffe_parameters_t *par,
    const struct coffe_background_t *bg,
    const struct coffe_integral_array_t *integral,
    const double z_mean,
    const double mu,
    const double r
);

double functions_single_integrated(
    const struct coffe_parameters_t *par,
    const struct coffe_background_t *bg,
    const struct coffe_integral_array_t *integral,
    const double z_mean,
    const double mu,
    const double r,
    const double x
);

double functions_double_integrated(
    const struct coffe_parameters_t *par,
    const struct coffe_background_t *bg,
    const struct coffe_integral_array_t *integral,
    const double z_mean,
    const double mu,
    const double r,
    const double x1,
    const double x2
);

double functions_flatsky_lensing_lensing_multipoles(
    const struct coffe_parameters_t *par,
    const struct coffe_background_t *bg,
    const struct coffe_integral_array_t *integral,
    const double z_mean,
    const double sep,
    const int l,
    const double x
);

#endif

