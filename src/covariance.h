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

#ifndef COFFE_COVARIANCE
#define COFFE_COVARIANCE

#include "background.h"
#include "common.h"

typedef struct coffe_covariance_t
{
    coffe_covariance_coords_t coords;
    double value;
} coffe_covariance_t;

typedef struct coffe_average_covariance_t
{
    coffe_average_covariance_coords_t coords;
    double value;
} coffe_average_covariance_t;

typedef struct coffe_covariance_array_t
{
    coffe_covariance_t *array;
    size_t size;
} coffe_covariance_array_t;

typedef struct coffe_average_covariance_array_t
{
    coffe_average_covariance_t *array;
    size_t size;
} coffe_average_covariance_array_t;

int coffe_covariance_init(
    coffe_parameters_t *par,
    coffe_background_t *bg,
    coffe_covariance_array_t *cov_mp,
    coffe_average_covariance_array_t *cov_ramp
);

int coffe_covariance_free(
    coffe_covariance_array_t *cov
);

int coffe_average_covariance_free(
    coffe_average_covariance_array_t *cov
);

coffe_covariance_t coffe_covariance_find(
    const coffe_covariance_array_t *cov,
    const double z_mean,
    const int l1,
    const int l2,
    const double sep1,
    const double sep2
);

double covariance_coefficient(
    double bias1,
    double bias2,
    double growth_rate,
    int ell
);

#endif
