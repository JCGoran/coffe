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

#ifndef COFFE_CORRFUNC_H
#define COFFE_CORRFUNC_H

enum coffe_corrfunc_coordinates_enum
{
    COFFE_COORDINATE_MEAN_REDSHIFT,
    COFFE_COORDINATE_SEPARATION,
    COFFE_COORDINATE_ANGLE_ALPHA,
    COFFE_COORDINATE_ANGLE_BETA,
    COFFE_COORDINATE_ANGLE_GAMMA,
    COFFE_COORDINATE_ANGLE_MU,
    COFFE_COORDINATE_SEPARATION_PARALLEL,
    COFFE_COORDINATE_SEPARATION_PERP
};

typedef struct coffe_corrfunc_coordinate_t
{
    double value;
    enum coffe_corrfunc_coordinates_enum name;
} coffe_corrfunc_coordinate_t;

/* simple wrapper for the above */
typedef struct coffe_corrfunc_coordinate_array_t
{
    coffe_corrfunc_coordinate_t value[3];
} coffe_corrfunc_coordinate_array_t;

typedef struct coffe_corrfunc_t
{
    double value;
    coffe_corrfunc_coordinate_array_t coordinates;
} coffe_corrfunc_t;

typedef struct coffe_corrfunc_array_t
{
    coffe_corrfunc_t *value;
    size_t size;
} coffe_corrfunc_array_t;

/* for switching between coordinates */
int coffe_corrfunc_coordinate_transform(
    const coffe_corrfunc_coordinate_array_t input,
    const enum coffe_corrfunc_coordinates_enum coord1,
    const enum coffe_corrfunc_coordinates_enum coord2,
    const enum coffe_corrfunc_coordinates_enum coord3,
    coffe_corrfunc_coordinate_array_t *output
);

int coffe_corrfunc_init(
    coffe_parameters_t *par,
    coffe_background_t *bg,
    coffe_integral_array_t *integral,
    coffe_corrfunc_array_t *cf
);

int coffe_corrfunc_free(
    coffe_corrfunc_array_t *cf
);

#endif
