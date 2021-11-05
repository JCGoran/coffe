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

#ifndef COFFE_MULTIPOLES_H
#define COFFE_MULTIPOLES_H

#include "common.h"

/**
    contains the value of the multipole for one point, i.e. at (l, r, z)
**/
typedef struct coffe_multipoles_t
{
    coffe_multipoles_coords_t coords;
    double value;
} coffe_multipoles_t;


/**
    container for the above structure with a size
**/
typedef struct coffe_multipoles_array_t
{
    coffe_multipoles_t *array;
    size_t size;
} coffe_multipoles_array_t;


/**
    main function for computing the multipoles
**/
int coffe_multipoles_init(
    coffe_parameters_t *par,
    coffe_background_t *bg,
    coffe_integral_array_t *integral,
    coffe_multipoles_array_t *mp
);

int coffe_multipoles_free(
    coffe_multipoles_array_t *mp
);

#endif
