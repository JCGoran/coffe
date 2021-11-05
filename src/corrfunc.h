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

#include "common.h"

typedef struct coffe_corrfunc_t
{
    coffe_corrfunc_coords_t coords;
    double value;
} coffe_corrfunc_t;

typedef struct coffe_corrfunc_array_t
{
    coffe_corrfunc_t *array;
    size_t size;
} coffe_corrfunc_array_t;

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
