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

#ifndef COFFE_INTEGRALS_H
#define COFFE_INTEGRALS_H

struct coffe_integrals_t
{
    int n, l;
    struct coffe_interpolation result;
    struct coffe_interpolation2d renormalization;
    struct coffe_interpolation renormalization0;
};

/**
    computes all the nonzero I^n_l integrals
**/

int coffe_integrals_init(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integrals_t integral[]
);

int coffe_integrals_free(
    struct coffe_integrals_t integral[]
);

#endif
