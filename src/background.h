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

#ifndef COFFE_BACKGROUND_H
#define COFFE_BACKGROUND_H


struct coffe_background_t
{
    struct coffe_interpolation z_as_chi; /* redshift as function of comoving distance */

    struct coffe_interpolation a; /* scale factor (normalized so that now a=1) */

    struct coffe_interpolation Hz; /* hubble parameter H(z) */

    struct coffe_interpolation conformal_Hz; /* conformal hubble parameter */

    struct coffe_interpolation conformal_Hz_prime; /* derivative of conformal hubble parameter wrt conformal time */

    struct coffe_interpolation D1; /* growth rate D_1(a) */

    struct coffe_interpolation f; /* growth function f=d(log D)/d(log a) */

    struct coffe_interpolation G1, G2;

    struct coffe_interpolation comoving_distance; /* comoving distance, dimensionless */

    int flag;

};


int coffe_background_init(
    const struct coffe_parameters_t *par,
    struct coffe_background_t *bg
);

int coffe_background_free(
    struct coffe_background_t *bg
);

#endif
