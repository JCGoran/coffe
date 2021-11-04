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


typedef struct coffe_background_t
{
    coffe_interpolation z_as_chi; /* redshift as function of comoving distance */

    coffe_interpolation a; /* scale factor (normalized so that now a=1) */

    coffe_interpolation Hz; /* hubble parameter H(z) */

    coffe_interpolation conformal_Hz; /* conformal hubble parameter */

    coffe_interpolation conformal_Hz_prime; /* derivative of conformal hubble parameter wrt conformal time */

    coffe_interpolation D1; /* growth rate D_1(a) */

    coffe_interpolation f; /* growth function f=d(log D)/d(log a) */

    coffe_interpolation G1, G2;

    coffe_interpolation comoving_distance; /* comoving distance, dimensionless */

    int flag;

} coffe_background_t;


int coffe_background_init(
    const coffe_parameters_t *par,
    coffe_background_t *bg
);

int coffe_background_free(
    coffe_background_t *bg
);

/**
    checks whether `separation` fits inside a bin of half-with `deltaz`
    centered around a mean redshift `z_mean`, using cosmology `bg`
    Returns 1 if yes, 0 if no
**/
int coffe_check_range(
    const double separation,
    const double z_mean,
    const double deltaz,
    coffe_background_t *bg
);

#endif
