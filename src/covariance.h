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

struct coffe_covariance_t
{
    double *z_mean, *deltaz, *density, *fsky;
    double *zmin, *zmax; /* placeholder for RAMPs */
    size_t list_len; /* length of z_mean, density and fsky */
    double *pixelsize;
    double **sep;
    int *l;
    size_t *sep_len, l_len; /* sep_len is actually number of pixels */
    double ***result; /* result[i][j*l_len + k][l*sep_len[i] + m] gives the covariance at z_mean[i], with f_sky[i], density[i], of multipoles l[j] and l[k], at separations sep[i][l] and sep[i][m] */
    int flag;
};

int coffe_covariance_init(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_covariance_t *cov_mp,
    struct coffe_covariance_t *cov_ramp
);

int coffe_covariance_free(
    struct coffe_covariance_t *cov
);

#endif
