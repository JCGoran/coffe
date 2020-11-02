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

struct coffe_corrfunc_ang_t
{
    double *result;
    double *theta;
    size_t theta_len;
    int flag;
};

struct coffe_corrfunc_t
{
    double **result;
    double *sep;
    double *mu;

    size_t sep_len;
    size_t mu_len;
    int flag;
};

struct coffe_corrfunc2d_t
{
    double **result;
    double *sep_parallel;
    double *sep_perpendicular;
    size_t sep_len;
    int flag;
};

int coffe_corrfunc_init(
    struct coffe_parameters_t *par,
    struct coffe_background_t *bg,
    struct coffe_integral_array_t *integral,
    struct coffe_corrfunc_ang_t *cf_ang,
    struct coffe_corrfunc_t *cf,
    struct coffe_corrfunc2d_t *cf2d
);

int coffe_corrfunc_ang_free(
    struct coffe_corrfunc_ang_t *cf_ang
);

int coffe_corrfunc_free(
    struct coffe_corrfunc_t *cf
);

int coffe_corrfunc2d_free(
    struct coffe_corrfunc2d_t *cf2d
);

#endif
