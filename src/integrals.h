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

/* for detection of half integers (spherical (int) vs. ordinary Bessel (half-int) functions) */
enum coffe_integer_state {COFFE_INTEGER, COFFE_HALF_INTEGER};

struct coffe_integral_t
{
    /* result I^n_l, or r^(n - l) I^n_l if n > l */
    struct coffe_interpolation result;
    /* result I^n_l(r, z), or r^(n - l) I^n_l if n > l */
    struct coffe_interpolation2d result2d;
    /* the suitable renormalization (if necessary) */
    struct coffe_interpolation2d renormalization;
    /* the suitable renormalization at zero separation (if necessary) */
    struct coffe_interpolation renormalization_zero_separation;
    /* the power and the degree of the spherical Bessel function, respectively */
    int n, l;
    /* if state = COFFE_integer_state, the passed value of n or l should be the NUMERATOR */
    enum coffe_integer_state state_n, state_l;
};


/* only has a pointer to multiple structures */
struct coffe_integral_array_t
{
    struct coffe_integral_t *integral;
    /* how many integrals are currently allocated */
    size_t size;
};


/* returns a pointer to the suitable array, or NULL if nothing is found */
struct coffe_integral_t *coffe_find_integral(
    const struct coffe_integral_array_t *array,
    const int n,
    const int l,
    const enum coffe_integer_state state_n,
    const enum coffe_integer_state state_l
);


int coffe_integrals_renormalizable(
    double **output_x,
    double **output_y,
    const size_t output_len,
    size_t *real_output_len,
    const struct coffe_interpolation *spectrum,
    const int n,
    const int l,
    const enum coffe_integer_state state_n,
    const enum coffe_integer_state state_l,
    const double x_min,
    const double x_max
);


/**
    computes all the nonzero I^n_l integrals
**/

int coffe_integrals_init(
    const struct coffe_parameters_t *par,
    const struct coffe_background_t *bg,
    struct coffe_integral_array_t *integral
);

int coffe_integrals_free(
    struct coffe_integral_array_t *integral
);

#endif
