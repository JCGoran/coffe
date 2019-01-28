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

#include <stdio.h>
#include "errors.h"

const char *error_msg(int code)
{
    switch (code) {
        COFFE_ERROR_CODES(COFFE_ERROR_TEXT)
    }

    return "Unknown error";
}

void print_error(int code)
{
    fprintf(stderr, "ERROR: %s\n", error_msg(code));
}

void print_error_verbose(int code, const char *message)
{
    fprintf(stderr,
        "ERROR: %s\n"
        "In file %s in function %s on line %d\n"
        "VALUE: %s\n",
        error_msg(code),
        __FILE__, __func__, __LINE__,
        message
    );
}
