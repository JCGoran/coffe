/*
 * This file is part of COFFE
 * Copyright (C) 2018 Goran Jelic-Cizmek
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

#ifndef COFFE_ERRORS_H
#define COFFE_ERRORS_H

#ifndef COFFE_ERROR_CODES
#define COFFE_ERROR_CODES(X)                                            \
    X(0,   PROG_SUCCESS,            "No error at all!")         \
    X(10,   PROG_FAIL,               "Unknown error occurred")   \
    X(20,   PROG_NULL_PARAM,         "Illegal null parameter")   \
    X(30,   PROG_ALLOC_ERROR,        "Allocation failed")        \
    X(40,   PROG_READ_ERROR,         "Cannot read file")         \
    X(45,   PROG_POS_ERROR,          "Cannot rewind file")       \
    X(50,   PROG_WRITE_ERROR,        "Cannot write to file")     \
    X(60,   PROG_OPEN_ERROR,         "Cannot open file")         \
    X(70,   PROG_CLOSE_ERROR,        "Cannot close file")        \
    X(80,   PROG_PARSE_ERROR,        "Cannot parse setting")     \
    X(90,   PROG_VALUE_ERROR,        "Cannot parse value")       \
    X(100,  PROG_MKDIR_ERROR,        "A file with the same name already exists")

#endif
#define ERROR_ENUM(ID, NAME, TEXT) NAME = ID,
#define COFFE_ERROR_TEXT(ID, NAME, TEXT) case ID: return TEXT;

enum {
    COFFE_ERROR_CODES(ERROR_ENUM)
};

const char *error_msg(int code);

void print_error(int code);

void print_error_verbose(int code, const char *message);

#endif
