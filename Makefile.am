# The main product
bin_PROGRAMS = coffe # make all

#include .c and .h in SOURCES so that both appear in dist
coffe_SOURCES = \
    src/common.c \
    src/covariance.c \
    src/errors.c \
    src/parser.c \
    src/twofast.c \
    src/integrals.c \
    src/background.c \
    src/functions.c \
    src/tanhsinh.c \
    src/signal.c \
    src/corrfunc.c \
    src/multipoles.c \
    src/average_multipoles.c \
    src/output.c \
    src/main.c

TESTS = background_test integrals_test multipoles_test

check_PROGRAMS = background_test integrals_test multipoles_test

background_test_SOURCES = \
    tests/test_background.c \
    src/common.c \
    src/errors.c \
    src/parser.c \
    src/tanhsinh.c \
    src/background.c
background_test_CPPFLAGS = -Isrc/ -Itests/ -Itests/benchmarks/ -Wall -Wextra -Wpedantic -std=c99
background_test_LDFLAGS = $(LDFLAGS)

integrals_test_SOURCES = \
    tests/test_integrals.c \
    src/common.c \
    src/errors.c \
    src/parser.c \
    src/tanhsinh.c \
    src/background.c \
    src/twofast.c \
    src/integrals.c
integrals_test_CPPFLAGS = -Isrc/ -Itests/ -Itests/benchmarks/ -Wall -Wextra -Wpedantic -std=c99
integrals_test_LDFLAGS = $(LDFLAGS)

multipoles_test_SOURCES = \
    tests/test_multipoles.c \
    src/common.c \
    src/errors.c \
    src/parser.c \
    src/tanhsinh.c \
    src/background.c \
    src/twofast.c \
    src/integrals.c \
    src/functions.c \
    src/signal.c
multipoles_test_CPPFLAGS = -Isrc/ -Itests/ -Itests/benchmarks/ -Wall -Wextra -Wpedantic -std=c99
multipoles_test_LDFLAGS = $(LDFLAGS)


