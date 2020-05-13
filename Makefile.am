# The main product
bin_PROGRAMS = coffe # make all

#include .c and .h in SOURCES so that both appear in dist
coffe_SOURCES = \
    src/common.c \
    src/covariance.c \
    src/utils.c \
    src/twobessel.c \
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

TESTS = test_background test_integrals test_multipoles test_covariance

check_PROGRAMS = test_background test_integrals test_multipoles test_covariance

test_background_SOURCES = \
    tests/test_background.c \
    src/common.c \
    src/errors.c \
    src/parser.c \
    src/tanhsinh.c \
    src/background.c
test_background_CPPFLAGS = -Isrc/ -Itests/ -Itests/benchmarks/ -Wall -Wextra -Wpedantic
test_background_LDFLAGS = $(LDFLAGS)
test_background_AM_CFLAGS = $(AM_CFLAGS)

test_integrals_SOURCES = \
    tests/test_integrals.c \
    src/common.c \
    src/errors.c \
    src/parser.c \
    src/tanhsinh.c \
    src/background.c \
    src/twofast.c \
    src/integrals.c
test_integrals_CPPFLAGS = -Isrc/ -Itests/ -Itests/benchmarks/ -Wall -Wextra -Wpedantic
test_integrals_LDFLAGS = $(LDFLAGS)
test_integrals_AM_CFLAGS = $(AM_CFLAGS)

test_multipoles_SOURCES = \
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
test_multipoles_CPPFLAGS = -Isrc/ -Itests/ -Itests/benchmarks/ -Wall -Wextra -Wpedantic
test_multipoles_LDFLAGS = $(LDFLAGS)
test_multipoles_AM_CFLAGS = $(AM_CFLAGS)

test_covariance_SOURCES = \
    tests/test_covariance.c \
    src/common.c \
    src/errors.c \
    src/parser.c \
    src/tanhsinh.c \
    src/background.c \
    src/twofast.c \
    src/integrals.c \
    src/functions.c \
    src/twobessel.c \
    src/utils.c \
    src/covariance.c
test_covariance_CPPFLAGS = -Isrc/ -Itests/ -Itests/benchmarks/ -Wall -Wextra -Wpedantic
test_covariance_LDFLAGS = $(LDFLAGS)
test_covariance_AM_CFLAGS = $(AM_CFLAGS)

# this is an expected failure since l=2,4 multipoles
# are extremely small for den, d1, d2, g1, g2, g3, g4, g5,
# and the integration has a lot of numerical cancellations,
# which causes rounding errors;
# this by no means reflects the general accuracy of COFFE,
# just a limitation of floating point operations
XFAIL_TESTS=test_multipoles
