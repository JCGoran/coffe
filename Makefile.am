# The main product
bin_PROGRAMS = coffe # make all

lib_LIBRARIES = libcoffe.a

libcoffe_adir = $(includedir)/coffe

# the non-source files that will be included
# when running `make dist`
dist_pkgdata_DATA = \
    WAVENUMBER_HEADER.dat \
    POWER_SPECTRUM_HEADER.dat \
    settings.cfg \
    separations.dat \
    PkL_CLASS.dat

EXTRA_DIST = \
    tests/benchmarks/benchmark_all_corrfunc0.dat \
    tests/benchmarks/benchmark_all_corrfunc1.dat \
    tests/benchmarks/benchmark_all_corrfunc2.dat \
    tests/benchmarks/benchmark_all_corrfunc3.dat \
    tests/benchmarks/benchmark_all_corrfunc4.dat \
    tests/benchmarks/benchmark_all_multipoles0.dat \
    tests/benchmarks/benchmark_all_multipoles2.dat \
    tests/benchmarks/benchmark_all_multipoles4.dat \
    tests/benchmarks/benchmark_background.dat \
    tests/benchmarks/benchmark_d1_corrfunc0.dat \
    tests/benchmarks/benchmark_d1_corrfunc1.dat \
    tests/benchmarks/benchmark_d1_corrfunc2.dat \
    tests/benchmarks/benchmark_d1_corrfunc3.dat \
    tests/benchmarks/benchmark_d1_corrfunc4.dat \
    tests/benchmarks/benchmark_d1_multipoles0.dat \
    tests/benchmarks/benchmark_d1_multipoles2.dat \
    tests/benchmarks/benchmark_d1_multipoles4.dat \
    tests/benchmarks/benchmark_d2_corrfunc0.dat \
    tests/benchmarks/benchmark_d2_corrfunc1.dat \
    tests/benchmarks/benchmark_d2_corrfunc2.dat \
    tests/benchmarks/benchmark_d2_corrfunc3.dat \
    tests/benchmarks/benchmark_d2_corrfunc4.dat \
    tests/benchmarks/benchmark_d2_multipoles0.dat \
    tests/benchmarks/benchmark_d2_multipoles2.dat \
    tests/benchmarks/benchmark_d2_multipoles4.dat \
    tests/benchmarks/benchmark_den_corrfunc0.dat \
    tests/benchmarks/benchmark_den_corrfunc1.dat \
    tests/benchmarks/benchmark_den_corrfunc2.dat \
    tests/benchmarks/benchmark_den_corrfunc3.dat \
    tests/benchmarks/benchmark_den_corrfunc4.dat \
    tests/benchmarks/benchmark_den_multipoles0.dat \
    tests/benchmarks/benchmark_den_multipoles2.dat \
    tests/benchmarks/benchmark_den_multipoles4.dat \
    tests/benchmarks/benchmark_flatsky_density_lensing_multipoles0.dat \
    tests/benchmarks/benchmark_flatsky_density_lensing_multipoles2.dat \
    tests/benchmarks/benchmark_flatsky_density_lensing_multipoles4.dat \
    tests/benchmarks/benchmark_flatsky_lensing_lensing_corrfunc0.dat \
    tests/benchmarks/benchmark_flatsky_lensing_lensing_corrfunc1.dat \
    tests/benchmarks/benchmark_flatsky_lensing_lensing_corrfunc2.dat \
    tests/benchmarks/benchmark_flatsky_lensing_lensing_corrfunc3.dat \
    tests/benchmarks/benchmark_flatsky_lensing_lensing_corrfunc4.dat \
    tests/benchmarks/benchmark_flatsky_lensing_lensing_multipoles0.dat \
    tests/benchmarks/benchmark_flatsky_lensing_lensing_multipoles2.dat \
    tests/benchmarks/benchmark_flatsky_lensing_lensing_multipoles4.dat \
    tests/benchmarks/benchmark_g1_corrfunc0.dat \
    tests/benchmarks/benchmark_g1_corrfunc1.dat \
    tests/benchmarks/benchmark_g1_corrfunc2.dat \
    tests/benchmarks/benchmark_g1_corrfunc3.dat \
    tests/benchmarks/benchmark_g1_corrfunc4.dat \
    tests/benchmarks/benchmark_g1_multipoles0.dat \
    tests/benchmarks/benchmark_g1_multipoles2.dat \
    tests/benchmarks/benchmark_g1_multipoles4.dat \
    tests/benchmarks/benchmark_g2_corrfunc0.dat \
    tests/benchmarks/benchmark_g2_corrfunc1.dat \
    tests/benchmarks/benchmark_g2_corrfunc2.dat \
    tests/benchmarks/benchmark_g2_corrfunc3.dat \
    tests/benchmarks/benchmark_g2_corrfunc4.dat \
    tests/benchmarks/benchmark_g2_multipoles0.dat \
    tests/benchmarks/benchmark_g2_multipoles2.dat \
    tests/benchmarks/benchmark_g2_multipoles4.dat \
    tests/benchmarks/benchmark_g3_corrfunc0.dat \
    tests/benchmarks/benchmark_g3_corrfunc1.dat \
    tests/benchmarks/benchmark_g3_corrfunc2.dat \
    tests/benchmarks/benchmark_g3_corrfunc3.dat \
    tests/benchmarks/benchmark_g3_corrfunc4.dat \
    tests/benchmarks/benchmark_g3_multipoles0.dat \
    tests/benchmarks/benchmark_g3_multipoles2.dat \
    tests/benchmarks/benchmark_g3_multipoles4.dat \
    tests/benchmarks/benchmark_g4_corrfunc0.dat \
    tests/benchmarks/benchmark_g4_corrfunc1.dat \
    tests/benchmarks/benchmark_g4_corrfunc2.dat \
    tests/benchmarks/benchmark_g4_corrfunc3.dat \
    tests/benchmarks/benchmark_g4_corrfunc4.dat \
    tests/benchmarks/benchmark_g4_multipoles0.dat \
    tests/benchmarks/benchmark_g4_multipoles2.dat \
    tests/benchmarks/benchmark_g4_multipoles4.dat \
    tests/benchmarks/benchmark_g5_corrfunc0.dat \
    tests/benchmarks/benchmark_g5_corrfunc1.dat \
    tests/benchmarks/benchmark_g5_corrfunc2.dat \
    tests/benchmarks/benchmark_g5_corrfunc3.dat \
    tests/benchmarks/benchmark_g5_corrfunc4.dat \
    tests/benchmarks/benchmark_g5_multipoles0.dat \
    tests/benchmarks/benchmark_g5_multipoles2.dat \
    tests/benchmarks/benchmark_g5_multipoles4.dat \
    tests/benchmarks/benchmark_integral0.dat \
    tests/benchmarks/benchmark_integral1.dat \
    tests/benchmarks/benchmark_integral2.dat \
    tests/benchmarks/benchmark_integral3.dat \
    tests/benchmarks/benchmark_integral4.dat \
    tests/benchmarks/benchmark_integral5.dat \
    tests/benchmarks/benchmark_integral6.dat \
    tests/benchmarks/benchmark_integral7.dat \
    tests/benchmarks/benchmark_integral8.dat \
    tests/benchmarks/benchmark_integral8_renormalization.dat \
    tests/benchmarks/benchmark_integral9.dat \
    tests/benchmarks/benchmark_len_corrfunc0.dat \
    tests/benchmarks/benchmark_len_corrfunc1.dat \
    tests/benchmarks/benchmark_len_corrfunc2.dat \
    tests/benchmarks/benchmark_len_corrfunc3.dat \
    tests/benchmarks/benchmark_len_corrfunc4.dat \
    tests/benchmarks/benchmark_len_multipoles0.dat \
    tests/benchmarks/benchmark_len_multipoles2.dat \
    tests/benchmarks/benchmark_len_multipoles4.dat \
    tests/benchmarks/benchmark_multipoles_covariance_00.dat \
    tests/benchmarks/benchmark_multipoles_covariance_02.dat \
    tests/benchmarks/benchmark_multipoles_covariance_04.dat \
    tests/benchmarks/benchmark_multipoles_covariance_20.dat \
    tests/benchmarks/benchmark_multipoles_covariance_22.dat \
    tests/benchmarks/benchmark_multipoles_covariance_24.dat \
    tests/benchmarks/benchmark_multipoles_covariance_40.dat \
    tests/benchmarks/benchmark_multipoles_covariance_42.dat \
    tests/benchmarks/benchmark_multipoles_covariance_44.dat \
    tests/benchmarks/benchmark_rsd_corrfunc0.dat \
    tests/benchmarks/benchmark_rsd_corrfunc1.dat \
    tests/benchmarks/benchmark_rsd_corrfunc2.dat \
    tests/benchmarks/benchmark_rsd_corrfunc3.dat \
    tests/benchmarks/benchmark_rsd_corrfunc4.dat \
    tests/benchmarks/benchmark_rsd_multipoles0.dat \
    tests/benchmarks/benchmark_rsd_multipoles2.dat \
    tests/benchmarks/benchmark_rsd_multipoles4.dat \
    tests/benchmarks/benchmark_std_corrfunc0.dat \
    tests/benchmarks/benchmark_std_corrfunc1.dat \
    tests/benchmarks/benchmark_std_corrfunc2.dat \
    tests/benchmarks/benchmark_std_corrfunc3.dat \
    tests/benchmarks/benchmark_std_corrfunc4.dat \
    tests/benchmarks/benchmark_std_d1_corrfunc0.dat \
    tests/benchmarks/benchmark_std_d1_corrfunc1.dat \
    tests/benchmarks/benchmark_std_d1_corrfunc2.dat \
    tests/benchmarks/benchmark_std_d1_corrfunc3.dat \
    tests/benchmarks/benchmark_std_d1_corrfunc4.dat \
    tests/benchmarks/benchmark_std_d1_multipoles0.dat \
    tests/benchmarks/benchmark_std_d1_multipoles2.dat \
    tests/benchmarks/benchmark_std_d1_multipoles4.dat \
    tests/benchmarks/benchmark_std_len_corrfunc0.dat \
    tests/benchmarks/benchmark_std_len_corrfunc1.dat \
    tests/benchmarks/benchmark_std_len_corrfunc2.dat \
    tests/benchmarks/benchmark_std_len_corrfunc3.dat \
    tests/benchmarks/benchmark_std_len_corrfunc4.dat \
    tests/benchmarks/benchmark_std_len_multipoles0.dat \
    tests/benchmarks/benchmark_std_len_multipoles2.dat \
    tests/benchmarks/benchmark_std_len_multipoles4.dat \
    tests/benchmarks/benchmark_std_multipoles0.dat \
    tests/benchmarks/benchmark_std_multipoles2.dat \
    tests/benchmarks/benchmark_std_multipoles4.dat \
    README.md \
    manual/coffe_manual.pdf

AM_CPPFLAGS = \
	-DDATADIR=\"$(srcdir)/\" \
	-Wall -Wextra -Wpedantic \
	-I$(srcdir)/src/ -I$(srcdir)/tests/

libcoffe_a_HEADERS = \
    src/common.h \
    src/parser.h \
    src/background.h \
    src/integrals.h \
    src/corrfunc.h \
    src/multipoles.h \
    src/average_multipoles.h \
    src/covariance.h

libcoffe_a_SOURCES = \
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
    src/average_multipoles.c

#include .c and .h in SOURCES so that both appear in dist
coffe_SOURCES = \
    $(libcoffe_a_HEADERS) \
    src/utils.h \
    src/twobessel.h \
    src/twofast.h \
    src/errors.h \
    src/functions.h \
    src/tanhsinh.h \
    src/signal.h \
    src/output.h \
    $(libcoffe_a_SOURCES) \
    src/output.c \
    src/main.c

TESTS = \
    test_background \
    test_integrals \
    test_corrfunc \
    test_multipoles \
    test_covariance

check_PROGRAMS = \
    test_background \
    test_integrals \
    test_corrfunc \
    test_multipoles \
    test_covariance

test_background_SOURCES = \
    tests/tools.h \
    tests/test_background.c \
    src/common.c \
    src/errors.c \
    src/parser.c \
    src/tanhsinh.c \
    src/background.c

test_integrals_SOURCES = \
    tests/tools.h \
    tests/test_integrals.c \
    src/common.c \
    src/errors.c \
    src/parser.c \
    src/tanhsinh.c \
    src/background.c \
    src/twofast.c \
    src/integrals.c

test_corrfunc_SOURCES = \
    tests/tools.h \
    tests/test_corrfunc.c \
    src/common.c \
    src/errors.c \
    src/parser.c \
    src/tanhsinh.c \
    src/background.c \
    src/twofast.c \
    src/integrals.c \
    src/functions.c \
    src/signal.c

test_multipoles_SOURCES = \
    tests/tools.h \
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

test_covariance_SOURCES = \
    tests/tools.h \
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

