#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([coffe], [2.0.4], [goran.jelic-cizmek@unige.ch])
AM_INIT_AUTOMAKE([foreign subdir-objects]) # Does not require NEWS, COPYING, AUTHORS, ChangeLog or README

# silent make https://autotools.io/automake/silent.html
# silent rules enabled by default with 'yes'
# disable silent runles with ./configure --disable-silent-rules
AM_SILENT_RULES([yes]) # less verbose make output

AC_CONFIG_SRCDIR([src/main.c])

AC_LANG([C]) # Use C not C++

: ${CFLAGS=""}

# Checks for programs.
AC_PROG_CC

AC_PROG_CC_C99

AC_PROG_RANLIB

AC_C_CONST

AC_C_INLINE

AC_OPENMP

AC_SUBST(AM_CFLAGS,"$OPENMP_CFLAGS")

AM_CONDITIONAL([FASTFLAG], false)
AX_CHECK_COMPILE_FLAG([-Ofast], [FASTFLAG=true], [FASTFLAG=false])

if $FASTFLAG; then
    AC_SUBST(AM_CFLAGS,"$AM_CFLAGS -Ofast")
else
    AC_SUBST(AM_CFLAGS,"$AM_CFLAGS -O3")
fi

dnl AM_CONDITIONAL([OPENACCFLAG], false)
dnl AX_CHECK_COMPILE_FLAG([-fopenacc], [OPENACCFLAG=true], [OPENACCFLAG=false])

dnl if $OPENACCFLAG; then
dnl    AC_SUBST(AM_CFLAGS,"$AM_CFLAGS -fopenacc")
dnl fi

AC_PROG_INSTALL


# Checks for libraries.

# Found libraries are automatically addded to LIBS
AC_CHECK_LIB([m],[cos],[], AC_MSG_FAILURE([Math library not found.]))
AC_CHECK_LIB([gslcblas],[cblas_dgemm], [], AC_MSG_FAILURE([GSL CBLAS not found.]))
AC_CHECK_LIB([gsl],[gsl_blas_dgemm], [], AC_MSG_FAILURE([GSL not found.]))

AC_SEARCH_LIBS([fftw_plan_dft_r2c_1d], [fftw3],[],[
                 AC_MSG_ERROR([You need to install the FFTW3 library.])
                 ])

AC_SEARCH_LIBS([config_lookup_float], [config],[],[
                 AC_MSG_ERROR([You need to install the libconfig library.])
                 ])


# Checks for optional libraries
AC_SEARCH_LIBS([Cuhre], [cuba],[],[
                 AC_MSG_WARN([CUBA library not found.])
                 ])

AC_SEARCH_LIBS([input_default_params], [class],[],[
                 AC_MSG_WARN([CLASS library not found.])
                 ], ["$OPENMP_CFLAGS"])

dnl AC_SEARCH_LIBS([flint_malloc], [flint], [], [
dnl                 AC_MSG_WARN([FLINT library not found.])
dnl                 ])

dnl AC_SEARCH_LIBS([acb_init], [flint flint-arb arb], [], [
dnl                 AC_MSG_WARN([ARB library not found.])
dnl                 ])

# Checks for header files.
AC_HEADER_ASSERT # ./configure --disable-assert to define NDEBUG
AC_CHECK_HEADERS([stdlib.h stdarg.h complex.h math.h string.h float.h limits.h stddef.h], [], [
    AC_MSG_ERROR([Certain standard headers not found.])
    ])

AC_CHECK_FUNCS([mkdir pow sqrt], [],[
    AC_MSG_ERROR([Certain standard functions not found.])
    ])

AC_ARG_ENABLE(class, [  --enable-class      Automatically includes CLASS library if found], CPPFLAGS="$CPPFLAGS -DHAVE_CLASS", [])
dnl AC_ARG_ENABLE(covariance, [--enable-covariance Automatically includes the covariance calculation if required libraries found], CPPFLAGS="$CPPFLAGS -DHAVE_COVARIANCE", [])
AC_ARG_ENABLE(cuba, [  --enable-cuba    Automatically includes the CUBA library if found], CPPFLAGS="$CPPFLAGS -DHAVE_CUBA", [])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
