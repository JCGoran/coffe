cmake_minimum_required(VERSION 3.10)
project(coffe VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Add the option to disable silent rules
option(DISABLE_SILENT_RULES "Disable silent rules" OFF)
if(NOT DISABLE_SILENT_RULES)
    set(CMAKE_VERBOSE_MAKEFILE OFF)
endif()

if (COFFE_ENABLE_PYTHON)
    find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
    include(cmake/FindCython.cmake)
endif()

# Define source directory
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

# Compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wpedantic")
set(CMAKE_C_STANDARD 99)

# Check for libraries
find_package(GSL 2.1 REQUIRED)
find_library(FFTW3_LIB fftw3 NAMES fftw fftw3 libfftw libfftw3 REQUIRED)
find_library(CONFIG_LIB config NAMES config config++ libconfig libconfig++ REQUIRED)
find_library(M_LIB m NAMES m libm REQUIRED)

# Check for optional libraries
if (COFFE_ENABLE_CUBA)
    find_library(CUBA_LIB cuba NAMES cuba libcuba REQUIRED)
endif()
if (COFFE_ENABLE_CLASS)
    find_library(CLASS_LIB class NAMES class libclass REQUIRED)
endif()

# OpenMP support
find_package(OpenMP)
if(OpenMP_C_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

# Optimization flags
include(CheckCCompilerFlag)
if(COFFE_ENABLE_MATHOPTS)
    if (CMAKE_C_COMPILER_ID STREQUAL "Clang" OR CMAKE_C_COMPILER_ID STREQUAL "GNU")
        set(COFFE_MATHOPTS "-ffast-math")
    else()
        message("Unsupported compiler ${CMAKE_C_COMPILER_ID_STREQUAL}" FATAL_ERROR)
    endif()
else()
    set(COFFE_MATHOPTS "")
endif()
# Include directories
include_directories(${SRC_DIR} ${CMAKE_SOURCE_DIR}/tests)

# Add definitions for conditional compilation
if(CLASS_LIB)
    add_compile_definitions(-DHAVE_CLASS)
endif()
if(CUBA_LIB)
    add_compile_definitions(-DHAVE_CUBA)
endif()

set(COFFE_HEADERS
    ${SRC_DIR}/common.h
    ${SRC_DIR}/covariance.h
    ${SRC_DIR}/utils.h
    ${SRC_DIR}/twobessel.h
    ${SRC_DIR}/errors.h
    ${SRC_DIR}/parser.h
    ${SRC_DIR}/twofast.h
    ${SRC_DIR}/integrals.h
    ${SRC_DIR}/background.h
    ${SRC_DIR}/functions.h
    ${SRC_DIR}/tanhsinh.h
    ${SRC_DIR}/signal.h
    ${SRC_DIR}/corrfunc.h
    ${SRC_DIR}/multipoles.h
    ${SRC_DIR}/average_multipoles.h
)

set(COFFE_LIB_SOURCES
    ${SRC_DIR}/common.c
    ${SRC_DIR}/covariance.c
    ${SRC_DIR}/utils.c
    ${SRC_DIR}/twobessel.c
    ${SRC_DIR}/errors.c
    ${SRC_DIR}/parser.c
    ${SRC_DIR}/twofast.c
    ${SRC_DIR}/integrals.c
    ${SRC_DIR}/background.c
    ${SRC_DIR}/functions.c
    ${SRC_DIR}/tanhsinh.c
    ${SRC_DIR}/signal.c
    ${SRC_DIR}/corrfunc.c
    ${SRC_DIR}/multipoles.c
    ${SRC_DIR}/average_multipoles.c
)

# Main program
add_executable(coffe-cli
    ${SRC_DIR}/main.c
    ${COFFE_LIB_SOURCES}
    ${SRC_DIR}/output.c
)

# Add definition
target_compile_definitions(coffe-cli PRIVATE PACKAGE_STRING="${VERSION}")

# Link libraries
target_link_libraries(coffe-cli ${M_LIB} ${GSL_LIBRARIES} ${FFTW3_LIB} ${CONFIG_LIB} ${CUBA_LIB} ${CLASS_LIB})

target_compile_options(coffe-cli PRIVATE ${COFFE_MATHOPTS})


# Tests
if (COFFE_ENABLE_TESTS)
    enable_testing()

    if (NOT DEFINED COFFE_TEST_DATADIR)
        set(COFFE_TEST_DATADIR ${CMAKE_SOURCE_DIR})
    endif()

    set(COFFE_TESTS background integrals corrfunc multipoles covariance)
    foreach(name IN LISTS COFFE_TESTS)
        # Add the test executable
        add_executable(test_${name}
            ${CMAKE_SOURCE_DIR}/tests/test_${name}.c
            ${COFFE_LIB_SOURCES}
        )
        # Add defines
        target_compile_definitions(test_${name} PRIVATE COFFE_TEST_DATADIR="${COFFE_TEST_DATADIR}")
        # Link libraries
        target_link_libraries(test_${name} ${M_LIB} ${GSL_LIBRARIES} ${FFTW3_LIB} ${CONFIG_LIB} ${CUBA_LIB} ${CLASS_LIB})
        # Register tests
        add_test(NAME test_${name} COMMAND test_${name})
    endforeach()
endif()

if(COFFE_ENABLE_PYTHON)
    add_subdirectory(${CMAKE_SOURCE_DIR}/python/coffe)
endif()

# Installation
install(TARGETS coffe-cli COMPONENT Runtime DESTINATION bin)

# Install headers and other files
install(FILES
    ${COFFE_HEADERS}
    DESTINATION include/coffe
)

install(FILES
    ${SRC_DIR}/WAVENUMBER_HEADER.dat
    ${SRC_DIR}/POWER_SPECTRUM_HEADER.dat
    ${CMAKE_SOURCE_DIR}/settings.cfg
    ${CMAKE_SOURCE_DIR}/separations.dat
    ${CMAKE_SOURCE_DIR}/PkL_CLASS.dat
    ${CMAKE_SOURCE_DIR}/README.md
    DESTINATION share/coffe
)
