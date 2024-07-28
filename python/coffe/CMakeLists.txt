include(${CMAKE_SOURCE_DIR}/cmake/UseCython.cmake)
Cython_transpile(${CMAKE_CURRENT_SOURCE_DIR}/coffe.pyx
    OUTPUT_VARIABLE _coffe_sources)

Python_add_library(
  coffe
  MODULE ${_coffe_sources} ${COFFE_LIB_SOURCES}
  WITH_SOABI
)

target_compile_options(
    coffe
    PRIVATE
    ${COFFE_MATHOPTS}
)

# Add all include directories
# Note that we include CLASS directories last due to naming conflicts
# (i.e. both COFFE and CLASS have a `parser.h` file)
target_include_directories(coffe PRIVATE ${CONFIG_INCLUDE_DIRS} ${FFTW3_INCLUDE_DIRS} ${CUBA_INCLUDE_DIRS} ${SRC_DIR} ${CLASS_INCLUDE_DIRS})

# Link the libraries properly
target_link_libraries(coffe PUBLIC ${M_LIB} ${GSL_LIBRARIES} ${FFTW3_LIB} ${FFTW3_LIBRARIES} ${CONFIG_LIB} ${libconfig_LIBRARIES} ${CUBA_LIB} ${CLASS_LIB})

# "install" the coffe Python library
# TODO: if we are not building a wheel, put it somewhere else,
# like `DIR/lib/pythonX.Y/site-packages/coffe`
install(TARGETS coffe COMPONENT Runtime DESTINATION coffe)
