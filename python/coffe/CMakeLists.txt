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

target_link_libraries(coffe PUBLIC ${M_LIB} ${GSL_LIBRARIES} ${FFTW3_LIB} ${CONFIG_LIB} ${CUBA_LIB} ${CLASS_LIB})

install(TARGETS coffe COMPONENT Runtime DESTINATION coffe)
