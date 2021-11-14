cmake_minimum_required(VERSION 3.14...3.18)

project(Classy)

find_package(Cython REQUIRED)
find_package(NumPy REQUIRED)
find_package(PythonExtensions REQUIRED)

message(STATUS "${Classy_SOURCE_DIR}")

add_cython_target(_classy)
add_library(_classy MODULE ${_classy})
target_link_libraries(_classy ${LIBCLASS} ${TSRLIB} OpenMP::OpenMP_C)
target_include_directories(_classy PUBLIC ${Classy_SOURCE_DIR} ${Classy_SOURCE_DIR}/../include ${NumPy_INCLUDE_DIRS})

python_extension_module(_classy)
install(TARGETS _classy LIBRARY DESTINATION classy)
install(TARGETS ts_rates LIBRARY DESTINATION classy)
