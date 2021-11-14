cmake_minimum_required(VERSION 3.14...3.18)

project(Class LANGUAGES C CXX)

find_package(OpenMP)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")

# ===================================================================
# ---- Build Scattering Rate Library --------------------------------
# ===================================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)


find_package(Boost REQUIRED)
find_package(Threads REQUIRED)
find_package(GSL REQUIRED)

find_package(Python COMPONENTS Interpreter Development)
find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)
find_package(pybind11)

set(TSR_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/scattering/src")
set(TSR_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/scattering/include")
set(TSR_SOURCES "${TSR_SOURCE_DIR}/MMS.cpp" "${TSR_SOURCE_DIR}/Tools.cpp")
set(TSRLIB tsr)

add_library(${TSRLIB} STATIC ${TSR_SOURCES})
target_link_libraries(
  ${TSRLIB} PRIVATE
  GSL::gsl
  GSL::gslcblas
  ${Boost_LIBRARIES}
)
target_include_directories(
  ${TSRLIB}
  PUBLIC "${TSR_INCLUDE_DIR}"
)

pybind11_add_module(
  ts_rates
  ${TSR_SOURCE_DIR}/ThermalScatteringRatePy.cpp
)
target_link_libraries(ts_rates PRIVATE ${TSRLIB})
set_target_properties(
  ts_rates
  PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/classy"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/classy"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/classy"
)

# ===================================================================
# ---- Class --------------------------------------------------------
# ===================================================================


set(CLASS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/source")
set(
  CLASS_SOURCES
  "${CLASS_SOURCE_DIR}/background.c"
  "${CLASS_SOURCE_DIR}/input.c"
  "${CLASS_SOURCE_DIR}/lensing.c"
  "${CLASS_SOURCE_DIR}/nonlinear.c"
  "${CLASS_SOURCE_DIR}/output.c"
  "${CLASS_SOURCE_DIR}/perturbations.c"
  "${CLASS_SOURCE_DIR}/primordial.c"
  "${CLASS_SOURCE_DIR}/spectra.c"
  "${CLASS_SOURCE_DIR}/thermodynamics.c"
  "${CLASS_SOURCE_DIR}/transfer.c"
)

set(TOOLS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tools")
set(TOOLS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tools")
set(
  TOOLS_SOURCES
  "${TOOLS_SOURCE_DIR}/arrays.c"
  "${TOOLS_SOURCE_DIR}/common.c"
  "${TOOLS_SOURCE_DIR}/dei_rkck.c"
  "${TOOLS_SOURCE_DIR}/evolver_ndf15.c"
  "${TOOLS_SOURCE_DIR}/evolver_rkck.c"
  "${TOOLS_SOURCE_DIR}/growTable.c"
  "${TOOLS_SOURCE_DIR}/hyperspherical.c"
  "${TOOLS_SOURCE_DIR}/parser.c"
  "${TOOLS_SOURCE_DIR}/quadrature.c"
  "${TOOLS_SOURCE_DIR}/sparse.c"
  "${TOOLS_SOURCE_DIR}/trigonometric_integrals.c"
)

set(HYREC_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/hyrec")
set(HYREC_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/hyrec")
set(
  HYREC_SOURCES
  "${HYREC_SOURCE_DIR}/helium.c"
  "${HYREC_SOURCE_DIR}/history.c"
  "${HYREC_SOURCE_DIR}/hydrogen.c"
  "${HYREC_SOURCE_DIR}/hyrec.c"
  "${HYREC_SOURCE_DIR}/hyrectools.c"
)

# ===================================================================
# ---- Build Class Library ------------------------------------------
# ===================================================================

set(LIBCLASS class)
set(LIBCLASS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(BBN_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bbn")
set(HYREC_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/hyrec")

set(HYREC_ALPHA_INF_DAT "${HYREC_DATA_DIR}/Alpha_inf.dat")
set(HYREC_R_INF_DAT "${HYREC_DATA_DIR}/R_inf.dat")
set(HYREC_TWO_PHOTON_TABLES_DAT "${HYREC_DATA_DIR}/two_photon_tables.dat")
set(BBN_SBBN_2017_DAT "${BBN_DATA_DIR}/sBBN_2017.dat")

configure_file(${LIBCLASS_INCLUDE_DIR}/data.h.in ${LIBCLASS_INCLUDE_DIR}/data.h)

add_library(
  ${LIBCLASS} STATIC
  ${TOOLS_SOURCES}
  ${CLASS_SOURCES}
  ${HYREC_SOURCES}
)
target_include_directories(${LIBCLASS} PUBLIC ${LIBCLASS_INCLUDE_DIR} "${TSR_INCLUDE_DIR}")
target_link_libraries(${LIBCLASS} PUBLIC m ${TSRLIB})

if(OpenMP_C_FOUND)
  target_link_libraries(${LIBCLASS} PUBLIC "OpenMP::OpenMP_C")
endif()

# ===================================================================
# ---- Build Class Executable ---------------------------------------
# ===================================================================

set(CLASSEXE class_main)
set(CLASSEXE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/main")
set(CLASSEXE_SOURCES "${CLASSEXE_SOURCE_DIR}/class.c")
add_executable(${CLASSEXE} ${CLASSEXE_SOURCES})
target_link_libraries(${CLASSEXE} PRIVATE ${LIBCLASS})
target_include_directories(${CLASSEXE} PUBLIC ${LIBCLASS_INCLUDE_DIR})
set_target_properties(
  ${CLASSEXE}
  PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)

# ===================================================================
# ---- Class Python -------------------------------------------------
# ===================================================================

# set(CLASSY classy)
# add_library(${CLASSY} ${})
# add_custom_command(
#   TARGET ${CLASSY} POST_BUILD
#   COMMAND cp python/setup.py python/autosetup.py
#   COMMAND cd python
#   COMMAND export CC=${CC}
#   COMMAND ${Python_Interpreter} autosetup.py install || ${Python_Interpreter} autosetup.py install --user
#   COMMAND rm python/autosetup.py
#   DEPENDS ${LIBCLASS}
#   VERBATIM
# )

# add_custom_target(
#   ReplicatePythonSourceTree ALL ${CMAKE_COMMAND} -P
#   ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ReplicatePythonSourceTree.cmake
#   ${CMAKE_CURRENT_BINARY_DIR}
#   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
# )
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/python)

# cython_add_module(classy ${CMAKE_CURRENT_SOURCE_DIR}/python/cclassy.pyx)


# set(CLASSY classy)
# set(CLASSY_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/python/src")
# set(CLASSY_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/python/include")


# set(CLASSY_SOURCES "${CLASSY_SOURCE_DIR}/Init.cpp" "${CLASSY_SOURCE_DIR}/Tools.cpp")

# pybind11_add_module(${CLASSY} ${CLASSY_SOURCES})
# target_link_libraries(${CLASSY} PUBLIC ${LIBCLASS} GSL)
# target_include_directories(${CLASSY} PRIVATE ${CLASSY_INCLUDE_DIR})

add_subdirectory(classy)

# ===================================================================
# ---- Class Tests --------------------------------------------------
# ===================================================================

set(MAKE_TESTS false)

set(CLASS_TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/test")
set(
  CLASS_TEST_FILES
  "test_loops"
  "test_loops_omp"
  "test_spectra"
  "test_transfer"
  "test_nonlinear"
  "test_perturbations"
  "test_thermodynamics"
  "test_background"
  # "test_scattering_rate"
  "test_hyperspherical"
)

if(${MAKE_TESTS})
  foreach(tfile ${CLASS_TEST_FILES})
    add_executable(${tfile} ${CLASS_TEST_DIR}/${tfile}.c)
    target_link_libraries(${tfile} PUBLIC ${LIBCLASS})
    if(OpenMP_C_FOUND)
      target_link_libraries(${tfile} PUBLIC OpenMP::OpenMP_C)
    endif()
    set_target_properties(
      ${tfile}
      PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )
  endforeach()
endif()
