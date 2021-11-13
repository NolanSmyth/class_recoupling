cmake_minimum_required(VERSION 3.4...3.18)

project(ThermalScatteringRates LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(Boost REQUIRED)
find_package(Threads REQUIRED)
find_package(pybind11)
find_package(GSL REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")

set(TSR_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/ThermalScatteringRate.cpp)
set(TSRLIB tsr)

add_library(${TSRLIB} SHARED ${TSR_SOURCES})
target_link_libraries(${TSRLIB} PRIVATE GSL::gsl GSL::gslcblas)
target_include_directories(${TSRLIB}
                           PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

pybind11_add_module(ts_rates
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/ThermalScatteringRatePy.cpp)
target_link_libraries(ts_rates PRIVATE ${TSRLIB})
