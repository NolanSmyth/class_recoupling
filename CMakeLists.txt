cmake_minimum_required(VERSION 3.11)
project(CLASS
	VERSION 2.9.4
	DESCRIPTION "CLASS Implementation of Dark Matter Kinetic Recoupling."
	LANGUAGES C)


#============================================================================
#---- Add preprocessor definitions ------------------------------------------
#============================================================================

add_compile_definitions(BBN_DATA_PATH="${CMAKE_SOURCE_DIR}/bbn")

# =================================
# ---- Construct CLASS library ----
# =================================

add_library(class STATIC
	"source/background.c"
	"source/input.c"
	"source/lensing.c"
	"source/nonlinear.c"
	"source/output.c"
	"source/perturbations.c"
	"source/primordial.c"
	"source/spectra.c"
	"source/thermodynamics.c"
	"source/transfer.c"
	"tools/arrays.c"
	"tools/common.c"
	"tools/dei_rkck.c"
	"tools/evolver_ndf15.c"
	"tools/evolver_rkck.c"
	"tools/growTable.c"
	"tools/hyperspherical.c"
	"tools/parser.c"
	"tools/quadrature.c"
	"tools/sparse.c"
	"tools/trigonometric_integrals.c"
	"hyrec/helium.c"
	"hyrec/history.c"
	"hyrec/hydrogen.c"
	"hyrec/hyrec.c"
	"hyrec/hyrectools.c")
target_include_directories(class PUBLIC 
	"include"
	"hyrec")


# ====================================
# ---- Construct Main Entry Point ----
# ====================================
# Make compile the class.c into an executable
add_executable(class_main2 "main/class.c")
# Link to the CLASS library
target_link_libraries(class_main2 class)
# Place the executable in the project-root directory
set_target_properties(class_main2 PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
