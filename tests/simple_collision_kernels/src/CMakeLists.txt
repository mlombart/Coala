cmake_minimum_required(VERSION 3.10)

project(coala)
enable_language (Fortran)

option(COALA_DEBUG "Enable Coala debug features" OFF)


if(COALA_DEBUG)
   message(STATUS "Debug mode")
   set(CMAKE_Fortran_FLAGS "-g -ffree-line-length-none  -fcheck=all -fbacktrace -std=f2008 -fbounds-check -Wall -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan")
else()
   message(STATUS "normal mode")
   set(CMAKE_Fortran_FLAGS "-O3 -ffree-line-length-none -std=f2008 ")
endif()

set(DIR_COALA ../../../src)
include_directories(${DIR_COALA})

set(SOURCES
   ${DIR_COALA}/precision.f90
   ${DIR_COALA}/progress_bar.f90
   ${DIR_COALA}/tension_module.f90
   ${DIR_COALA}/init_massgrid.f90
   ${DIR_COALA}/polynomials_legendre.f90
   ${DIR_COALA}/limiter.f90
   ${DIR_COALA}/reconstruction_g.f90
   ${DIR_COALA}/GQ_legendre_nodes_weights.f90
   ${DIR_COALA}/coag_term/coagflux_functions_GQ.f90
   ${DIR_COALA}/coag_term/coagintflux_functions_GQ.f90
   ${DIR_COALA}/collision_kernels.f90
   ${DIR_COALA}/functions_tabflux_tabintflux.f90
   ${DIR_COALA}/functions_flux_intflux.f90
   ${DIR_COALA}/solver_DG.f90
   ${DIR_COALA}/initial_condition_L2proj_GQ.f90
   ${DIR_COALA}/compute_solver_coag.f90
   ../../utils.f90
   setup.f90
   iterate_coag_k0.f90
   iterate_coag.f90
   coala.f90
)
add_executable(coala ${SOURCES})




